#include "adaptation/adapt.h"
#include "adaptation/filter.h"
#include "adaptation/metric.h"
#include "adaptation/parameters.h"

#include "mesh/topology.h"

#include "common/tools.h"
#include "common/types.h"

#include <unordered_set>

namespace luma
{

template<typename type>
Insert<type>::Insert( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("inserter");
}

template<typename type>
bool
Insert<type>::visibleParameterSpace( real_t* x , real_t* params , Entity* ep0 )
{
  EGADS::Object* ep = (EGADS::Object*) ep0;

  luma_assert( ep->number()==2 );
  if (!this->curved_) return true;
  if (ep->interior()) return true;

  // based on the type of face, we may need to flip the sign for the volume calculation
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*ep->object(), &oclass, &mtype,&ref, &prev, &next) );
  if (mtype==SREVERSE)
    this->gcavity_.sign() = -1;
  else
    this->gcavity_.sign() = 1;

  // extract the geometry cavity
  this->extractGeometry( ep ); // no facet information, will assign the cavity to non-ghosts
  luma_assert(this->G_.nb()>0);

  // add the parameter coordinates for the inserted point along with the mapped index
  this->u_.create( params );
  this->u2v_.push_back( this->points_.nb()-1 );

  // the insertion should be visible in parameter space
  bool accept = this->gcavity_.compute( this->u_.nb()-1 , params , this->S_ );
  luma_assert(accept);

  // check the orientation of the facets
  GeometryOrientationChecker checker(this->points_,this->u_,this->u2v_,ep);
  int s = checker.signof( this->gcavity_ );
  if (s<0) return false;

  luma_assert( checker.hasPositiveVolumes(this->gcavity_,mtype) );

  return true;
}


template<typename type>
bool
Insert<type>::apply( const index_t e0 , const index_t e1 , real_t* x , real_t* u , const std::vector<index_t>& shell )
{
  //std::vector<index_t> elems;
  elems_.clear();

  // insertions should not enlarge the initial set of cavity elements
  this->enlarge_ = false;
  this->checkVisibility_ = true;
  enlarged_ = false; // no enlargement (for geometry) yet

  // determine the cavity elements if none were provided
  if (shell.size()==0)
  {
    // likely an interior edge, cavity is the shell around the edge
    this->topology_.intersect( {e0,e1} , elems_ );
  }
  else
  {
    // likely a geometry insertion, use the provided shell
    elems_ = shell;
  }

  // index of the vertex to be inserted
  index_t ns = this->topology_.points().nb();

  // we need to add the vertex
  this->topology_.points().create(x);

  // this might be a boundary insertion,
  // compute intersection of geometry entities of each vertex
  Entity* entity0 = this->topology_.points().entity(e0);
  Entity* entity1 = this->topology_.points().entity(e1);
  Entity* entitys;
  int bodys;
  if (entity0==NULL || entity1==NULL)
  {
    // interior split
    entitys = NULL;
    bodys   = 0;
  }
  else
  {
    // set the geometry of the inserted vertex as the intersection
    // of the endpoint geometries
    entitys = this->geometry(e0,e1);

    // determine the body
    if (entitys==NULL)
    {
      // interior split
      bodys = 0;
    }
    else
    {
      int body0 = this->topology_.points().body(e0);
      int body1 = this->topology_.points().body(e1);

      if (body0==body1) bodys = body0;
      else bodys = 0; // not sure how this can even happen but just being safe
    }
  }

  if (entitys!=NULL)
  {
    // enlarge the cavity for boundary insertions
    // if the shell was given, then this should have been precomputed
    std::vector<index_t> elems0 = elems_;
    if (!this->findGeometry( x , elems_ ))
    {
      // we could not enlarge to find the inserted point
      this->topology_.remove_point(ns);
      return false;
    }
    if (elems_.size()!=elems0.size())
    {
      enlarged_ = true;
    }
  }

  // attempt the operator
  bool accept = this->compute( ns , x , elems_ );
  if (!accept)
  {
    // the cavity requested enlargment which is possible because of a
    // minimum volume constraint in the cavity operator
    // remove the vertex and inform the caller the insertion is not allowed
    this->topology_.remove_point(ns);
    return false;
  }

  // check if all produce elements have a positive determinant of implied metric
  #if 0
  // pcaplan REMOVE THIS in master branch (this ruins timing and should only be used for curved=true)
  if (!this->positiveImpliedMetrics())
  {
    this->topology_.remove_point(ns);
    return false;
  }
  #endif


  if (entitys!=NULL && entitys->number()==2)
  {
    // check if the insertion is visible to the boundary of the geometry
    // cavity in the parameter space of the geometry entity
    accept = visibleParameterSpace( x , u , entitys );
    if (!accept)
    {
      this->topology_.remove_point(ns);
      return false;
    }
  }

  if (entitys!=NULL && this->topology_.number()==3 && !entitys->interior())
  {
    index_t nb_ghost = 0;
    for (index_t k=0;k<this->nb();k++)
    {
      if (this->ghost(k))
        nb_ghost++;
    }
    PRIMITIVE_CHECK( nb_ghost==4 );
    if (nb_ghost!=4)
    {
      this->topology_.remove_point(ns);
      return false;
    }
  }

  if (this->nb_removedNodes()>0)
  {
    // when we are enlarging, do not allow points to be removed
    this->topology_.remove_point(ns);
    return false;
  }

  // set the geometry entity and body
  this->topology_.points().set_entity( ns , entitys );
  this->topology_.points().set_param( ns , u );
  this->topology_.points().body(ns) = bodys;

  // apply the operator to the topology if the caller did not request a delay
  if (!this->delay_)
    this->topology_.apply(*this);
  return true;
}

template<typename type>
void
AdaptThread<type>::split_edges( real_t lt, bool limitlength , bool swapout )
{
  luma_assert( metric_.check(topology_) );

  real_t dof_factor = params_.insertion_volume_factor();

  index_t nb_swaps;
  index_t nb_inserted,nb_inserted_total = 0;
  index_t nb_length_rejected = 0;
  index_t nb_quality_rejected = 0;
  index_t nb_visiblity_rejected = 0;
  index_t nb_count_rejected = 0;
  real_t lmin1,lmax1;

  inserter_.nb_parameter_rejections() = 0;
  inserter_.nb_parameter_tests() = 0;

  //topology_.evaluate(metric);
  real_t Q0 = worst_quality(topology_,metric_);

  // don't be too restrictive with insertions when the quality is already good
  if (Q0>0.4) Q0 = 0.1;

  std::vector<real_t> lens;
  metric_.lengths(topology_,lens);
  real_t lmin0 = *std::min_element( lens.begin() , lens.end() );
  real_t lmax0 = *std::max_element( lens.begin() , lens.end() );

  index_t pass = 0;
  bool done = false;
  printf("-> performing edge splits on edges with lt > %1.3f and dof_factor %g:\n",lt,dof_factor);
  while (!done)
  {
    // anything beyond 20 passes is borderline ridiculous
    // the metric is probably way off from the current mesh
    if (pass>20)
    {
      printf("warning: too many insertions, metric too far from current mesh.\n");
      break;
    }

    nb_inserted = 0;
    nb_swaps = 0;
    nb_length_rejected = 0;
    nb_quality_rejected = 0;
    nb_visiblity_rejected = 0;
    nb_count_rejected = 0;

    done = true;

    // setup the insertion filter
    Filter filter( topology_.points().dim() );

    // add the current points in the topology
    for (index_t k=0;k<topology_.points().nb();k++)
      filter.createPermanent( topology_.points()[k] );

    // create the candidates on edges longer than lt in the target space
    filter.generateCandidates( topology_ , metric_ , lt , inserter_ );

    printf("\t pass %lu: long = %lu, l = [%3.4f,%3.4f] -> insert %lu\n",
                pass,filter.nb_long(),filter.minlength(),filter.maxlength(),
                filter.nb_candidates());

    // build the kd-tree search structure only if we need to check the
    // distance to the nearest neighbours
    filter.buildTree();

    // insert the points
    std::unordered_set<index_t> removed;
    std::unordered_set<index_t> flagged;
    std::vector<index_t> shell;
    std::vector<index_t> N;
    for (index_t k=0;k<filter.nb_candidates();k++)
    {
      // index of the vertex stored in the filter
      index_t idx = filter.candidate(k);

      // get the edge along which this insertion takes place
      index_t n0 = filter.edge( idx , 0 );
      index_t n1 = filter.edge( idx , 1 );

      real_t lk = metric_.length( topology_.points() , n0 , n1 );

      // insertions on the edges with fixed nodes are not allowed
      // as these are partition boundaries
      if (topology_.points().fixed(n0) && topology_.points().fixed(n1))
        continue;

      // do not insert on ghost edges
      if (n0<topology_.points().nb_ghost() ||
          n1<topology_.points().nb_ghost())
        continue;

      // check if the points were removed or flagged
      if (removed.find(n0)!=removed.end() || removed.find(n1)!=removed.end())
        continue;

      if (flagged.find(n0)!=flagged.end() || flagged.find(n1)!=flagged.end())
        continue;

      // the metric needs to be interpolated for the filter to evaluate lengths
      // so we need to add the vertex and also add the interpolated metric
      index_t ns = topology_.points().nb();
      topology_.points().create(filter[idx]);
      metric_.add(n0,n1,filter[idx]);

      // notify the inverse topology that we want to store extra data
      topology_.inverse().create(1);

      // check the distance to nearby points
      shell.clear();
      topology_.intersect( {n0,n1} , shell );
      N.clear();
      for (index_t j=0;j<shell.size();j++)
      {
        for (index_t i=0;i<topology_.nv(shell[j]);i++)
          N.push_back( topology_(shell[j],i) );
      }
      uniquify(N);
      bool bad = false;

      // set the minimum length any insertion can create
      real_t Lmin = sqrt(0.5);

      // if the current length is greater than 4.0, we need to be more flexible
      if (lk>4.0) Lmin = 0.0;

      // also relax the insertion criterion when we insert on geometry Edges
      Entity* ge = inserter_.geometry(n0,n1);
      if (ge!=NULL && ge->number()==1)
        Lmin = 0.25;

      for (index_t j=0;j<N.size();j++)
      {
        real_t lj = metric_.length( topology_.points() , N[j] , ns );
        if (lj<Lmin)
        {
          bad = true;
          break;
        }
      }

      bool swapped;
      if (bad && limitlength)
      {
        // remove the interpolated tensor with its associated vertex
        topology_.points().remove(ns);
        metric_.remove(ns);
        topology_.inverse().remove(ns);
        nb_length_rejected++;

        // option to swap out of the rejected configuration
        if (swapout)
        {
          swapped = swap_edge(n0,n1,Q0,lmin0,lmax0);
          if (swapped) nb_swaps++;
        }
        continue;
      }
      luma_assert( metric_.check(topology_) );

      // the vertex needs to be removed because the inserter will add it again
      // TODO clean up this inefficiency
      topology_.points().remove(ns);

      // apply the insertion
      inserter_.clear();
      inserter_.restart();
      inserter_.delay() = true;
      bool result = inserter_.apply( n0 , n1 , filter[idx] , filter.u(idx) );
      if (!result)
      {
        // the metric needs to be removed because the vertex was rejected.
        // for curverd geometries, the cavity might need to be enlarged.
        // the vertex was already removed (above)
        // so only the metric requires removing
        metric_.remove(ns);
        nb_visiblity_rejected++;
        continue;
      }

      // if the inserter was enlarged, don't be too restrictive with quality
      //inserter_.evaluate( metric );
      real_t qwi = worst_quality(inserter_,metric_);
      if (qwi<Q0)
      {
        topology_.points().remove(ns);
        metric_.remove(ns);
        topology_.inverse().remove(ns);
        nb_quality_rejected++;
        continue;
      }

      // check if the metric volume is respected
      real_t vol = 0.0;
      for (index_t k=0;k<inserter_.nb();k++)
        vol += metric_.volume( inserter_ , k );
      index_t count = index_t(vol/topology_.master().reference().vunit());
      if (dof_factor>0 && dof_factor*count<inserter_.nb_real())
      {
        if (ge==NULL || ge->number()>2)
        {
          // it's going to be really hard to come back from something like this
          topology_.points().remove(ns);
          metric_.remove(ns);
          topology_.inverse().remove(ns);
          nb_count_rejected++;
          continue;
        }
      }

      // apply the insertion into the topology
      topology_.apply(inserter_);
      luma_assert( metric_.check(topology_) );

      if (swapout)
      {
        // if we're doing swaps then we need to update the quality information
        //inserter_.evaluate(metric);
        //topology_.update( inserter_ , metric );
      }

      // determine if any points were removed
      for (index_t j=0;j<inserter_.nb_removedNodes();j++)
      {
        index_t removed_node = inserter_.removedNode(j);
        if (removed.find(removed_node)==removed.end())
        {
          printf("vertex %lu was removed!\n",removed_node);
          removed.insert( removed_node );
        }
      }

      // check if the cavity was enlarged to turn off some of the existing edges
      if (inserter_.enlarged())
      {
        // turn off points in the cavity, we can attempt them on the next pass
        const std::vector<index_t>& nodes = inserter_.nodes();
        for (index_t j=0;j<nodes.size();j++)
          flagged.insert( nodes[j] );
      }

      // the filter is happy
      filter.accept( idx , ns );
      done = false; // since an insertion was performed

      nb_inserted++;
    }

    // delete the points that were removed during cavity enlargements
    std::vector<index_t> removed_indices;
    for (std::unordered_set<index_t>::iterator it=removed.begin();it!=removed.end();it++)
      removed_indices.push_back( *it );

    std::sort( removed_indices.begin() , removed_indices.end() );
    std::reverse( removed_indices.begin() , removed_indices.end() );
    for (index_t j=0;j<removed_indices.size();j++)
    {
      topology_.remove_point( j );
      metric_.remove(j);
      topology_.inverse().remove(j);
    }

    printf("\t\tinserted %lu, swapped %lu, lrej = %lu, qrej = %lu, vrej = %lu, prej = %lu/%lu, dof_rej = %lu\n",
                nb_inserted,nb_swaps,nb_length_rejected,nb_quality_rejected,nb_visiblity_rejected,
                inserter_.nb_parameter_rejections(),inserter_.nb_parameter_tests(),nb_count_rejected);
    nb_inserted_total += nb_inserted;

    pass++;
  }

  // analyze the resulting edge lengths
  std::vector<index_t> edges;
  topology_.get_edges(edges);
  std::vector<real_t> lengths(edges.size()/2);
  for (index_t k=0;k<edges.size()/2;k++)
    lengths[k] = metric_.length( topology_.points() , edges[2*k] , edges[2*k+1] );
  lmin1 = *std::min_element( lengths.begin() , lengths.end() );
  lmax1 = *std::max_element( lengths.begin() , lengths.end() );
  printf("done insertions: total insert = %lu, l = [%3.4g,%3.4g]\n",
    nb_inserted_total,lmin1,lmax1);
  printf("\tnb_parameter_rejections = (%lu/%lu)\n",
            inserter_.nb_parameter_rejections(),inserter_.nb_parameter_tests());

}

template class Insert<Simplex>;
template class AdaptThread<Simplex>;

} // luma
