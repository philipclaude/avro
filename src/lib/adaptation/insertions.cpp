//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/adapt.h"
#include "adaptation/filter.h"
#include "adaptation/metric.h"

#include "geometry/entity.h"

#include "mesh/topology.h"

#include "common/tools.h"
#include "avro_types.h"

#include <time.h>
#include <unordered_set>

namespace avro
{

template<typename type>
Insert<type>::Insert( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("inserter");
}

template<typename type>
bool
Insert<type>::visible_geometry( real_t* x , real_t* params , Entity* ep , const std::vector<index_t>& edge ) {

  avro_assert( ep->number() == 2 );
  if (!this->curved_) return true;
  if (ep->interior()) return true;

  this->geometry_cavity_.set_entity(ep);

  // based on the type of face, we may need to flip the sign for the volume calculation
  this->geometry_cavity_.sign() = ep->sign();

  // extract the geometry cavity
  if (edge.size() == 2)
    this->extract_geometry( ep , edge ); // extract shell
  else
    this->extract_geometry( ep ); // no facet information, will assign the cavity to non-ghosts
  avro_assert(this->geometry_topology_.nb() > 0);

  // add the parameter coordinates for the inserted point along with the mapped index
  if (this->topology_.element().parameter()) {
    std::vector<real_t> xu(this->u_.dim());
    this->u_.create( xu.data() );
    this->u_.set_param( this->u_.nb()-1 , params );
  }
  else
    this->u_.create( params );
  this->u2v_.push_back( this->points_.nb()-1 );

  // assign the geometry entities (in case this hasn't already been done)
  for (index_t j = 0; j < this->u_.nb(); j++) {
    if (j < this->u_.nb_ghost()) continue;
    this->u_.set_entity( j , this->points_.entity( this->u2v_[j] ) );
  }

  if (this->topology_.element().parameter()) {
    this->convert_to_parameter(ep);
    this->geometry_cavity_.sign() = 1;
  }

  // the insertion should be visible in parameter space
  index_t ns = this->u_.nb()-1;
  bool accept = this->geometry_cavity_.compute( ns , this->u_[ns] , this->S_ );
  if (!accept) {
    avro_assert( this->topology_.element().parameter());
    return false;
  }
  avro_assert(accept);

  if (this->topology_.element().parameter()) {
    this->convert_to_physical();
  }

  // reset the geometry checker to this face (which also computes the normals at the vertices of the geometry topology)
  this->geometry_inspector_.reset(ep);
  int s = this->geometry_inspector_.signof( this->geometry_cavity_ );
  if (s < 0) {
    return false;
  }

  avro_assert( this->geometry_inspector_.positive_volumes(this->geometry_cavity_,this->geometry_cavity_.sign()) );

  return true;
}


template<typename type>
bool
Insert<type>::apply( const index_t e0 , const index_t e1 , real_t* x , real_t* u , const std::vector<index_t>& shell , int idx ) {

  elems_.clear();

  // insertions should not enlarge the initial set of cavity elements
  this->enlarge_ = false;
  this->check_visibility_ = true;
  enlarged_ = false; // no enlargement (for geometry) yet

  // determine the cavity elements if none were provided
  if (shell.size() == 0) {
    // likely an interior edge, cavity is the shell around the edge
    this->topology_.intersect( {e0,e1} , elems_ );
  }
  else {
    // likely a geometry insertion, use the provided shell
    elems_ = shell;
  }
  if (elems_.size() == 0) return false;

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
  if (entity0 == NULL || entity1 == NULL) {
    // interior split
    entitys = NULL;
    bodys   = 0;
  }
  else {
    // set the geometry of the inserted vertex as the intersection
    // of the endpoint geometries
    entitys = this->geometry(e0,e1);

    // determine the body
    if (entitys == NULL) {
      // interior split
      bodys = 0;
    }
    else {
      int body0 = this->topology_.points().body(e0);
      int body1 = this->topology_.points().body(e1);

      if (body0 == body1) bodys = body0;
      else bodys = 0; // not sure how this can even happen but just being safe
    }
  }

  #if 1 // philip may 28 (moved up from below)
  this->set_entity( entitys );
  this->topology_.points().set_entity( ns , entitys );
  this->topology_.points().set_param( ns , u );
  this->topology_.points().body(ns) = bodys;
  #endif

  if (entitys != NULL && !this->topology_.element().parameter()) {
    // enlarge the cavity for boundary insertions
    // if the shell was given, then this should have been precomputed
    std::vector<index_t> elems0 = elems_;
    if (!this->find_geometry( x , elems_ )) {
      // we could not enlarge to find the inserted point
      this->topology_.remove_point(ns);
      return false;
    }
    if (elems_.size() != elems0.size()) {
      enlarged_ = true;
    }
  }

  bool accept;
  if (this->topology_.element().parameter()) {

    avro_assert( entitys != nullptr );
    avro_assert( elems_.size() == 2 );

    if (entitys->number() == 2) {
      // check if the insertion is visible to the boundary of the geometry
      // cavity in the parameter space of the geometry entity
      this->Cavity<type>::clear();
      for (index_t k = 0; k < elems_.size(); k++)
        this->add_cavity(elems_[k]);
      accept = visible_geometry( x , u , entitys );
      if (!accept) {
        this->topology_.remove_point(ns);
        return false;
      }

      this->check_visibility(false);
      accept = this->compute( ns , x , elems_ );
      avro_assert(accept);
      this->check_visibility(true);
    }
    else {

      // loop through the parents of this geometry Edge
      for (index_t k = 0; k < entitys->nb_parents(); k++) {

        if (entitys->parents(k)->number() != 2 || !entitys->parents(k)->tessellatable())
          continue;

        // determine all elements on this parent Face
        this->Cavity<type>::clear();
        for (index_t j = 0; j < elems_.size(); j++) {
          Entity* entityp = BoundaryUtils::geometryFacet( this->topology_.points() , this->topology_(elems_[j]) , this->topology_.nv(elems_[j]) );
          avro_assert( entityp != nullptr );
          avro_assert( entityp->number() == 2 );
          if (entityp != entitys->parents(k)) continue;
          this->add_cavity( elems_[j] );
        }

        Entity* parent = entitys->parents(k);
        if (!visible_geometry( x , u , parent , {e0,e1} )) {
          this->topology_.remove_point(ns);
          return false;
        }
      }

      // we need to insert the correct topology
      // so recompute the cavity (visibility check off) for all
      // original cavity elements
      this->check_visibility(false);
      bool accept = this->compute( ns , x , elems_ );
      avro_assert(accept);
      this->check_visibility(true);
    }

    return true;

  } // if parameter space

  // attempt the operator
  accept = this->compute( ns , x , elems_ );
  if (!accept) {
    // the cavity requested enlargment which is possible because of a
    // minimum volume constraint in the cavity operator
    // remove the vertex and inform the caller the insertion is not allowed
    this->topology_.remove_point(ns);
    return false;
  }

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positive_implied_metrics()) {
    this->topology_.remove_point(ns);
    return false;
  }

  if (entitys != NULL && entitys->number() == 2 && this->topology_.number() > 2) {
    // check if the insertion is visible to the boundary of the geometry
    // cavity in the parameter space of the geometry entity
    // TODO check visibility on all faces of this is an edge entity
    avro_assert( !this->topology_.element().parameter() );
    accept = visible_geometry( x , u , entitys );
    if (!accept) {
      this->topology_.remove_point(ns);
      return false;
    }
  }

  if (entitys != NULL && this->topology_.number() == 3 && !entitys->interior()) {
    index_t nb_ghost = 0;
    for (index_t k = 0; k < this->nb(); k++) {
      if (this->ghost(k))
        nb_ghost++;
    }
    if (nb_ghost != 4) {
      this->topology_.remove_point(ns);
      return false;
    }
  }

  if (this->nb_removed_nodes() > 0) {
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

#if AVRO_MPI
static index_t MAX_VALENCY[5] = { 1 , 3 , 20 , 50 , 350 };
#endif

template<typename type>
void
AdaptThread<type>::split_edges( real_t lt, bool limitlength , bool swapout ) {

  avro_assert( metric_.check(topology_) );

  real_t dof_factor = params_["insertion volume factor"];

  index_t nb_inserted,nb_inserted_total = 0;
  index_t nb_length_rejected = 0;
  index_t nb_quality_rejected = 0;
  index_t nb_visiblity_rejected = 0;
  index_t nb_bad_interp = 0;
  index_t nb_high_valency = 0;
  real_t lmin1,lmax1;

  inserter_.nb_parameter_rejections() = 0;
  inserter_.nb_parameter_tests() = 0;

  // don't be too restrictive with insertions when the quality is already good
  real_t Q0 = worst_quality(topology_,metric_);
  if (Q0 > 0.4)  Q0 = 0.1;
  if (Q0 < 1e-3) Q0 = 1e-3;

  index_t pass = 0;
  bool done = false;
  bool problem = false;
  printf("-> performing edge splits on edges with lt > %1.3f and dof_factor %g:\n",lt,dof_factor);
  while (!done && !problem) {

    // anything beyond 20 passes is a bit much,
    // the metric is probably way off from the current mesh
    if (pass > 20) {
      printf("warning: too many insertions, metric too far from current mesh.\n");
      break;
    }

    nb_inserted = 0;
    nb_length_rejected = 0;
    nb_quality_rejected = 0;
    nb_visiblity_rejected = 0;
    nb_bad_interp = 0;
    nb_high_valency = 0;

    clock_t PASS_T0 = clock();

    done = true;

    topology_.inverse().build();
    topology_.inverse().use_ball(true);

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

    // insert the points
    std::unordered_set<index_t> flagged;
    std::vector<index_t> shell;
    std::vector<index_t> N;
    for (index_t k = 0; k < filter.nb_candidates(); k++) {

      // index of the vertex stored in the filter
      index_t idx = filter.candidate(k);

      // get the edge along which this insertion takes place
      index_t n0 = filter.edge( idx , 0 );
      index_t n1 = filter.edge( idx , 1 );

      // insertions on the edges with fixed nodes are not allowed
      // as these are partition boundaries
      if (topology_.points().fixed(n0) && topology_.points().fixed(n1)) {
        topology_.points().age( n0 )++;
        topology_.points().age( n1 )++;
        continue;
      }

      // do not insert on ghost edges
      if (n0 < topology_.points().nb_ghost() ||
          n1 < topology_.points().nb_ghost())
        continue;

      // do not insert on edges with vertices that have been flagged
      if (flagged.find(n0) != flagged.end() || flagged.find(n1) != flagged.end())
        continue;

      // the metric needs to be interpolated for the filter to evaluate lengths
      // so we need to add the vertex and also add the interpolated metric
      index_t ns = topology_.points().nb();
      topology_.points().create(filter[idx]);
      topology_.points().set_entity(ns, inserter_.geometry(n0,n1) );
      topology_.points().set_param( ns , filter.u(idx) );
      bool interp_result = metric_.add(n0,n1,ns,filter[idx]);
      if (!interp_result) {
        // something bad happened in the metric interpolation...
        // probably because of a curved geometry
        topology_.points().remove(ns);
        nb_bad_interp++;
        continue;
      }

      // notify the inverse topology that we want to store extra data
      topology_.inverse().create(1);

      // check the distance to nearby points
      shell.clear();
      topology_.intersect( {n0,n1} , shell );
      N.clear();
      for (index_t j = 0; j < shell.size(); j++)
      {
        for (index_t i = 0; i < topology_.nv(shell[j]); i++)
          N.push_back( topology_(shell[j],i) );
      }
      uniquify(N);
      bool too_short = false;

      // set the minimum length any insertion can create
      real_t Lmin = sqrt(0.5);

      // if the current length is greater than 4.0, we need to be more flexible
      //if (lk>4.0) Lmin = 0.0; // philip: commented june 15 2021

      // also relax the insertion criterion when we insert on geometry Edges
      Entity* ge = inserter_.geometry(n0,n1);
      if (ge != NULL && ge->number() == 1)
        Lmin = 0.25;

      // set the geometry entity for the inserted point
      topology_.points().set_entity(ns,ge);

      for (index_t j = 0; j < N.size(); j++) {
        real_t lj = metric_.length( topology_.points() , N[j] , ns );
        if (lj < Lmin) {
          too_short = true;
          break;
        }
      }

      if (too_short && limitlength) {
        // remove the interpolated tensor with its associated vertex
        topology_.points().remove(ns);
        metric_.remove(ns);
        topology_.inverse().remove(ns);
        nb_length_rejected++;
        continue;
      }
      avro_assert( metric_.check(topology_) );

      // the vertex needs to be removed because the inserter will add it again
      // note: this might seem inefficient because we are calling std::erase for every candidate edge
      // but I (philip) profiled this and it's not that bad - plus, std::erase isn't too bad when deleting
      // things from the back of a list (which we are doing here)
      topology_.points().remove(ns);

      // apply the insertion
      inserter_.clear();
      inserter_.restart();
      inserter_.delay() = true;
      bool result = inserter_.apply( n0 , n1 , filter[idx] , filter.u(idx) , shell );
      if (!result) {
        // the metric needs to be removed because the vertex was rejected.
        // for curverd geometries, the cavity might need to be enlarged.
        // the vertex was already removed (above)
        // so only the metric requires removing
        metric_.remove(ns);
        nb_visiblity_rejected++;
        continue;
      }

      // check that quality didn't degrade below the requested limit
      inserter_.cavity_quality().resize( inserter_.nb() );
      real_t qwi = worst_quality(inserter_,metric_ , inserter_.cavity_quality().data() );
      if (qwi < Q0) {
        topology_.points().remove(ns);
        metric_.remove(ns);
        topology_.inverse().remove(ns);
        nb_quality_rejected++;
        continue;
      }

      // we always need to maintain a closed boundary
      if (!inserter_.closed_boundary()) {
        // it's going to be really hard to come back from something like this
        //avro_assert_not_reached; // is this even possible?
        topology_.points().remove(ns);
        metric_.remove(ns);
        topology_.inverse().remove(ns);
        problem = true;
        continue;
      }

      // check that we do not create valencies which are too large
      #if AVRO_MPI
      if (mpi::size() > 1) {
        bool high_valency = false;
        std::vector<index_t> surrounding;
        for (index_t j = 0; j < N.size(); j++) {
          if (N[j] < topology_.points().nb_ghost()) continue;
          surrounding.clear();
          topology_.inverse().ball( N[j] , surrounding );
          if (surrounding.size() > MAX_VALENCY[topology_.number()]) {
            high_valency = true;
            break;
          }
        }
        if (high_valency) {
          topology_.points().remove(ns);
          metric_.remove(ns);
          topology_.inverse().remove(ns);
          nb_high_valency++;
          continue;
        }
      }
      #endif

      // apply the insertion into the topology
      topology_.apply(inserter_);
      avro_assert( metric_.check(topology_) );

      // set the age of the endpoint vertices to 0 since a split was performed
      // this information is needed to penalize interfaces which border elements
      // that want to be adapted (when migrating interfaces in parallel)
      topology_.points().set_age( n0 , 0 );
      topology_.points().set_age( n1 , 0 );

      // ensure no vertices were removed (even if the cavity was enlarged)
      avro_assert( inserter_.nb_removed_nodes() == 0);

      // check if the cavity was enlarged to turn off some of the existing edges
      if (inserter_.enlarged()) {
        // turn off points in the cavity, we can attempt them on the next pass
        const std::vector<index_t>& nodes = inserter_.nodes();
        for (index_t j = 0; j < nodes.size(); j++)
          flagged.insert( nodes[j] );
      }

      // the filter is happy
      filter.accept( idx , ns );
      done = false; // since an insertion was performed

      nb_inserted++;
    }

    printf("\t\tinserted %lu, lrej = %lu, qrej = %lu, vrej = %lu, prej = %lu/%lu, val_rej = %lu, interp_rej = %lu\n",
                nb_inserted,nb_length_rejected,nb_quality_rejected,nb_visiblity_rejected,
                inserter_.nb_parameter_rejections(),inserter_.nb_parameter_tests(),nb_high_valency,nb_bad_interp);
    nb_inserted_total += nb_inserted;
    printf("\t\tpass time = %g\n",real_t(clock()-PASS_T0)/real_t(CLOCKS_PER_SEC));
    pass++;
  }

  topology_.inverse().use_ball(false);


  // analyze the resulting edge lengths
  std::vector<index_t> edges;
  topology_.get_edges(edges);
  std::vector<real_t> lengths(edges.size()/2);
  for (index_t k = 0; k < edges.size()/2; k++)
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

} // avro
