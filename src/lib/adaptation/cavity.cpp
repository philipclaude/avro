#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "adaptation/cavity.h"
#include "adaptation/geometry.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/geometry.h"

#include <egads.h>

#include <algorithm>
#include <cstdio>
#include <array>
#include <cmath>
#include <limits>

namespace avro
{

template<typename type>
Cavity<type>::Cavity( Topology<type>& _topology ) :
  Topology<type>(_topology.points(),_topology.number()),
  topology_(_topology),
  node_removal_allowed_(true),
  enlarge_(true),
  check_visibility_(true),
  star_(0),
  boundary_( topology_.points() , topology_.number()-1 ),
  minvol_(1e-12),
  sign_(1.0),
  nb_error_(0),
  rethrow_(false),
  mk_(topology_.master()),
  ignore_(false),
  entity_(nullptr)
{
  // allocate enough space for the re-inserted point coordinates
  //point_.resize( topology_.points().dim() );
  point_.resize( topology_.number() );

  // we do not want to store the elements because they will be oriented
  // to maintain a positive volume
  boundary_.set_sorted(false);
  this->set_sorted(false);

  this->master().set_parameter( topology_.master().parameter() );
}

template<typename type>
bool
Cavity<type>::positive_implied_metrics()
{
  for (index_t k=0;k<this->nb();k++)
  {
    // retrieve the implied metric of the element
    if (mk_.determinant(this->points_,this->operator()(k),this->nv(k))<0.0)
      return false;
  }
  return true;
}

template<typename type>
bool
Cavity<type>::compute( const index_t p , real_t* x , const std::vector<index_t>& C0 )
{
  cavity_.clear();
  boundary_.clear();
  nodes_.clear();
  removed_nodes_.clear();

  // assign the coordinates and index of the star
  for (coord_t d=0;d<point_.size();d++)
    point_[d] = x[d];
  star_ = p;

  // compute the cavity about p with elements initialized to C0
  for (index_t k=0;k<C0.size();k++)
    add_cavity(C0[k]);

  // enlarge the cavity until the point is visible from every face of the boundary
  if (check_visibility_)
  {
    bool accept = enlarge();
    if (!accept)
    {
      // cavity was rejected, point cannot be visible
      return false;
    }
  }
  else
  {
    avro_assert(C0.size()>0);
  }

  // accumulate the nodes of the cavity elements
  // this is needed to detect if nodes are removed
  compute_nodes();

  // compute the boundary of the cavity elements
  try
  {
    if (!compute_boundary())
    {
      avro_assert_not_reached;
      return false;
    }
  }
  catch(...)
  {
    nb_error_++;
    if (rethrow_) // option for debugging
      avro_assert_not_reached; // any kind of failed assertion
    return false;
  }

  // sort the elements in the cavity,
  // necessary for the application of the operator by the mesh topology
  std::sort( cavity_.begin() , cavity_.end() );

  // apply the cavity to the candidate (also computes removed nodes)
  apply();

  if (!node_removal_allowed_ && nb_removed_nodes()>0)
  {
    printf("node removal not allowed\n");
    return false;
  }

  return true;
}

template<typename type>
bool
Cavity<type>::compute_boundary()
{
  avro_assert( boundary_.nb()==0 );
  avro_assert_msg( topology_.closed() , "requires implementation + testing without closed topologies" );

  idx_.clear();

  const index_t nf = topology_.number()+1;
  std::vector<index_t> b(nf-1);

  // loop through all cavity elements
  for (index_t k=0;k<nb_cavity();k++)
  {
    for (index_t j=0;j<nf;j++)
    {
      int neighbour = topology_.neighbours()( cavity_[k] , j );

      // for closed topologies there should be no boundary
      if (topology_.closed())
      {
        if (neighbour<0)
        {
          printf("topology is closed but there is a boundary neighbour!\n");
          printf("element %lu, facet %lu has empty neighbour\n",cavity_[k],j);
          return false;
        }
      }

      if (topology_.closed())
        avro_assert( neighbour>=0 );

      // skip neighbours that are already in the cavity
      if (contains(index_t(neighbour)) && topology_.closed())
      {
        continue;
      }

      if (!topology_.closed() && neighbour>=0 && contains(index_t(neighbour)))
        continue;

      // neighbour is not in cavity which means we hit a boundary facet
      // set the indices
      index_t count = 0;
      for (index_t i=0;i<nf;i++)
      {
        if (i==j)
        {
          // index of the re-inserted point relative to the boundary facet
          idx_.push_back(i);
          continue;
        }
        b[count++] = topology_( cavity_[k] , i );
      }

      // store the boundary facet
      boundary_.add( b.data() , b.size() );
    }
  }
  return true;
}


template<typename type>
void
Cavity<type>::compute_geometry( Entity* entity0 , Topology<type>& geometry , std::map<index_t,index_t>& v2u , std::vector<index_t>& u2v )
{
  // this function should be called after the cavity operator has been applied
  // because we need to find ghost cavity elements to compute the boundary.
  // first check if we are computing coordinates in xyz- or uv-space
  bool uvcoords = false;
  EGADS::Object* entity = (EGADS::Object*) entity0;
  if (geometry.points().dim()==topology_.points().dim()-1)
  {
    // the points of the geometry topology are requested to be in uv-space!
    uvcoords = true;
  }
  avro_assert_msg( topology_.closed() , "requires implementation + testing without closed topologies" );

  // loop through all facets on the boundary of the cavity and retrieve
  // facets which are on the requested geometry
  Entity* e;
  std::vector<index_t> f;
  index_t elem;
  for (index_t k=0;k<nb_cavity();k++)
  {
    if (topology_.master().parameter()) break;
    elem = cavity_[k];
    if (!topology_.ghost(elem)) continue;

    // the first index should be 0 (ghost)
    avro_assert( topology_(elem,0) < topology_.points().nb_ghost() );

    // first check if this is a facet along the requested entity
    // by counting how many points it has on the entity or lower in the hierarchy
    coord_t count = 0;
    for (index_t j=1;j<topology_.nv(elem);j++)
    {
      e = topology_.points().entity( topology_(elem,j) );
      if (!entity->above(e) && e!=entity)
        continue;
      count++;
    }
    // we may not retrieve all geometry facets of the original cavity
    // because there might be situations in which it touches multiple faces
    if (count!=geometry.number()+1)
      continue;


    // get the opposite neighbour
    int k1 = -1;
    for (index_t j=0;j<topology_.neighbours().nfacets();j++)
    {
      if (topology_.ghost(topology_.neighbours()(elem,j)))
        continue;
      k1 = topology_.neighbours()(elem,j);
      break;
    }
    avro_assert( k1 >= 0 );

    // get the index of k1
    int j = topology_.neighbours().oppositeIndex(k1,elem);
    avro_assert( j>=0 );

    // get the oriented facet (TODO, not restricted to tetrahedra...)
    f.resize(3);
    avro_assert_msg(topology_.number()==3,"implement...");
    if (j==0)
    {
      f[0] = topology_(k1,1);
      f[1] = topology_(k1,2);
      f[2] = topology_(k1,3);
    }
    else if (j==1)
    {
      f[0] = topology_(k1,0);
      f[1] = topology_(k1,3);
      f[2] = topology_(k1,2);
    }
    else if (j==2)
    {
      f[0] = topology_(k1,0);
      f[1] = topology_(k1,1);
      f[2] = topology_(k1,3);
    }
    else if (j==3)
    {
      f[0] = topology_(k1,0);
      f[1] = topology_(k1,2);
      f[2] = topology_(k1,1);
    }
    geometry.add(f.data(),f.size());

  } // loop over cavity elements

  if (topology_.master().parameter())
  {
    // add all elements in the cavity
    avro_assert(geometry.nb()==0);
    for (index_t k=0;k<nb_cavity();k++)
    {
      f.resize(3);
      for (index_t j=0;j<3;j++)
        f[j] = topology_( cavity_[k],j );
      geometry.add(f.data(),f.size());
    }
  }

  avro_assert(geometry.nb()>0);
  if (geometry.nb()==0) return;

  // get all the points in the geometry topology
  std::vector<index_t> N = geometry.data();
  uniquify( N );

  avro_assert( v2u.size()==0 );
  avro_assert( u2v.size()==0 );

  // get the parameter coordinates of each vertex
  real_t params[2] = {0,0};;
  for (index_t k=0;k<N.size();k++)
  {
    // create a map from topology points to the geometry topology ones
    v2u.insert( std::pair<index_t,index_t>( N[k] , k ) );
    u2v.push_back( N[k] );

    if (uvcoords)
    {
      // we were requested to create the topology in uv-space
      geometry.points().create( params );
      geometry.points().set_entity(k , topology_.points().entity(N[k]) );
    }
    else
    {
      // we were requested to create the topology in the physical space
      geometry.points().create( this->points_[N[k]] );
      geometry.points().set_param(k , topology_.points().u(N[k]) );
      geometry.points().set_entity(k , topology_.points().entity(N[k]) );
    }
  }

  if (uvcoords)
  {
    coord_t udim = this->points_.udim();
    avro_assert(udim == geometry.points().dim());
    std::vector<real_t> params;
    for (index_t k=0;k<geometry.nb();k++)
    {
      params.resize(udim*geometry.nv(k));
      geometry_params( entity , this->points_ , geometry(k) , geometry.nv(k) , params.data() );
      for (index_t j=0;j<geometry.nv(k);j++)
      {
        if (v2u.find(geometry(k,j))==v2u.end()) avro_assert(false);
        index_t m = v2u.at( geometry(k,j) );
        for (index_t i=0;i<udim;i++)
          geometry.points()[m][i] = params[ udim*j+i ];
      }
    }
  }

  // map the geometry topology indices
  for (index_t k=0;k<geometry.nb();k++)
  for (index_t j=0;j<geometry.nv(k);j++)
  {
    avro_assert( v2u.find(geometry(k,j))!=v2u.end() );
    geometry(k,j) = v2u[ geometry(k,j) ];
  }

  // close the geometry
  geometry.close();

  // compute the neighbours
  geometry.neighbours().compute();

  // increment the map indices to account for the ghost
  avro_assert( geometry.points().nb_ghost()==1 );
  std::map<index_t,index_t>::iterator it;
  for (it=v2u.begin();it!=v2u.end();it++)
    it->second++;
  u2v.insert( u2v.begin() , this->points_.nb()+1 );

  // check consistency between the maps
  for (index_t k=0;k<geometry.points().nb();k++)
  {
    if (k<geometry.points().nb_ghost()) continue;
    index_t v = u2v[k];
    if (v2u.find(v)==v2u.end()) avro_assert_not_reached;
    avro_assert( v2u[v]==k );
    avro_assert( v2u[v]>=geometry.points().nb_ghost() );
  }

}


template<typename type>
bool
Cavity<type>::enlarge( bool verbose )
{
  //const coord_t dim = topology_.points().dim();
  index_t nf = this->number_+1;
  std::vector<const real_t*> xk(nf);

  std::vector<index_t> C; // new cavity elements
  avro_assert_msg( topology_.closed() , "requires implementation + testing without closed topologies" );

  // loop through the current cavity
  for (index_t k=0;k<nb_cavity();k++)
  {
    // dont check facets of ghost elements
    if (topology_.ghost( cavity_[k] ) ) continue;

    // loop through the neighbours (facets) of the cavity element
    for (index_t j=0;j<nf;j++)
    {
      // retrieve the neighbour bordering this facet
      int neighbour = topology_.neighbours()( cavity_[k] , j );

      // for closed topologies there should be no boundary
      // so neighbours should have non-negative indices in the topology
      if (topology_.closed())
        avro_assert( neighbour>=0 );

      // skip neighbours that are already in the cavity
      // because then this is not a boundary facet
      if (contains(index_t(neighbour)) && topology_.closed())
      {
        if (verbose) printf("topology is closed and neighbour %d already in cavity\n",neighbour);
        continue;
      }

      if (!topology_.closed() && neighbour>=0 && contains(index_t(neighbour)))
      {
        if (verbose) printf("topology is not closed and hit interior neighbour %d\n",neighbour);
        continue;
      }

      // check if this facet contains the proposed point
      bool nopenope = false;
      for (index_t i=0;i<nf;i++)
      {
        if (j==i) continue;
        if (topology_( cavity_[k],i )==star_)
        {
          nopenope = true;
          break;
        }
      }
      if (nopenope) continue; // don't count boundary facets of the cavity

#if 0
      // neighbour is not in cavity which means we hit a boundary facet
      // set the coordinates
      for (index_t i=0;i<nf;i++)
        xk[i] = topology_.points()[ topology_( cavity_[k] , i ) ];

      // set the last coordinate to the proposed point
      xk[j] = point_.data();

      // check the orientation
      real_t vol = numerics::simplex_volume(xk,dim);
#else
      real_t vol = get_volume( topology_ , entity_ , cavity_[k] , j , point_.data() );
#endif
      if (sign_*vol<minvol_) // slivers?
      {
        // the facet cannot see the point on the geometry, no hope
        if (topology_.points().entity(star_)!=NULL)
        {
          // the proposed point is on a geometry entity but the facet
          // cannot see the point so the operator should be rejected
          if (!ignore_)
            return false; // TODO is this even necessary anymore?? added on july 2nd for geometry cavities
        }

        // add the neighbour
        if (neighbour>=0)
          C.push_back( index_t(neighbour) );
      }
    }
  }

  // if no elements were added on this pass then we are done
  if (C.size()==0)
  {
    if (verbose) printf("cavity was not enlarged\n");
    return true;
  }

  // specific operator may disallow enlarging further
  if (C.size()>0 && !enlarge_)
  {
    return false;
  }

  // add the new elements to the cavity
  for (index_t k=0;k<C.size();k++)
    add_cavity( C[k] );

  uniquify(cavity_);

  // recursively enlarge
  if (!enlarge()) return false;

  return true;
}

template<typename type>
bool
Cavity<type>::find_geometry( real_t* x , std::vector<index_t>& C0 )
{
  if (topology_.master().parameter()) return true;

  const coord_t dim = topology_.points().dim();
  index_t nf = this->number_+1;
  std::vector<const real_t*> xk(nf);

  std::vector<index_t> C; // new cavity elements

  // loop through the current set of elements
  for (index_t k=0;k<C0.size();k++)
  {

    // dont check facets of ghost elements
    if (topology_.ghost( C0[k] ) ) continue;

    // loop through the neighbours (facets) of the k'th element
    for (index_t j=0;j<nf;j++)
    {
      // get the j'th neighbour
      int neighbour = topology_.neighbours()( C0[k] , j );

      // for closed topologies there should be no boundary
      if (topology_.closed())
        avro_assert( neighbour>=0 );

      // skip neighbours that are already in the cavity
      // because then this is not a boundary facet
      // TODO move this to Set::contains
      bool has = false;
      for (index_t i=0;i<C0.size();i++)
      {
        if (C0[i]==index_t(neighbour))
        {
          has = true;
          break;
        }
      }
      if (has) continue;

      // neighbour is not in cavity which means we hit a boundary facet
      // set the coordinates
      for (index_t i=0;i<nf;i++)
        xk[i] = topology_.points()[ topology_( C0[k] , i ) ];

      // set the j'th coordinate to that of the requested point
      xk[j] = x;

      // check the orientation
      real_t vol = numerics::simplex_volume(xk,dim);
      if (sign_*vol<=minvol_)
      {
        // if the volume is negative for a facet bounding a ghost,
        // then the point is not visible
        if (topology_.ghost( index_t(neighbour) ) )
        {
          return false;
        }

        // add the neighbour
        C.push_back( index_t(neighbour) );
      }
    }
  }

  // if no elements were added on this pass then we are done
  if (C.size()==0) return true;

  // add the new elements to the cavity
  for (index_t k=0;k<C.size();k++)
    C0.push_back( C[k] );
  uniquify(C0);

  // recursively enlarge
  if (!find_geometry(x,C0)) return false;

  return true;
}

template<typename type>
void
Cavity<type>::compute_removed_nodes()
{
  removed_nodes_.clear();
  for (index_t k=0;k<nb_nodes();k++)
  {
    if (this->has(nodes_[k])) continue;
    add_removed_node( nodes_[k] );
  }
}

template<typename type>
bool
Cavity<type>::contains( const index_t c ) const
{
  for (index_t k=0;k<nb_cavity();k++)
  {
    if (cavity_[k]==c)
      return true;
  }
  return false;
}

template<typename type>
void
Cavity<type>::apply()
{
  // clear the mesh topology from a previous candidate (or the original cavity)
  Topology<type>::clear();

  // the inserted element to fill
  std::vector<index_t> t(this->number_+1);

  // connect the boundary da to the candidate star_
  avro_assert( idx_.size()==nb_bnd() );
  for (index_t k=0;k<nb_bnd();k++)
  {
    // skip cavity boundary facets which do not contain the proposed star_
    if (boundary_.has(k,star_)) continue;

    // get the cavity boundary facet and star it to the star_
    index_t count = 0;
    for (index_t j=0;j<index_t(this->number_+1);j++)
    {
      if (j==idx_[k])
        t[j] = star_;
      else
        t[j] = boundary_(k,count++);
    }

    // add the starred element which is oriented correctly
    add( t.data() , t.size() );
  }

  inserted_.resize( nb_insert() );
  removes_.resize( nb_cavity() );
  std::fill( removes_.begin() , removes_.end() , true ); // default

  // determine if any nodes get deleted
  compute_removed_nodes();
}

template<typename type>
void
Cavity<type>::compute_nodes()
{
  avro_assert(nodes_.size()==0);

  // compute the nodes which are initially in the elements of the cavity
  for (index_t k=0;k<nb_cavity();k++)
  {
    for (index_t j=0;j<topology_.nv( cavity_[k] ); j++)
      add_node( topology_(cavity_[k],j) );
  }
  uniquify(nodes_);
}

template<typename type>
void
Cavity<type>::print() const
{
  printf("star = %lu\n",star_);
  printf("info = %s\n",info_.c_str());
  topology_.points().print(star_,true);
  if (entity_!=nullptr) entity_->print();
  for (index_t k=0;k<nb_cavity();k++)
    print_inline( topology_.get(cavity(k)) , "cavity["+stringify(cavity(k))+"]" );
  print_inline( nodes_ , "nodes" );
  print_inline( removed_nodes_ , "removed_nodes" );

  for (index_t k=0;k<nb_insert();k++)
    print_inline( this->get(k) , "insertion["+stringify(k)+"]");

  for (index_t k=0;k<nb_bnd();k++)
    print_inline( boundary_.get(k) , "bnd["+stringify(k)+"]");

  for (index_t k=0;k<nb_nodes();k++)
    topology_.points().print(nodes_[k],true);

  // save the inserted elements (b)
  Topology<Simplex> b(this->points_,this->number_);
  for (index_t k=0;k<nb_cavity();k++)
    b.add( topology_(cavity(k)) , topology_.nv(cavity(k)) );
}

template<typename type>
bool
Cavity<type>::has_unique_elems()
{
  // check that any of the insertions are not contained in the ball of the star
  std::vector<index_t> ball;
  topology_.intersect( {star_} , ball );

  // only track elements which are not in the cavity
  Topology<type> T(topology_.points(),topology_.number());
  for (index_t k=0;k<ball.size();k++)
  {
    if (!contains(ball[k]))
      T.add( topology_(ball[k]) , topology_.nv(ball[k]) );
  }

  for (index_t k=0;k<T.nb();k++)
  {
    if (this->cardinality( T(k) , T.nv(k) )>0)
      return false;
  }
  return true;
}

template<typename type>
void
Cavity<type>::add( const index_t* v , const index_t nv )
{
  Topology<type>::add(v,nv);
}

template<typename type>
void
Cavity<type>::add_node( const index_t node )
{
  nodes_.push_back(node);
}

template<typename type>
void
Cavity<type>::add_removed_node( const index_t node )
{
  removed_nodes_.push_back(node);
}

template<typename type>
void
Cavity<type>::add_cavity( const index_t elem )
{
  cavity_.push_back(elem);
}

template<typename type>
void
Cavity<type>::clear()
{
  Topology<type>::clear();
  cavity_.clear();
  nodes_.clear();
  boundary_.clear();
}

template<typename type>
void
Cavity<type>::copy( const Cavity<type>& cavity )
{
  for (index_t k=0;k<cavity.nb();k++)
    this->add( cavity(k) , cavity.nv(k) );

  for (index_t k=0;k<cavity.nb_cavity();k++)
    this->add_cavity( cavity.cavity(k) );

  //print_inline(cavity_ , "CAVITY" );

  for (index_t k=0;k<cavity.nodes().size();k++)
    add_node( cavity.nodes()[k] );
  inserted_.resize( nb_insert() );
}

template class Cavity<Simplex>;

} // avro
