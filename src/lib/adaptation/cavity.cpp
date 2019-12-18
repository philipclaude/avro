#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "adaptation/cavity.h"

#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/geometry.h"

#include <egads.h>

#include <algorithm>
#include <cstdio>
#include <array>
#include <cmath>
#include <limits>

namespace luna
{

template<typename type>
Cavity<type>::Cavity( Topology<type>& _topology ) :
  Topology<type>(_topology.points(),_topology.number()),
  topology_(_topology),
  nodeRemovalAllowed_(true),
  enlarge_(true),
  checkVisibility_(true),
  star_(0),
  boundary_( topology_.points() , topology_.number()-1 ),
  minvol_(1e-12),
  sign_(1.0),
  nb_error_(0),
  rethrow_(false),
  mk_(topology_.master()),
  ignore_(false)
{
  // allocate enough space for the re-inserted point coordinates
  point_.resize( topology_.points().dim() );

  // we do not want to store the elements because they will be oriented
  // to maintain a positive volume
  boundary_.set_sorted(false);
  this->set_sorted(false);
}

template<typename type>
bool
Cavity<type>::positiveImpliedMetrics()
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
  removedNodes_.clear();

  // assign the coordinates and index of the star
  for (coord_t d=0;d<point_.size();d++)
    point_[d] = x[d];
  star_ = p;

  // compute the cavity about p with elements initialized to C0
  for (index_t k=0;k<C0.size();k++)
    addCavity(C0[k]);

  // enlarge the cavity until the point is visible from every face of the boundary
  if (checkVisibility_)
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
    luna_assert(C0.size()>0);
  }

  // accumulate the nodes of the cavity elements
  // this is needed to detect if nodes are removed
  computeNodes();

  // compute the boundary of the cavity elements
  try
  {
    if (!computeBoundary())
    {
      luna_assert_not_reached;
      return false;
    }
  }
  catch(...)
  {
    nb_error_++;
    if (rethrow_) // option for debugging
      luna_assert_not_reached; // any kind of failed assertion
    return false;
  }

  // sort the elements in the cavity,
  // necessary for the application of the operator by the mesh topology
  std::sort( cavity_.begin() , cavity_.end() );

  // apply the cavity to the candidate (also computes removed nodes)
  apply();

  if (!nodeRemovalAllowed_ && nb_removedNodes()>0)
  {
    printf("node removal not allowed\n");
    return false;
  }

  return true;
}

template<typename type>
bool
Cavity<type>::computeBoundary()
{
  luna_assert( boundary_.nb()==0 );
  luna_assert_msg( topology_.closed() , "requires implementation + testing without closed topologies" );

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
        luna_assert( neighbour>=0 );

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


// retrieves the geometric parameter coordinates given a facet on the entity
void
geometryParams( Entity* e0 , const Points& points , const index_t* v , const index_t nv , real_t* params )
{
  luna_assert( nv>0 );

  coord_t udim = points.udim();

  EGADS::Object* e = (EGADS::Object*) e0;

  // get all the entities these points are on
  index_t count = 0;
  std::vector<EGADS::Object*> entities( nv );
  std::vector<bool> ufound( nv, false );
  for (index_t k=0;k<nv;k++)
  {
    entities[k] = (EGADS::Object*) points.entity( v[k] );

    // if the vertex entity is that of the requested one, set the uv coords
    if (entities[k]==e)
    {
      for (int i=0;i<udim;i++)
         params[k*udim+i] = points.u( v[k] )[i];
      ufound[k] = true;
      count++;
    }

    luna_assert(entities[k] != nil); // all points must have entities
  }

  if (count==nv) return; // all parameters have been retrieved

  //---------------------------------------//
  if (e->object_class()==FACE)
  {
    luna_assert(udim == 2);

    std::vector<std::array<double,2>> uvp(nv), uvm(nv);

    EGADS::Object* face = e;
    for (index_t k=0;k<nv;k++)
    {
      // skip any parametric coordinates already set
      if (ufound[k]) continue;

      // try to find an edge and t-value to get uv
      double t = 0;
      EGADS::Object* edge = nil;

      if (entities[k]->object_class()==NODE)
      {
        // find an edge that is a parent of the node and a child of the face
        EGADS::Object* node = entities[k];
        for (index_t i=0;i<node->nb_parents();i++)
        {
          edge = (EGADS::Object*) node->parents(i);

          if (edge->object_class()!=EDGE) continue;

          // skip edges that define periodicity in parameter space if possible
          if (edge->sense_required() && i<node->nb_parents()-1 ) continue;

          // if the face is a parent of the this edge
          if (edge->has_parent(face))
          {
            int periodic;
            double trange[2];
            EGADS_ENSURE_SUCCESS( EG_getRange(*edge->object(), trange, &periodic) );

            #if 0
            if (periodic==1)
            {
              t = trange[0]; // TODO: Is this right?!?!?

              // @marshall: changed to the same logic as non-periodic
              if (edge->egchild(0) == node->object()) t = trange[0];
              if (edge->egchild(1) == node->object()) t = trange[1];
            }
            else
            {
              if (edge->egchild(0) == node->object()) t = trange[0];
              if (edge->egchild(1) == node->object()) t = trange[1];
            }
            #else
            if (edge->egchild(0) == node->object()) t = trange[0];
            if (edge->egchild(1) == node->object()) t = trange[1];
            #endif

            break;
          }
          else
            edge = nil;
        }
      }
      else if (entities[k]->object_class()==EDGE)
      {
        edge = entities[k];
        t = points.u( v[k] )[0];
      }

      if (edge->sense_required())
      {
        // need to sort out which uv should be used later
        EGADS_ENSURE_SUCCESS( EG_getEdgeUV( *face->object() , *edge->object() , -1 , t , uvm[k].data() ) );
        EGADS_ENSURE_SUCCESS( EG_getEdgeUV( *face->object() , *edge->object() ,  1 , t , uvp[k].data() ) );
      }
      else
      {
        // get the uv value from the edge t value. No need to worry about periodicity in uv
        double uv[2] = {0,0};
        int status = EG_getEdgeUV( *face->object() , *edge->object() , 0 , t , uv );

        // special treatment of interior edges that might not be connected to the face
        if (status == EGADS_NOTFOUND && edge->interior())
        {
          bool foundGuess = false;
          for (index_t i=0;i<nv;i++)
          {
            if (!ufound[i]) continue;
            uv[0] = params[i*udim  ];
            uv[1] = params[i*udim+1];
            foundGuess = true;
            break;
          }

          double x_inv[3];
          double* x = const_cast<double*>(points[ v[k] ]);
          if (foundGuess)
          {
            EGADS_ENSURE_SUCCESS( EG_invEvaluateGuess( *face->object() , x , uv , x_inv ) );
          }
          else
          {
            EGADS_ENSURE_SUCCESS( EG_invEvaluate( *face->object() , x , uv , x_inv ) );
          }

          real_t d = numerics::distance2(x,x_inv,3);
          real_t tol_edge=0, tol_face=0;
          EGADS_ENSURE_SUCCESS( EG_tolerance( *edge->object(), &tol_edge ) );
          EGADS_ENSURE_SUCCESS( EG_getTolerance( *face->object(), &tol_face ) );
          luna_assert_msg( d <= std::min(tol_edge,tol_face) ,
                           "d = %1.16e, tol_edge = %1.16e, tol_face = %1.16e" , d , tol_edge , tol_face );
        }
        else
          EGADS_ENSURE_SUCCESS( status );

        params[k*udim  ] = uv[0];
        params[k*udim+1] = uv[1];

        ufound[k] = true;
        count++;
      }
    }
    luna_assert(count > 0); // at least one uv value must be set...

    for (index_t k=0;k<nv;k++)
    {
      if (ufound[k]) continue;

      // get a uv-value that has already been set
      double uv[2] = {0,0};
      for (index_t j=0;j<nv;j++)
      {
        if (j==k || !ufound[k]) continue;
        uv[0] = params[j*udim  ];
        uv[1] = params[j*udim+1];
        break;
      }

      double mindist = std::numeric_limits<double>::max();
      index_t jmin = nv;
      int sense = 0;
      for (index_t j=0;j<nv;j++)
      {
        if (j == k) continue;

        double distm = sqrt(pow(uvm[j][0] - uv[0], 2) + pow(uvm[j][1] - uv[1], 2));
        double distp = sqrt(pow(uvp[j][0] - uv[0], 2) + pow(uvp[j][1] - uv[1], 2));
        if (distm < mindist)
        {
          mindist = distm;
          jmin = j;
          sense = -1;
        }
        if (distp < mindist)
        {
          mindist = distp;
          jmin = j;
          sense = 1;
        }
      }

      if (sense==1)
      {
        params[k*udim  ] = uvp[jmin][0];
        params[k*udim+1] = uvp[jmin][1];
      }
      else
      {
        params[k*udim  ] = uvm[jmin][0];
        params[k*udim+1] = uvm[jmin][1];
      }
      ufound[k] = true;
      count++;

      if (count == nv) break; // all parameters found
    }
  }

  //---------------------------------------//
  else if (e->object_class()==EDGE)
  {
    luna_assert(udim >= 1);

    EGADS::Object* edge = e;

    for (index_t k=0;k<nv;k++)
    {
      // skip any parametric coordinates already set
      if (ufound[k]) continue;

      // try to find an edge and t-value to get uv
      double t = 0;

      if (entities[k]->object_class()==NODE)
      {
        // find an edge that is a parent of the node and a child of the face
        EGADS::Object* node = entities[k];

        int periodic;
        double trange[2];
        EGADS_ENSURE_SUCCESS( EG_getRange(*edge->object(), trange, &periodic) );
        #if 0
        if (periodic==1)
        {
          t = trange[0]; // TODO: Is this right?!?!?

          // @marshall: changed to the same logic as non-periodic
          if (edge->egchild(0) == node->object()) t = trange[0];
          if (edge->egchild(1) == node->object()) t = trange[1];
        }
        else
        {
          if (edge->egchild(0) == node->object()) t = trange[0];
          if (edge->egchild(1) == node->object()) t = trange[1];
        }
        #else
        if (edge->egchild(0) == node->object()) t = trange[0];
        if (edge->egchild(1) == node->object()) t = trange[1];
        #endif

        // set the parametric value
        params[k*udim] = t;

        ufound[k] = true;
        count++;
      }
      else
        luna_assert(false); // this should not happen...

      if (count == nv) break; // all parameters found
    }
  }

  luna_assert(count == nv);

  coord_t dim  = points.dim();

  for (index_t k=0;k<nv;k++)
  {
    // evaluate the coordinates for the parameters we found
    double x_eval[18];
    EGADS_ENSURE_SUCCESS( EG_evaluate( *e->object() , &params[k*udim] , x_eval ) );

    // check the evaluated coordinates are close to the true ones
    const real_t* x = points[ v[k] ];
    real_t d = numerics::distance2(x,x_eval,dim);
    real_t tol = 0;
    EGADS_ENSURE_SUCCESS( EG_getTolerance( *e->object(), &tol ) );

    if (d>tol)
    {
      printf("Failure on facet vertex: %lu\n",v[k]);
      if (nv == 2)
        printf("facet = (%lu,%lu)\n",v[0],v[1]);
      else if (nv == 3)
        printf("facet = (%lu,%lu,%lu)\n",v[0],v[1],v[2]);
      printf("x_true = (%g,%g,%g), x_eval = (%g,%g,%g) with param coordinates (%g,%g)\n",
                x[0],x[1],x[2],x_eval[0],x_eval[1],x_eval[2],params[k*udim],params[k*udim+1]);
      printf("Getting parameters for entity:\n");
      e->print();
      printf("Entity storing vertex with param = (%g,%g):\n", points.u( v[k] )[0], points.u( v[k] )[1]);
      entities[k]->print();
    }

    luna_assert_msg( d <= tol , "d = %1.16e, tol = %1.16e" , d , tol );
  }
}

template<typename type>
void
Cavity<type>::computeGeometry( Entity* entity0 , Topology<type>& geometry , std::map<index_t,index_t>& v2u , std::vector<index_t>& u2v )
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
  luna_assert_msg( topology_.closed() , "requires implementation + testing without closed topologies" );

  // loop through all facets on the boundary of the cavity and retrieve
  // facets which are on the requested geometry
  Entity* e;
  std::vector<index_t> f;
  index_t elem;
  for (index_t k=0;k<nb_cavity();k++)
  {
    elem = cavity_[k];
    if (!topology_.ghost(elem)) continue;

    // the first index should be 0 (ghost)
    luna_assert( topology_(elem,0) < topology_.points().nb_ghost() );

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
    luna_assert( k1 >= 0 );

    // get the index of k1
    int j = topology_.neighbours().oppositeIndex(k1,elem);
    luna_assert( j>=0 );

    // get the oriented facet (TODO, not restricted to tetrahedra...)
    f.resize(3);
    luna_assert_msg(topology_.number()==3,"implement...");
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

  luna_assert(geometry.nb()>0);
  if (geometry.nb()==0) return;

  // get all the points in the geometry topology
  std::vector<index_t> N = geometry.data();
  uniquify( N );

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
    }
  }

  if (uvcoords)
  {
    coord_t udim = this->points_.udim();
    luna_assert(udim == geometry.points().dim());
    std::vector<real_t> params;
    for (index_t k=0;k<geometry.nb();k++)
    {
      params.resize(udim*geometry.nv(k));
      geometryParams( entity , this->points_ , geometry(k) , geometry.nv(k) , params.data() );
      for (index_t j=0;j<geometry.nv(k);j++)
      {
        index_t m = v2u[ geometry(k,j) ];
        for (index_t i=0;i<udim;i++)
          geometry.points()[m][i] = params[ udim*j+i ];
      }
    }
  }

  // map the geometry topology indices
  for (index_t k=0;k<geometry.nb();k++)
  for (index_t j=0;j<geometry.nv(k);j++)
    geometry(k,j) = v2u[ geometry(k,j) ];

  // close the geometry
  geometry.close();

  // compute the neighbours
  geometry.neighbours().compute();

  // increment the map indices to account for the ghost
  luna_assert( geometry.points().nb_ghost()==1 );
  std::map<index_t,index_t>::iterator it;
  for (it=v2u.begin();it!=v2u.end();it++)
    it->second++;
  u2v.insert( u2v.begin() , this->points_.nb()+1 );

  // check consistency between the maps
  for (index_t k=0;k<geometry.points().nb();k++)
  {
    if (k<geometry.points().nb_ghost()) continue;
    index_t v = u2v.at(k);
    luna_assert( v2u[v]==k );
    luna_assert( v2u[v]>=geometry.points().nb_ghost() );
  }

}

template<typename type>
bool
Cavity<type>::enlarge( bool verbose )
{
  const coord_t dim = topology_.points().dim();
  index_t nf = this->number_+1;
  std::vector<const real_t*> xk(nf);

  std::vector<index_t> C; // new cavity elements
  luna_assert_msg( topology_.closed() , "requires implementation + testing without closed topologies" );

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
        luna_assert( neighbour>=0 );

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

      // neighbour is not in cavity which means we hit a boundary facet
      // set the coordinates
      for (index_t i=0;i<nf;i++)
        xk[i] = topology_.points()[ topology_( cavity_[k] , i ) ];

      // set the last coordinate to the proposed point
      xk[j] = point_.data();

      // check the orientation
      real_t vol = numerics::simplex_volume(xk,dim);
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
    addCavity( C[k] );

  uniquify(cavity_);

  // recursively enlarge
  if (!enlarge()) return false;

  return true;
}

template<typename type>
bool
Cavity<type>::findGeometry( real_t* x , std::vector<index_t>& C0 )
{

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
        luna_assert( neighbour>=0 );

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
  if (!findGeometry(x,C0)) return false;

  return true;
}

template<typename type>
void
Cavity<type>::computeRemovedNodes()
{
  removedNodes_.clear();
  for (index_t k=0;k<nb_nodes();k++)
  {
    if (this->has(nodes_[k])) continue;
    addRemovedNode( nodes_[k] );
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
  luna_assert( idx_.size()==nb_bnd() );
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
  computeRemovedNodes();
}

template<typename type>
void
Cavity<type>::computeNodes()
{
  luna_assert(nodes_.size()==0);

  // compute the nodes which are initially in the elements of the cavity
  for (index_t k=0;k<nb_cavity();k++)
  {
    for (index_t j=0;j<topology_.nv( cavity_[k] ); j++)
      addNode( topology_(cavity_[k],j) );
  }
  uniquify(nodes_);
}

template<typename type>
void
Cavity<type>::print() const
{
  printf("star = %lu\n",star_);
  printf("info = %s\n",info_.c_str());
  topology_.points().print(star_);
  for (index_t k=0;k<nb_cavity();k++)
    print_inline( topology_.get(cavity(k)) , "cavity["+stringify(cavity(k))+"]" );
  print_inline( nodes_ , "nodes" );
  print_inline( removedNodes_ , "removedNodes" );

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
Cavity<type>::philipcondition()
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
Cavity<type>::addNode( const index_t node )
{
  nodes_.push_back(node);
}

template<typename type>
void
Cavity<type>::addRemovedNode( const index_t node )
{
  removedNodes_.push_back(node);
}

template<typename type>
void
Cavity<type>::addCavity( const index_t elem )
{
  cavity_.push_back(elem);
}

template<typename type>
void
Cavity<type>::clear()
{
  Topology<type>::clear();
}

template class Cavity<Simplex>;

} // luna
