#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "library/metric.h"

#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include <egads.h>

#include <set>

namespace luna
{

template<typename type>
Entity*
Primitive<type>::geometry( index_t p0 , index_t p1 )
{
  Entity* e0 = this->topology_.points().entity(p0);
  Entity* e1 = this->topology_.points().entity(p1);
  if (e0==NULL || e1==NULL) return NULL;

  Entity* g = e0->intersect(e1);
  if (g==NULL) return NULL;

  if (g->interior()) return g; // skip the ghost check

  // we need to make sure the edge is attached to some ghosts
  std::vector<index_t> shell;
  this->topology_.intersect( {p0,p1} , shell );

  for (index_t k=0;k<shell.size();k++)
  {
    if (this->topology_.ghost(shell[k]))
      return g;
  }
  // there are no ghosts, cannot be a geometry edge
  return NULL;
}

template<typename type>
void
Primitive<type>::extractGeometry( Entity* e , const std::vector<index_t>& f )
{
  luna_assert( e->number()==2 );
  u_.clear();
  G_.clear();
  v2u_.clear();
  u2v_.clear();
  gcavity_.clear();
  S_.clear();
  G_.neighbours().forceCompute(); // forces the neighbours to be recomputed
  G_.set_closed(false); // forces the re-closing of the mesh

  this->computeGeometry( e , G_ , v2u_ , u2v_ );
  if (G_.nb()==0) return;

  G_.inverse().build();
  if (f.size()==0)
  {
    // a size of zero is a code for extracting non-ghost elements
    for (index_t k=0;k<G_.nb();k++)
    {
      if (G_.ghost(k)) continue;
      S_.push_back(k);
    }
  }
  else if (f.size()==1)
  {
    // requested a vertex v, look up the u value
    index_t u = v2u_.at(f[0]);
    G_.inverse().ball(u,S_);
  }
  else if (f.size()==2)
  {
    // requested an edge, lookup non-ghost elements
    index_t u0 = v2u_.at(f[0]);
    index_t u1 = v2u_.at(f[1]);
    G_.inverse().shell(u0,u1,S_);
    luna_assert( S_.size()==2 );
  }
  else
  {
    print_inline(f,"unsupported facet: ");
    luna_assert_not_reached;
  }
}

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

  luna_assert( ep->number()==2 );
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
  luna_assert(this->G_.nb()>0);

  // add the parameter coordinates for the inserted point along with the mapped index
  this->u_.create( params );
  this->u2v_.push_back( this->points_.nb()-1 );

  // the insertion should be visible in parameter space
  bool accept = this->gcavity_.compute( this->u_.nb()-1 , params , this->S_ );
  luna_assert(accept);

  // check the orientation of the facets
  GeometryOrientationChecker checker(this->points_,this->u_,this->u2v_,ep);
  int s = checker.signof( this->gcavity_ );
  if (s<0) return false;

  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype) );

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
EdgeSwap<type>::EdgeSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("edge swapper");
  nb_geometry_rejections_.resize( _topology.number() );
}

template<typename type>
bool
EdgeSwap<type>::visibleParameterSpace( index_t p , index_t e0 , index_t e1 , Entity* face0 )
{
  EGADS::Object* face = (EGADS::Object*) face0;

  // check if the edge swap is valid in the parameter space
  // 2d swaps area already checks, 4d is limited to linear geometries for now
  if (this->topology_.number()!=3) return true;
  if (!this->curved_) return true;
  if (face->interior()) return true;

  // based on the type of face, we may need to flip the sign for the volume calculation
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*face->object(), &oclass, &mtype,&ref, &prev, &next) );
  if (mtype==SREVERSE)
    this->gcavity_.sign() = -1;
  else
    this->gcavity_.sign() = 1;

  // extract the geometry cavity
  this->extractGeometry( face, {e0,e1} );
  luna_assert( this->G_.nb()>0 );

  if (!this->G_.closed())
    luna_assert_msg( this->G_.nb()==2 ,
                    "|G| = %lu, |swap ghosts| = %lu" , this->G_.nb() , this->nb_ghost() );

  if (this->G_.closed())
    luna_assert( this->G_.nb_real()==2 );

  // ensure the e0 is visible to the cavity boundary
  bool accept = this->gcavity_.compute( this->v2u_[e0] , this->u_[this->v2u_[e0]] , this->S_ );
  luna_assert( accept );

  GeometryOrientationChecker checker( this->points_ , this->u_ , this->u2v_ , face );
  int s = checker.signof( this->gcavity_ );
  luna_assert_msg( s > 0 , "negative orientation for edge (%lu,%lu) with vertex %lu" , e0,e1,e0 );
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  // ensure e1 is visible to the cavity boundary
  accept = this->gcavity_.compute( this->v2u_[e1] , this->u_[this->v2u_[e1]] , this->S_ );
  luna_assert( accept );
  s = checker.signof( this->gcavity_ );
  luna_assert_msg( s > 0 , "negative orientation for edge (%lu,%lu) with vertex %lu" , e0,e1,e1 );
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  // check for visibility of p
  nb_parameter_tests_++;
  accept = this->gcavity_.compute( this->v2u_[p] , this->u_[this->v2u_[p]] , this->S_ );
  luna_assert( this->gcavity_.nb_real()==2 );
  if (!accept)
  {
    //printf("swap not visible in parameter space!\n");
    nb_parameter_rejections_++;
    return false;
  }

  s = checker.signof( this->gcavity_ );
  if (s<0)
  {
    //printf("swap creates negative normals!\n");
    return false;
  }

  if (checker.createsBadGeometry(this->gcavity_))
  {
    //printf("swap creates bad geometry!\n");
    nb_invalid_geometry_++;
    return false;
  }

  // the geometry cavity should have positive volumes
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  return true;
}

template<typename type>
bool
EdgeSwap<type>::valid( const index_t p , const index_t e0 , const index_t e1 )
{
  // topology checks
  if (this->topology_.points().fixed(e0) ||
      this->topology_.points().fixed(e1))
    return false;
  if (p<this->topology_.points().nb_ghost()) return false;
  if (e0<this->topology_.points().nb_ghost()) return false;
  if (e1<this->topology_.points().nb_ghost()) return false;

  std::vector<index_t> edge = {e0,e1};
  if (this->C_.empty())
    this->topology_.intersect(edge,this->C_);

  std::vector<index_t>& elems = this->C_;

  for (index_t k=0;k<elems.size();k++)
  for (index_t j=0;j<this->topology_.nv(elems[k]);j++)
  {
    if (this->topology_.points().fixed(this->topology_(elems[k],j)))
      return false;
  }

  // check if the two points lie on an Edge entity
  Entity* ge = this->geometry(e0,e1);

  // check if this edge is in the volume
  if (ge==NULL)
  {
    // edge is in the volume, always valid
    return true;
  }

  // check if this edge is on a geometry Edge (can't swap these)
  if (ge!=NULL && ge->number()<2)
  {
    nb_geometry_rejections_[ge->number()]++;
    return false;
  }

  // check the entity on the re-inserted vertex
  if (ge!=NULL)
  {
    Entity* ep = this->topology_.points().entity(p);
    if (ep==NULL) return false;
    if (ep!=ge)
    {
      if (ge->above(ep))
      {} // the swap is okay
      else
      {
        nb_geometry_rejections_[ge->number()]++;
        return false;
      }
    }
  }

  if (ge->interior())
  {
    nb_wake_++;
  }

  return true;
}

template<typename type>
bool
EdgeSwap<type>::apply( const index_t p , const index_t e0 , const index_t e1 )
{
  // check if the swap is valid in terms of geometry
  if (!valid(p,e0,e1)) return false; // computes cavity C_

  this->info_ = "trying to swap edge (" + stringify<index_t>(e0) + "/"
                + stringify<index_t>(e1) +")"+
                " with reinsertion " + stringify<index_t>(p);

  // compute the original number of ghost elements
  index_t nb_ghost0 = 0;
  for (index_t j=0;j<this->C_.size();j++)
    if (this->topology_.ghost(this->C_[j]))
      nb_ghost0++;

  // apply the operator, checking visibility
  this->enlarge_ = false;
  bool accept = this->compute( p , this->topology_.points()[p] , this->C_ );
  if (!accept) return false;

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
    return false;

  // check visibility in the parametric space
  Entity* ge = this->geometry(e0,e1);
  if (ge!=NULL && ge->number()==2)
  {
    if (!visibleParameterSpace(p,e0,e1,ge))
    {
      //printf("swap not visible in parameter space!\n");
      return false;
    }
  }

  // count the resulting number of ghost elements
  index_t nb_ghost = 0;
  for (index_t j=0;j<this->nb();j++)
    if (this->ghost(j))
      nb_ghost++;

  // check if only ghosts are created
  index_t only_ghost = true;
  for (index_t k=0;k<this->nb();k++)
  {
    if (!this->ghost(k))
    {
      only_ghost = false;
      break;
    }
  }
  if (only_ghost) return false;

  // do not allow swaps which change the number of ghosts for 3-simplices
  if (this->topology_.number()==3 && nb_ghost0!=nb_ghost)
    return false;

  if (!accept) return false;

  return accept;
}

template<typename type>
RidgeSwap<type>::RidgeSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("ridge swapper");
  nb_geometry_rejections_.resize( _topology.number() );
}

template<typename type>
bool
RidgeSwap<type>::valid( index_t p , const index_t t0 , const index_t t1 , const index_t t2 )
{
  Entity* ep = this->topology_.points().entity(p);
  Entity* e0 = this->topology_.points().entity(t0);
  Entity* e1 = this->topology_.points().entity(t1);
  Entity* e2 = this->topology_.points().entity(t2);

  // for now disregard any ridge swap touch a boundary
  if (ep!=NULL || e0!=NULL || e1!=NULL || e2!=NULL) return false;
  return true;
}

template<typename type>
bool
RidgeSwap<type>::apply( index_t p , const index_t t0 , const index_t t1 , const index_t t2 )
{
  luna_assert( this->topology_.number()==4 ); // should only be used in 4d though technically a face swap in 3d
  if (!valid(p,t0,t1,t2)) return false;

  // get all elements touching this triangle
  this->topology_.intersect( {t0,t1,t2} , this->C_ );

  this->enlarge_ = false;
  bool accept = this->compute( p , this->topology_.points()[p] , this->C_ );

  // check if any points are deleted
  if (this->nb_removedNodes()>0) return false;

  return accept;
}

template<typename type>
FacetSwap<type>::FacetSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("facet swapper");
}

template<typename type>
bool
FacetSwap<type>::valid( const index_t p , const index_t k0 , const index_t k1 )
{
  if (p<=this->topology_.points().nb_ghost()) return false;
  if (this->topology_.ghost(k0) || this->topology_.ghost(k1))
    return false;
  return true;
}

template<typename type>
bool
FacetSwap<type>::apply( const index_t p , const index_t k1 , const index_t k2 )
{
  // check if the facet swap is valid
  if (!valid(p,k1,k2)) return false;
  this->enlarge_ = false;

  // check if the swap is visible
  std::vector<index_t> elems = {k1,k2};
  bool accept = this->compute( p , this->topology_.points()[p] , elems );

  // check if all produced elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
    return false;
  return accept;
}

template<typename type>
Smooth<type>::Smooth( Topology<type>& _topology ) :
  Primitive<type >(_topology),
  delta_(0.0),
  delta_min_(1e20),
  delta_max_(-1),
  M0_(_topology.number())
{
  this->setName("smoother");
  resetRejections();
  nb_accepted_ = 0;
}

template<typename type>
bool
Smooth<type>::visibleParameterSpace( index_t p , real_t* x , real_t* params , Entity* ep0 )
{
  EGADS::Object* ep = (EGADS::Object*) ep0;

  luna_assert( ep->number()==2 );
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
  this->extractGeometry( ep , {p} );
  luna_assert( this->G_.nb()>0 );
  //luna_assert_msg( this->G_.nb_real()==this->nb_ghost() , "|g| = %lu, |c| = %lu\n" , this->G_.nb_real() , this->nb_ghost() );

  for (coord_t d=0;d<2;d++)
    this->u_[this->v2u_[p]][d] = params[d];

  // the following assertion won't hold when p is on an Edge
  if (this->topology_.points().entity(p)->number()==2 && !this->G_.closed())
  {
    luna_assert_msg( this->G_.nb()==this->nb_ghost() ,
                      "|G| = %lu, |smooth ghosts| = %lu" , this->G_.nb() , this->nb_ghost() );
  }

  GeometryOrientationChecker checker( this->points_ , this->u_ , this->u2v_ , ep );

  // check for visibility of the new coordinates
  nb_parameter_tests_++;
  bool accept = this->gcavity_.compute( this->v2u_[p] , params , this->S_ );
  if (!accept)
  {
    nb_parameter_rejections_++;
    return false;
  }

  // check the orientation of facets produced w.r.t. the geometry
  int s = checker.signof( this->gcavity_ );
  if (s<0) return false;

  // check the produced volumes are positive
  luna_assert( checker.hasPositiveVolumes(this->gcavity_,mtype));

  return true;
}

template<typename type>
bool
Smooth<type>::apply( const index_t p , MetricField<type>& metric , real_t Q0 )
{
  // compute the cavity around p
  this->C_.clear();
  this->topology_.intersect( {p} , this->C_ );

  const coord_t dim = this->topology_.points().dim();
  std::vector<real_t> F( dim,0. );
  std::vector<real_t> u( dim,0. );
  std::vector<real_t> x( dim,0. );
  std::vector<real_t> x0 (dim,0. );

  const coord_t udim = this->topology_.points().udim();
  std::vector<real_t> params0( udim,0. );
  std::vector<real_t> params( udim,0. );

  // get all the nodes in the cavity (excluding this one)
  std::set<index_t> Ns;
  index_t q;
  Entity* ep = this->topology_.points().entity(p);
  if (ep!=NULL && ep->number()==0) return false; // don't move points on Nodes!
  Entity* eq;
  for (index_t k=0;k<this->C_.size();k++)
  {
    for (index_t j=0;j<this->topology_.nv(this->C_[k]);j++)
    {
      q = this->topology_(this->C_[k],j);
      if (q==p) continue;
      if (q<this->topology_.points().nb_ghost()) continue;
      eq = this->topology_.points().entity(q);
      if (ep!=NULL)
      {
        // for geometry points, the smoothing can only be affected by
        // geometry points which are equal to or higher in the hierarchy
        // than the geometry of the smoothed vertex
        // todo: use "geometryEdge"
        if (eq==NULL) continue;
        if (ep!=eq)
        {
          // the entities are different
          if (!ep->above(eq))
          {
            // ep must at least be higher in the topological hierarchy
            continue; // skip vertex q
          }
        }
      }
      Ns.insert(q);
    }
  }

  if (Ns.size()==0)
  {
    // no points were found to do the smoothing
    // this really should not happen but just in case
    nb_zero_valency_++;
    return false;
  }

  // copy the set of points into the vector of points N
  std::vector<index_t> N;
  for (std::set<index_t>::iterator it=Ns.begin();it!=Ns.end();++it)
    N.push_back(*it);
  Ntot_ += N.size();

  // compute the tensile force
  std::vector<real_t> lens0;
  real_t len;
  real_t f;
  for (index_t k=0;k<N.size();k++)
  {
    // compute the metric length
    len = metric.length(this->topology_.points(),p,N[k]);
    lens0.push_back(len);

    // physical vector from the vertex to the neighbour
    numerics::vector( this->topology_.points()[p] ,
                        this->topology_.points()[N[k]] , dim , u.data() );

    // compute the force on the vertex
    f = std::pow(len,4);
    f = (1. -len)*std::exp(-len);

    for (coord_t d=0;d<dim;d++)
      F[d] += -f*u[d];
  }

  for (coord_t d=0;d<dim;d++)
    objective_ += std::pow( F[d] , 2. );

  // compute the new position of the vertex
  real_t omega = 0.2; // relaxation factor
  for (index_t d=0;d<dim;d++)
  {
    x0[d] = this->topology_.points()[p][d];
    x[d]  = this->topology_.points()[p][d] +omega*F[d];
  }
  real_t delta_p = numerics::distance( x0.data() , x.data() , dim );
  delta_ += delta_p;

  if (delta_p<delta_min_) delta_min_ = delta_p;
  if (delta_p>delta_max_) delta_max_ = delta_p;

  this->enlarge_ = true;

  // check if the cavity needs to be enlarged
  if (ep!=NULL)
  {
    for (index_t d=0;d<udim;d++)
      params0[d] = params[d] = this->topology_.points().u(p,d);

    // use the previous parameter value as the initial guess of the projection
    /*if (!this->curved_)
      ep->project(x);
    else
      ep->projectGuess(x,params);
    */
    if (!this->curved_)
      ep->inverse(x,params);
    else
    {
      ep->inverse_guess(x,params);

      // check the re-evaluated parameter coordinates are close to the new x coordinates
      #if 0
      real_t result[18];
      EGADS_ENSURE_SUCCESS( EG_evaluate(*ep->object(),params.data(),result));
      real_t dg = numerics::distance(result,x.data(),dim);
      real_t gtol = 0;
      EGADS_ENSURE_SUCCESS( EG_tolerance(*ep->object(),&gtol) );
      luna_assert_msg( dg < gtol , "x = (%g,%g,%g), x(u) = (%g,%g,%g), distance = %g\n",
                        x[0],x[1],x[2],result[0],result[1],result[2],dg);
      #endif
    }

    // make sure the cavity is not enlarged
    std::vector<index_t> C0 = this->C_;
    this->findGeometry( x.data() , this->C_ );
    if (this->C_.size()!=C0.size())
    {
      nb_enlarged_rejections_++;
      return false;
    }

    // there are more efficient ways of doing this but this is okay for now...
    std::sort( this->C_.begin() , this->C_.end() );
    std::sort( C0.begin() , C0.end() );
    for (index_t j=0;j<C0.size();j++)
    {
      if (this->C_[j]!=C0[j])
      {
        nb_enlarged_rejections_++;
        return false;
      }
    }

    nb_geometry_++;
  }

  this->enlarge_ = false;
  bool accept = this->compute( p , x.data() , this->C_ );
  if (this->nb_removedNodes()>0)
  {
    nb_removed_rejections_++;
    return false;
  }
  if (!accept)
  {
    // the point was not visible within the cavity
    nb_visibility_rejections_++;
    return false;
  }

  // apply the physical coordinates
  for (coord_t d=0;d<dim;d++)
    this->topology_.points()[p][d] = x[d];

  // apply the parameter space coordinates if necesary
  if (ep!=NULL)
  {
    for (coord_t d=0;d<udim;d++)
      this->topology_.points().u(p,d) = params[d];
  }

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positiveImpliedMetrics())
  {
    nb_implied_metric_rejections_++;

    // revert the physical coordinates
    for (index_t d=0;d<dim;d++)
      this->topology_.points()[p][d] = x0[d];

    // revert the parameter space coordinates
    for (coord_t d=0;d<udim;d++)
      this->topology_.points().u(p,d) = params0[d];

    // signal the smoothing swas not applied
    return false;
  }

  // check if the point is visible on the cavity boundary
  if (ep!=NULL && ep->number()==2 && this->curved_)
  {
    if (!visibleParameterSpace( p , x.data(), params.data() , ep ))
    {
      // revert the physical coordinates
      for (index_t d=0;d<dim;d++)
        this->topology_.points()[p][d] = x0[d];

      // revert the parameter space coordinates
      for (coord_t d=0;d<udim;d++)
        this->topology_.points().u(p,d) = params0[d];

      // signal the smoothing was not applied
      return false;
    }
  }

  if (ep!=NULL && this->topology_.number()==3 && ep->number()==1 && this->curved_)
  {
    // we need to check for visibility on all Faces that parent this Edge
    for (index_t k=0;k<ep->nb_parents();k++)
    {
      Entity* parent = ep->parents(k);
      if (parent->number()!=2 || !parent->tessellatable())
        continue;
      if (!visibleParameterSpace(p, x.data() , params.data() ,parent))
      {
        // revert the physical coordinates
        for (index_t d=0;d<dim;d++)
          this->topology_.points()[p][d] = x0[d];

        // revert the parameter space coordinates
        for (coord_t d=0;d<udim;d++)
          this->topology_.points().u(p,d) = params0[d];

        // signal the smoothing was not applied
        return false;
      }
    }
  }

  this->enlarge_ = false;

  // save the previous metric
  for (index_t i=0;i<dim;i++)
  for (index_t j=i;j<dim;j++)
    M0_(i,j) = metric(this->topology_.points(),p)(i,j);
  index_t elem0 = metric.value(p).elem();

  // recompute the metric at the new point
  bool success = metric.recompute( p , x.data() );

  // count this metric interpolation as "outside"
  if (!success)
    nb_interpolated_outside_++;

  if (!success)
  {
    for (index_t d=0;d<dim;d++)
      this->topology_.points()[p][d] = x0[d];

    if (ep!=NULL)
      for (coord_t d=0;d<udim;d++)
        this->topology_.points().u(p,d) = params0[d];

    metric.attachment().assign( p , M0_ , elem0 );
    return false;
  }

  // evaluate the new quality
  if (Q0>0.0)
  {
    if (worst_quality(*this,metric)<Q0)
    {
      // revert the coordinates and re-assign the metric
      for (index_t d=0;d<dim;d++)
        this->topology_.points()[p][d] = x0[d];

      if (ep!=NULL)
        for (coord_t d=0;d<udim;d++)
          this->topology_.points().u(p,d) = params0[d];

      metric.attachment().assign( p , M0_ , elem0 );
    }
  }
  nb_accepted_++;

  return true;
}
template class Insert<Simplex>;
template class EdgeSwap<Simplex>;
template class FacetSwap<Simplex>;
template class RidgeSwap<Simplex>;
template class Smooth<Simplex>;
template class Primitive<Simplex>;

} // luna
