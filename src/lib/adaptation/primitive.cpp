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

template class Insert<Simplex>;
template class FacetSwap<Simplex>;
template class RidgeSwap<Simplex>;
template class Primitive<Simplex>;

} // luna
