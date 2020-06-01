#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "library/metric.h"

#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include <egads.h>

#include <set>

namespace avro
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

  if (this->topology_.master().parameter()) return g;

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
Primitive<type>::extract_geometry( Entity* e , const std::vector<index_t>& f )
{
  avro_assert( e->number()==2 );
  u_.clear();
  G_.clear();
  v2u_.clear();
  u2v_.clear();
  gcavity_.clear();
  S_.clear();
  G_.neighbours().forceCompute(); // forces the neighbours to be recomputed
  G_.set_closed(false); // forces the re-closing of the mesh

  this->compute_geometry( e , G_ , v2u_ , u2v_ );
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
    avro_assert( S_.size()==2 );
  }
  else
  {
    print_inline(f,"unsupported facet: ");
    avro_assert_not_reached;
  }
}

template<typename type>
bool
SurfaceCavity<type>::visible( index_t p )
{
  geometry_params( geometry_ , topology_.points() , &p , 1 , topology_.points()[p] );
  return this->compute( p , topology_.points()[p] , cavity_ );
}

template<typename type>
bool
SurfaceCavity<type>::check_normals()
{
  params_.clear();

  // determine the parameter coordinates along the geometry
  params0_.resize( 2*this->nodes().size() , 0.0 );
  geometry_params( geometry_ , topology_.points() , this->nodes().data() , this->nodes().size() , params0_.data() );

  real_t u0[2] = {0,0}; // dummy coordinates for ghost
  params_.create(u0);
  for (index_t k=0;k<this->nodes().size();k++)
  {
    params_.create( &params0_[2*k] );
    params_.set_entity( k+1 , topology_.points().entity(this->nodes()[k]) );

    U_[0] = params0_[2*k];
    U_[1] = params0_[2*k+1];

    // evaluate the coordinates for the orientation check
    geometry_->evaluate( U_ , X_ );
    for (coord_t d=0;d<3;d++)
      topology_.points()[ this->nodes()[k] ][d] = X_[d];
  }

  for (index_t k=0;k<this->geometry().points().nb();k++)
    this->geometry().points().set_entity( k , topology_.points().entity( u2v_[k] ) );

  this->extract_geometry(geometry_);

  GeometryOrientationChecker checker( topology_.points() , params_ , u2v_ , geometry_  );
  int s = checker.signof( this->geometry() );
  return s > 0;
}

template<typename type>
void
SurfaceCavity<type>::compute_coordinates()
{
  std::vector<real_t> U(2);
  std::vector<real_t> X(3);
  for (index_t k=0;k<this->nodes().size();k++)
  {
    if (this->nodes()[k] < topology_.points().nb_ghost()) continue;
    U[0] = topology_.points().u( this->nodes()[k] , 0 );
    U[1] = topology_.points().u( this->nodes()[k] , 1 );
    topology_.points().entity( this->nodes()[k] )->evaluate(U,X);
    for (coord_t d=0;d<3;d++)
      topology_.points()[ this->nodes()[k] ][d] = X[d];
  }
}

template class Primitive<Simplex>;
template class SurfaceCavity<Simplex>;

} // avro
