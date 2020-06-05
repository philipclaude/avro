//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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
  geometry_topology_.clear();
  v2u_.clear();
  u2v_.clear();
  geometry_cavity_.clear();
  S_.clear();
  geometry_topology_.neighbours().forceCompute(); // forces the neighbours to be recomputed
  geometry_topology_.set_closed(false); // forces the re-closing of the mesh

  this->compute_geometry( e , geometry_topology_ , v2u_ , u2v_ );
  if (geometry_topology_.nb()==0) return;

  geometry_topology_.inverse().build();
  if (f.size()==0)
  {
    // a size of zero is a code for extracting non-ghost elements
    for (index_t k=0;k<geometry_topology_.nb();k++)
    {
      if (geometry_topology_.ghost(k)) continue;
      S_.push_back(k);
    }
  }
  else if (f.size()==1)
  {
    // requested a vertex v, look up the u value
    index_t u = v2u_.at(f[0]);
    geometry_topology_.inverse().ball(u,S_);
  }
  else if (f.size()==2)
  {
    // requested an edge, lookup non-ghost elements
    index_t u0 = v2u_.at(f[0]);
    index_t u1 = v2u_.at(f[1]);
    geometry_topology_.inverse().shell(u0,u1,S_);
    avro_assert( S_.size()==2 );
  }
  else
  {
    print_inline(f,"unsupported facet: ");
    avro_assert_not_reached;
  }
}

template<typename type>
void
Primitive<type>::convert_to_parameter( Entity* entity )
{
  // convert the parameter coordinates
  avro_assert( this->topology_.master().parameter() );
  for (index_t k=0;k<this->geometry_cavity_.points().nb();k++)
  {
    if (k < this->geometry_cavity_.points().nb_ghost()) continue;
    geometry_params( entity , this->geometry_cavity_.points() , &k , 1 , this->geometry_cavity_.points()[k] );
  }
}

template<typename type>
void
Primitive<type>::convert_to_physical( const std::vector<index_t>& N )
{
  if (N.size()==0)
  {
    for (index_t k=0;k<geometry_cavity_.points().nb();k++)
    {
      if (k < geometry_cavity_.points().nb_ghost()) continue;
      index_t m = this->u2v_.at(k);
      std::vector<real_t> U( this->points_.u(m) , this->points_.u(m)+2 );
      std::vector<real_t> X(3);
      this->points_.entity(m)->evaluate(U,X);
      for (coord_t d=0;d<3;d++)
        this->points_[m][d] = X[d];
    }
  }
  else
  {
    for (index_t k=0;k<N.size();k++)
    {
      avro_assert( N[k] >= this->topology_.points().nb_ghost() );
      std::vector<real_t> U( this->points_.u(N[k]) , this->points_.u(N[k])+2 );
      std::vector<real_t> X(3);
      this->points_.entity(N[k])->evaluate(U,X);
      for (coord_t d=0;d<3;d++)
        this->points_[N[k]][d] = X[d];
    }
  }
}

template class Primitive<Simplex>;

} // avro
