//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/geometry.h"

#include "common/tools.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "library/metric.h"

#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include "numerics/geometry.h"

#include <egads.h>

#include <set>

namespace avro
{

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

template<typename type>
Primitive<type>::Primitive( Topology<type>& _topology ) :
  Cavity<type>(_topology),
  u_(_topology.master().parameter()?(_topology.points().dim()):(_topology.points().dim()-1)),
  geometry_topology_(u_,_topology.master().parameter()?(_topology.number()):(_topology.number()-1)),
  geometry_cavity_( geometry_topology_ ),
  geometry_inspector_( _topology.points() , u_ , u2v_ , nullptr ),
  delay_(false),
  curved_(true),
  debug_(true)
{
  // do not allow the geometry cavity to enlarge
  geometry_cavity_.set_enlarge(false);
  geometry_cavity_.set_ignore(true);
  geometry_cavity_.master().set_parameter( _topology.master().parameter() );
  geometry_topology_.master().set_parameter( _topology.master().parameter() );
}

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

template<typename type>
GeometryInspector<type>::GeometryInspector( Points& v , Points& u, const std::vector<index_t>& _u2v, Entity* entity ) :
  v_(v),
  u_(u),
  u2v_(_u2v),
  entity_((EGADS::Object*)entity),
  normal_(TableLayout_Rectangular,v.dim())
{}

template<typename type>
void
GeometryInspector<type>::reset( Entity* entity )
{
  avro_assert_msg( u2v_.size()==u_.nb() , "|u2v| = %lu, |u| = %lu", u2v_.size(),u_.nb() );
  entity_ = entity;
  compute_normals();
}

template<typename type>
void
GeometryInspector<type>::compute_normals()
{
  normal_.clear();

  // compute the normal to every vertex
  std::vector<real_t> N(v_.dim());
  for (index_t k=0;k<u_.nb();k++)
  {
    if (k<u_.nb_ghost())
    {
      std::fill(N.begin(),N.end(),0);
      normal_.add(N.data(),N.size());
      continue;
    }
    normal( k , N );
    normal_.add( N.data() , N.size() );
  }
}

template<typename type>
int
GeometryInspector<type>::signof( Topology<type>& topology , bool verbose )
{
  std::vector<int> sign( topology.nb() );

  // loop through all the cells (facets) in the topology
  std::vector<real_t> N(v_.dim());
  int s = 1;
  index_t counted = 0;
  for (index_t k=0;k<topology.nb();k++)
  {
    if (topology.ghost(k)) continue;
    counted++;

    const real_t* v0 = v_[ u2v_.at(topology(k,0)) ];
    const real_t* v1 = v_[ u2v_.at(topology(k,1)) ];
    const real_t* v2 = v_[ u2v_.at(topology(k,2)) ];

    if (verbose)
    {
      printf("checking facet (%lu,%lu,%lu):\n",
              u2v_.at(topology(k,0)),u2v_.at(topology(k,1)),u2v_.at(topology(k,2)));
      printf("v0 = (%g,%g,%g)\n",v0[0],v0[1],v0[2]);
      printf("v1 = (%g,%g,%g)\n",v1[0],v1[1],v1[2]);
      printf("v2 = (%g,%g,%g)\n",v2[0],v2[1],v2[2]);
    }

    // compute the normal of this facet in physical space
    normal( v0 , v1 , v2 , N );

    sign[k] = 1;
    for (index_t j=0;j<3;j++)
    {
      avro_assert( topology(k,j)>=topology.points().nb_ghost() );
      const real_t* n = normal_(topology(k,j));
      real_t dp = DOT(N,n);
      if (verbose)
      {
        printf("vertex %lu dot product = %g\n",topology(k,j),dp);
        std::vector<real_t> n0(n,n+3);
        print_inline(N,"nt");
        print_inline(n0,"ng");
      }
      if (dp<0.1)
      {
        sign[k] = -1;
        s = -1;
      }
    }
  }
  avro_assert_msg( counted>0 , "topology.nb = %lu" , topology.nb() );
  return s;
}

template<typename type>
bool
GeometryInspector<type>::positive_volumes( const Topology<type>& topology , int sign ) const
{
  std::vector<real_t> vg(topology.nb());
  topology.get_volumes(vg);
  for (index_t k=0;k<vg.size();k++)
  {
    if (topology.ghost(k)) continue;
    if (sign*vg[k]<=0.0)
    {
      print_inline(vg,"vols");
      return false;
    }
  }
  return true;
}

template<typename type>
bool
GeometryInspector<type>::invalidates_geometry( const Topology<type>& topology ) const
{
  for (index_t k=0;k<topology.nb();k++)
  {
    if (topology.ghost(k)) continue;

    // check for weird geometry configurations that cannot be recovered from
    Entity* e0 = topology.points().entity( topology(k,0) );
    Entity* e1 = topology.points().entity( topology(k,1) );
    Entity* e2 = topology.points().entity( topology(k,2) );

    if (e0==NULL || e1==NULL || e2==NULL) continue;

    Entity* face = e0->intersect(e1,e2,true); // only check, don't error!
    if (face==NULL)
    {
      return true;
    }
  }
  return false;
}

template<typename type>
void
GeometryInspector<type>::normal( index_t k , std::vector<real_t>& N )
{
  avro_assert( N.size()==3 );

  EGADS::Object* entity = (EGADS::Object*)entity_;

  // evaluate the normal at the vertex using the stored parameter coordinates
  real_t result[18];
  EGADS_ENSURE_SUCCESS( EG_evaluate( *entity->object() , u_[k] , result ) );
  real_t dx_du[3],dx_dv[3];
  dx_du[0] = result[3];
  dx_du[1] = result[4];
  dx_du[2] = result[5];
  dx_dv[0] = result[6];
  dx_dv[1] = result[7];
  dx_dv[2] = result[8];
  CROSS(N, dx_du, dx_dv);
  numerics::normalize( N.data() , N.size() );

  // flip the normal depending on the type of face (SFORWARD or SREVERSE)
  int oclass,mtype;
  ego ref,prev,next;
  EGADS_ENSURE_SUCCESS( EG_getInfo(*entity->object(), &oclass, &mtype,&ref, &prev, &next) );
  N[0] *= mtype;
  N[1] *= mtype;
  N[2] *= mtype;
}

template<typename type>
void
GeometryInspector<type>::normal( const real_t* v0 , const real_t* v1 , const real_t* v2 ,
             std::vector<real_t>& N )
{
  // evaluate the normal of the oriented triangle
  avro_assert( N.size()==3 );
  real_t X01[3],X02[3];
  for (coord_t d=0;d<3;d++)
  {
    X01[d] = v1[d] -v0[d];
    X02[d] = v2[d] -v0[d];
  }
  CROSS(N,X01,X02);
  real_t ln = std::sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
  avro_assert(ln!=0.0);
  numerics::normalize( N.data() , N.size() );
}

template class Primitive<Simplex>;
template class GeometryInspector<Simplex>;

} // avro
