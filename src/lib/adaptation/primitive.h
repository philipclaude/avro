//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_PRIMITIVE_H_
#define avro_LIB_ADAPTATION_PRIMITIVE_H_

#include "adaptation/cavity.h"
#include "adaptation/geometry.h"

#include "common/types.h"

#include "geometry/egads/object.h"

#include "mesh/points.h"

#include "numerics/geometry.h"


#define ALL_WITH_EDGE 1

namespace avro
{

class Entity;
class Metric;
template<typename type> class MetricField;

template<typename type> class Topology;

template<typename type>
class Primitive : public Cavity<type>
{
public:
  Primitive( Topology<type>& _topology ) :
    Cavity<type>(_topology),
    u_(_topology.master().parameter()?(_topology.points().dim()):(_topology.points().dim()-1)),
    geometry_topology_(u_,_topology.master().parameter()?(_topology.number()):(_topology.number()-1)),
    geometry_cavity_( geometry_topology_ ),
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

  void convert_to_parameter( Entity* entity );
  void convert_to_physical( const std::vector<index_t>& N = std::vector<index_t>() );

  bool invalidatesTopology() const
  {
    //avro_implement; // is validBodies even needed???
    //if (!this->topology_.validBodies(*this)) return true;
    return false;
  }

  // cavity assignment and reset functions
  void set_cavity( const std::vector<index_t>& C ) { C_ = C; }
  void restart() { C_.clear(); }

  // geometry-related functions
  Entity* geometry( index_t e0 , index_t e1 );
  void extract_geometry( Entity* entity , const std::vector<index_t>& f=std::vector<index_t>() );

  real_t worst_quality_geometry( MetricField<type>& metric ) const;

  bool& delay() { return delay_; }
  bool& curved() { return curved_; }

  bool& debug() { return debug_; }

  Cavity<type>& geometry_cavity() { return geometry_cavity_; } // formerly gcavity_
  Topology<type>& geometry_topology() { return geometry_topology_; } // formerly G_
  Points& u() { return u_; }
  const std::vector<index_t>& S() const { return S_; }
  const std::vector<index_t>& u2v() const { return u2v_; }

protected:
  std::vector<index_t> C_; // saved cavity

  // data used for geometry checks
  Points u_;                         // points in parameter space
  Topology<type> geometry_topology_; // local geometry topology
  Cavity<type> geometry_cavity_;     // geometry cavity (for visibility check)
  std::map<index_t,index_t> v2u_;    // map from mesh points to u_
  std::vector<index_t> u2v_;         // map from u_ to mesh points
  std::vector<index_t> S_;           // list of elements in local_topology_ forming the geometry_cavity_
  bool delay_;                       // delay the application of the primitive operator to the topology
  bool curved_;                      // if the problem is curved so that we check visiblity in parameter space
  bool debug_;                       // whether we should debug
};

template<typename type>
class Collapse : public Primitive<type>
{
public:
  Collapse( Topology<type>& _topology );

  bool apply( const index_t p , const index_t q , bool delay=false );
  bool valid( const index_t p , const index_t q );
  bool visibleParameterSpace( index_t p , index_t q , Entity *g , bool edge=false );

  index_t& nb_parameter_tests() { return nb_parameter_tests_; }
  index_t& nb_parameter_rejections() { return nb_parameter_rejections_; }
  index_t& nb_invalid_geometry() { return nb_invalid_geometry_; }

  index_t& nb_accepted( const index_t d ) { return nb_accepted_[d]; }
  index_t& nb_rejected( const index_t d ) { return nb_rejected_[d]; }

  index_t& nb_rej_vis_Edge() { return nb_rej_vis_Edge_; }
  index_t& nb_rej_sgn_Edge() { return nb_rej_sgn_Edge_; }
  index_t& nb_rej_geo_Edge() { return nb_rej_geo_Edge_; }

private:
  index_t nb_parameter_tests_;
  index_t nb_parameter_rejections_;
  index_t nb_invalid_geometry_;

  std::vector<index_t> nb_accepted_;
  std::vector<index_t> nb_rejected_;

  index_t nb_rej_vis_Edge_;
  index_t nb_rej_sgn_Edge_;
  index_t nb_rej_geo_Edge_;
};

template<typename type>
class Insert : public Primitive<type>
{
public:
  Insert( Topology<type>& _topology );

  // application of the mesh operator to edge e0-e1 with insertion x
  // initial shell of elements on edge may be provided
  // option to delay the operator in the full mesh topology
  bool apply( const index_t e0 , const index_t e1 , real_t* x , real_t* u ,
              const std::vector<index_t>& shell=std::vector<index_t>() );

  // check if the insertion x is visible to the cavity boundary in
  // the parameter space
  bool visibleParameterSpace( real_t* x , real_t* u , Entity* e );

  bool enlarged() const { return enlarged_; }

  index_t& nb_parameter_tests() { return nb_parameter_tests_; }
  index_t& nb_parameter_rejections() { return nb_parameter_rejections_; }

private:
  index_t nb_parameter_tests_;
  index_t nb_parameter_rejections_;
  bool enlarged_;
  std::vector<index_t> disabled_;

  std::vector<index_t> elems_;
};

template<typename type>
class EdgeSwap : public Primitive<type>
{
public:
  EdgeSwap( Topology<type>& _topology );

  bool apply( const index_t p , const index_t e0 , const index_t e1 );
  bool valid( const index_t p , const index_t e0 , const index_t e1 );
  bool visibleParameterSpace(index_t p , index_t e0 , index_t e1 , Entity* g);

  index_t& nb_parameter_tests() { return nb_parameter_tests_; }
  index_t& nb_parameter_rejections() { return nb_parameter_rejections_; }
  index_t& nb_geometry_rejections( index_t d ) { return nb_geometry_rejections_[d]; }
  index_t& nb_interior() { return nb_interior_; }
  index_t& nb_invalid_geometry() { return nb_invalid_geometry_; }

private:
  index_t nb_parameter_tests_;
  index_t nb_parameter_rejections_;
  index_t nb_interior_;
  std::vector<index_t> nb_geometry_rejections_;
  index_t nb_invalid_geometry_;
};

template<typename type>
class FacetSwap : public Primitive<type>
{
public:
  FacetSwap( Topology<type>& _topology );

  bool apply( const index_t p , const index_t k0 , const index_t k1 );
  bool valid( const index_t p , const index_t k0 , const index_t k1 );
};

template<typename type>
class Smooth : public Primitive<type>
{
public:
  Smooth( Topology<type>& _topology);
  bool visibleParameterSpace( index_t p , real_t* x , real_t* params , Entity* ep );
  bool apply( const index_t p , MetricField<type>& metric , real_t Q0=-1 );

  index_t& nb_parameter_tests() { return nb_parameter_tests_; }
  index_t& nb_parameter_rejections() { return nb_parameter_rejections_; }
  real_t& delta() { return delta_; }
  real_t& delta_min() { return delta_min_; }
  real_t& delta_max() { return delta_max_; }
  real_t& objective() { return objective_; }

  index_t& nb_visibility_rejections() { return nb_visibility_rejections_; }
  index_t& nb_implied_metric_rejections() { return nb_implied_metric_rejections_; }
  index_t& nb_removed_rejections() { return nb_removed_rejections_; }
  index_t& nb_geometry() { return nb_geometry_; }
  index_t& nb_enlarged_rejections() { return nb_enlarged_rejections_; }
  index_t& nb_accepted() { return nb_accepted_; }
  index_t& Ntot() { return Ntot_; }
  index_t& nb_zero_valency() { return nb_zero_valency_; }
  index_t& nb_interpolated_outside() { return nb_interpolated_outside_; }

  void resetRejections()
  {
    nb_visibility_rejections_ = nb_parameter_rejections_ = nb_parameter_tests_ = 0;
    nb_implied_metric_rejections_ = nb_removed_rejections_ = nb_enlarged_rejections_ = 0;
    nb_geometry_ = 0;
    Ntot_ = 0;
    nb_zero_valency_ = 0;
    nb_interpolated_outside_ = 0;
  }

private:
  real_t delta_,delta_min_,delta_max_;
  index_t nb_parameter_tests_;
  index_t nb_parameter_rejections_;
  numerics::SymMatrixD<real_t> M0_;
  real_t objective_;

  index_t nb_visibility_rejections_;
  index_t nb_implied_metric_rejections_;
  index_t nb_removed_rejections_;
  index_t nb_enlarged_rejections_;
  index_t nb_geometry_;
  index_t nb_accepted_;
  index_t Ntot_; // average no. points attached to vertex used in smoothing computation
  index_t nb_zero_valency_;
  index_t nb_interpolated_outside_;
};

#define CROSS(a,b,c)      a[0] = (b[1]*c[2]) - (b[2]*c[1]);\
                          a[1] = (b[2]*c[0]) - (b[0]*c[2]);\
                          a[2] = (b[0]*c[1]) - (b[1]*c[0])
#define DOT(a,b)         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])

#define FLIP_MTYPE 1

class GeometryOrientationChecker
{
public:
  GeometryOrientationChecker( Points& v , Points& u,
                              const std::vector<index_t>& _u2v,
                              Entity* entity ) :
    v_(v),
    u_(u),
    u2v_(_u2v),
    entity_((EGADS::Object*)entity),
    normal_(TableLayout_Rectangular,3)
  {
    #if 0 // philip april 23
    avro_assert( u_.dim()==v_.dim()-1 );
    avro_assert( v_.dim()==3 );
    #endif
    avro_assert_msg( u2v_.size()==u_.nb() , "|u2v| = %lu, |u| = %lu", u2v_.size(),u_.nb() );

    // compute the normal to every vertex
    std::vector<real_t> N(v_.dim());
    for (index_t k=0;k<u_.nb();k++)
    {
      if (k<u_.nb_ghost())
      {
        std::fill(N.begin(),N.end(),0);
        normal_.add(N.data(),3);
        continue;
      }
      normal( k , N );
      normal_.add( N.data() , 3 );
    }
  }

  template<typename type> int signof( Topology<type>& topology , bool verbose=false )
  {
    //avro_assert_msg( topology.number()==u_.dim() , "topology number = %u, u_.dim = %u" , topology.number() , u_.dim() );

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

  template<typename type> bool hasPositiveVolumes( const Topology<type>& topology , int sign=1 ) const
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

  template<typename type> bool createsBadGeometry( const Topology<type>& topology ) const
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

  void normal( index_t k , std::vector<real_t>& N )
  {
    avro_assert( N.size()==3 );

    // evaluate the normal at the vertex using the stored parameter coordinates
    real_t result[18];
    EGADS_ENSURE_SUCCESS( EG_evaluate( *entity_->object() , u_[k] , result ) );
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
    EGADS_ENSURE_SUCCESS( EG_getInfo(*entity_->object(), &oclass, &mtype,&ref, &prev, &next) );
    N[0] *= mtype;
    N[1] *= mtype;
    N[2] *= mtype;
  }

  void normal( const real_t* v0 , const real_t* v1 , const real_t* v2 ,
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

private:
  Points& v_; // original points in physical space
  Points& u_; // points in parameter space
  const std::vector<index_t>& u2v_; // map from parameter to physical coordinates
  EGADS::Object* entity_; // the geometry entity in question

  Table<real_t> normal_; // vertex normals (using the geometry)
};

} // avro

#endif
