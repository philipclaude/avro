//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_PRIMITIVE_H_
#define avro_LIB_ADAPTATION_PRIMITIVE_H_

#include "adaptation/cavity.h"

#include "avro_types.h"

#include "mesh/points.h"

namespace avro
{

class Entity;
class Metric;
template<typename type> class MetricField;

template<typename type> class Topology;

template<typename type>
class GeometryInspector
{
public:
  GeometryInspector( Points& v , Points& u, const std::vector<index_t>& _u2v, Entity* entity );

  void reset( Entity* entity );
  void compute_normals();

  int signof( Topology<type>& topology , bool verbose=false );
  bool positive_volumes( const Topology<type>& topology , int sign=1 ) const;
  bool invalidates_geometry( const Topology<type>& topology ) const;

  void normal( index_t k , std::vector<real_t>& N );
  void normal( const real_t* v0 , const real_t* v1 , const real_t* v2 ,
               std::vector<real_t>& N );

private:
  Points& v_; // original points in physical space
  Points& u_; // points in parameter space
  const std::vector<index_t>& u2v_; // map from parameter to physical coordinates
  Entity* entity_; // the geometry entity in question

  Table<real_t> normal_; // vertex normals (using the geometry)
};

template<typename type>
class Primitive : public Cavity<type>
{
public:
  Primitive( Topology<type>& _topology );

  void convert_to_parameter( Entity* entity );
  void convert_to_physical( const std::vector<index_t>& N = std::vector<index_t>() );

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

  Cavity<type>& geometry_cavity() { return geometry_cavity_; }
  Topology<type>& geometry_topology() { return geometry_topology_; }
  Points& u() { return u_; }
  const std::vector<index_t>& S() const { return S_; }
  const std::vector<index_t>& u2v() const { return u2v_; }

  std::vector<real_t>& cavity_quality() { return cavity_quality_; }

protected:
  std::vector<index_t> C_; // saved cavity

  // data used for geometry checks
  Points u_;                         // points in parameter space
  Topology<type> geometry_topology_; // local geometry topology
  Cavity<type> geometry_cavity_;     // geometry cavity (for visibility check)
  GeometryInspector<type> geometry_inspector_; // cheecker geometry orientation, volumes, invalidations
  std::map<index_t,index_t> v2u_;    // map from mesh points to u_
  std::vector<index_t> u2v_;         // map from u_ to mesh points
  std::vector<index_t> S_;           // list of elements in local_topology_ forming the geometry_cavity_
  bool delay_;                       // delay the application of the primitive operator to the topology
  bool curved_;                      // if the problem is curved so that we check visiblity in parameter space
  bool debug_;                       // whether we should debug

  std::vector<real_t> cavity_quality_;

};

template<typename type>
class Collapse : public Primitive<type>
{
public:
  Collapse( Topology<type>& _topology );

  bool apply( const index_t p , const index_t q , bool delay=false );
  bool valid( const index_t p , const index_t q );
  bool visible_geometry( index_t p , index_t q , Entity *g , bool edge=false );

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
              const std::vector<index_t>& shell=std::vector<index_t>() , int ns=-1 );

  // check if the insertion x is visible to the cavity boundary in
  // the parameter space
  bool visible_geometry( real_t* x , real_t* u , Entity* e , const std::vector<index_t>& edge=std::vector<index_t>() );

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
  bool visible_geometry(index_t p , index_t e0 , index_t e1 , Entity* g);

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
  bool visible_geometry( index_t p , real_t* x , real_t* params , Entity* ep );
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

  index_t& exponent() { return exponent_; }

  void set_inverse( std::vector<std::vector<index_t>>* inverse) { inverse_ = inverse; }

private:
  real_t delta_,delta_min_,delta_max_;
  index_t nb_parameter_tests_;
  index_t nb_parameter_rejections_;
  symd<real_t> M0_;
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

  index_t exponent_; // 1 for avro-style, 4 for bossen-heckbert style

  std::vector< std::vector<index_t> >* inverse_;
};

} // avro

#endif
