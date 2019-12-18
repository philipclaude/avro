// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_MESH_LOCAL_PRIMITIVE_H_
#define AVRO_MESH_LOCAL_PRIMITIVE_H_

#include "common/types.h"

#include "mesh/points.h"

#include "adaptation/cavity.h"

#define ALL_WITH_EDGE 1

#define PRIMITIVE_CHECK(X) if(!(X)) { this->logError( (#X) ); }

namespace luna
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
    u_(_topology.points().dim()-1),
    G_( u_ , _topology.number()-1 ),
    gcavity_( G_ ),
    delay_(false),
    curved_(true),
    debug_(true)
  {
    // do not allow the geometry cavity to enlarge
    gcavity_.setEnlarge(false);
    gcavity_.setIgnore(true);
  }

  bool invalidatesTopology() const
  {
    luna_implement; // is validBodies even needed???
    //if (!this->topology_.validBodies(*this)) return true;
    return false;
  }

  // cavity assignment and reset functions
  void setCavity( const std::vector<index_t>& C ) { C_ = C; }
  void restart() { C_.clear(); }

  // geometry-related functions
  Entity* geometry( index_t e0 , index_t e1 );
  void extractGeometry( Entity* entity , const std::vector<index_t>& f=std::vector<index_t>() );

  real_t worstQualityGeometry( MetricField<type>& metric ) const;

  bool& delay() { return delay_; }
  bool& curved() { return curved_; }

  bool& debug() { return debug_; }

  index_t nb_error() const { return errors_.size(); }
  void logError( const std::string& s )
  {
    errors_.push_back(s);
  }

  void printErrors() const
  {
    if (nb_error()==0) return;
    printf("%s:\n",this->name_.c_str());
    for (index_t k=0;k<nb_error();k++)
    {
      printf("error [%lu]: %s\n",k,errors_[k].c_str());
    }
  }

protected:
  std::vector<index_t> C_; // saved cavity

  // data used for geometry checks
  Points u_;                     // points in parameter space
  Topology<type> G_;               // geometry cavity topology
  Cavity<type> gcavity_;           // geometry cavity (for visibility check)
  std::map<index_t,index_t> v2u_;  // map from mesh points to u_
  std::vector<index_t> u2v_;       // map from u_ to mesh points
  std::vector<index_t> S_;         // list of elements in G_ (0 to G_.nb()-1)
  bool delay_;                     // delay the application of the primitive operator to the topology
  bool curved_;                    // if the problem is curved so that we check visiblity in parameter space
  bool debug_;                     // whether we should debug

  // string of error messages for debugging
  std::vector<std::string> errors_;

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
  index_t& nb_wake() { return nb_wake_; }
  index_t& nb_invalid_geometry() { return nb_invalid_geometry_; }

private:
  index_t nb_parameter_tests_;
  index_t nb_parameter_rejections_;
  index_t nb_wake_;
  std::vector<index_t> nb_geometry_rejections_;
  index_t nb_invalid_geometry_;
};

template<typename type>
class RidgeSwap : public Primitive<type>
{
public:
  RidgeSwap( Topology<type>& _topology );

  bool apply( const index_t p , const index_t t0 , const index_t t1 , const index_t t2 );
  bool valid( const index_t p , const index_t t0 , const index_t t1 , const index_t t2 );

  index_t& nb_geometry_rejections( index_t d ) { return nb_geometry_rejections_[d]; }

private:
  std::vector<index_t> nb_geometry_rejections_;
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

} // luna

#endif
