//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_GEOMETRY_TESSELLATION_H_
#define AVRO_LIB_GEOMETRY_TESSELLATION_H_

#include "avro_params.h"

#include "mesh/mesh.h"

namespace avro
{

class Body;
class Model;

class TessellationParameters : public ParameterSet
{
public:
  TessellationParameters() {
    register_parameter( "type" , "simplex" , "element type for tessellation" );
    register_parameter( "min size" , 0.25 , "minimum size" );
    register_parameter( "min length" , 0.01 , "minimum length" );
    register_parameter( "min angle" , 30.0 , "minimum angle in degrees" );
  }
};

class BodyTessellation : public Mesh
{
public:
  BodyTessellation( Points& model_points , Body& body , TessellationParameters& params );

  const real_t* internal_point( const index_t k ) const { return internal_points_[k]; }
  index_t nb_internal() const { return internal_points_.nb(); }


  void make_internal_points();
  real_t volume();

  TessellationParameters& parameters() { return params_; }

  Points& model_points() { return model_points_; }

private:
  Body& body_;
  TessellationParameters& params_;
  Points& model_points_;

  Points internal_points_;
};

class ModelTessellation : public Mesh
{
public:
  enum ModelTessellationType { Simplicial , QuadrilateralType };

  ModelTessellation(Model& model , TessellationParameters& params );

  void copy_mesh( BodyTessellation& body_tess );
  Model& model() const;

  void get_body_internal_points( const BodyTessellation& body_tess );
  const real_t* internal_vertex( const index_t k ) const { return internal_points_[k]; }
  index_t nb_internal() const { return internal_points_.nb(); }

  index_t nb_bodies() const { return nb_topologies(); }
  //Topology<Simplex>& body( const index_t k ) const { return *root_->child(k); }
  bool interior( const index_t k ) const { return interior_[k]; }

  real_t volume() const { return volume_; }

private:
  Points internal_points_;
  Model& model_;
  real_t volume_;

  std::vector<bool> interior_;
};

} // avro


#endif
