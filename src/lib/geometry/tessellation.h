#ifndef AVRO_LIB_GEOMETRY_TESSELLATION_H_
#define AVRO_LIB_GEOMETRY_TESSELLATION_H_

#include "common/parameters.h"

#include "mesh/mesh.h"

namespace avro
{

class Body;
class Model;

class TessellationParameters : public Parameters
{
public:
  TessellationParameters()
  {
    names_ = {"type","min-size","min-length","min-angle"};
  }

  std::string& type() { return stringParams_["type"]; }

  void standard()
  {
    if (!stringParams_.has("type")) stringParams_["type"] = "simplex";
    if (!realParams_.has("min-size")) realParams_["min-size"] = 0.25;
    if (!realParams_.has("min-length")) realParams_["min-length"] = 0.01;
    if (!realParams_.has("min-angle")) realParams_["min-angle"] = 30;
  }
};

class BodyTessellation : public Mesh
{
public:
  BodyTessellation( Body& body , TessellationParameters& params );

  const real_t* internal_point( const index_t k ) const { return internal_points_[k]; }
  index_t nb_internal() const { return internal_points_.nb(); }


  void make_internal_points();
  real_t volume();

  TessellationParameters& parameters() { return params_; }

private:
  Body& body_;
  TessellationParameters& params_;

  Points internal_points_;
};

class ModelTessellation : public Mesh
{
public:
  enum ModelTessellationType { Simplicial , QuadrilateralType };

  ModelTessellation(Model& model , TessellationParameters& params );

  void copy_mesh( const BodyTessellation& body_tess );
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
