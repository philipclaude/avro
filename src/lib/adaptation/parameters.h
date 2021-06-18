//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_PARAMETERS_H_
#define avro_LIB_ADAPTATION_PARAMETERS_H_

#include "common/parameters.h"
#include "common/types.h"

#include <map>
#include <string>
#include <vector>

namespace avro
{

class AdaptationParameters : public Parameters
{
public:
  AdaptationParameters()
  {
    names_ = {"directory","prefix","adapt_iter","output_redirect","write_mesh","write_conformity",
    "has_interior_boundaries","curved","has_uv","partitioned", "balanced","max_passes",
    "insertion_volume_factor","limit_insertion_length","swapout","lt_min",
    "lt_max","use_smoothing","fefloa","limit_metric"};
  }

  // i/o parameters
  std::string& directory() { return stringParams_["directory"]; }
  std::string& prefix() { return stringParams_["prefix"]; }
  int& adapt_iter() { return intParams_["adapt_iter"]; }
  std::string& output_redirect() { return stringParams_["output_redirect"]; }
  bool& write_mesh() { return boolParams_["write_mesh"]; }
  bool& write_conformity() { return boolParams_["write_conformity"]; }
  bool& export_boundary() { return boolParams_["export_boundary"]; }

  // properties of the incoming mesh
  bool& curved() { return boolParams_["curved"]; }
  bool& has_interior_boundaries() { return boolParams_["has_interior_boundaries"]; }
  bool& has_uv() { return boolParams_["has_uv"]; }
  bool& partitioned() { return boolParams_["partitioned"]; }
  bool& balanced() { return boolParams_["balanced"]; }
  int& max_passes() { return intParams_["max_passes"]; }

  // algorithm parameters
  real_t& insertion_volume_factor() { return realParams_["insertion_volume_factor"]; }
  bool& limit_insertion_length() { return boolParams_["limit_insertion_length"]; }
  bool& swapout() { return boolParams_["swapout"]; }
  real_t& lt_min() { return realParams_["lt_min"]; }
  real_t& lt_max() { return realParams_["lt_max"]; }
  bool& use_smoothing() { return boolParams_["use_smoothing"]; }
  bool& fefloa() { return boolParams_["fefloa"]; }
  bool& limit_metric() { return boolParams_["limit_metric"]; }
  int& smoothing_exponent() { return intParams_["smoothing-exponent"]; }

  std::string& parallel_method() { return stringParams_["parallel method"]; }
  int& elems_per_processor() { return intParams_["elems_per_processor"]; }

  // set the default
  void standard();
};

} // avro

#endif
