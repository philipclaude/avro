//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/parameters.h"

#include <math.h>

namespace avro
{

void
AdaptationParameters::standard()
{
  // i/o parameters
  if (!stringParams_.has("prefix")) stringParams_["prefix"] = "mesh";
  if (!intParams_.has("adapt_iter")) intParams_["adapt_iter"] = 1;
  if (!stringParams_.has("directory")) stringParams_["directory"] = "./";
  if (!boolParams_.has("write_mesh")) boolParams_["write_mesh"] = true;
  if (!boolParams_.has("write_conformity")) boolParams_["write_conformity"] = true;
  if (!stringParams_.has("output_redirect")) stringParams_["output_redirect"] = std::string();
  if (!boolParams_.has("export_boundary")) stringParams_["export_boundary"] = true;

  // incoming mesh proeperties
  if (!boolParams_.has("curved")) boolParams_["curved"] = true;
  if (!boolParams_.has("has_interior_boundaries")) boolParams_["has_interior_boundaries"] = false;
  if (!boolParams_.has("has_uv")) boolParams_["has_uv"] = false;
  if (!boolParams_.has("partitioned")) boolParams_["partitioned"] = false;
  if (!boolParams_.has("balanced")) boolParams_["balacned"] = true;
  if (!intParams_.has("max_passes")) intParams_["max_passes"] = 10;

  // algorithm parameters
  if (!realParams_.has("insertion_volume_factor")) realParams_["insertion_volume_factor"] = sqrt(2.);
  if (!boolParams_.has("limit_insertion_length")) boolParams_["limit_insertion_length"] = true;
  if (!boolParams_.has("swapout")) boolParams_["swapout"] = true;
  if (!realParams_.has("lt_min")) realParams_["lt_min"] = sqrt(2.0);
  if (!realParams_.has("lt_max")) realParams_["lt_max"] = 2.0;
  if (!boolParams_.has("use_smoothing")) boolParams_["use_smoothing"] = true;
  if (!boolParams_.has("fefloa")) boolParams_["fefloa"] = false;
  if (!boolParams_.has("limit_metric")) boolParams_["limit_metric"] = false;
  if (!stringParams_.has("smoothing-equuation")) stringParams_["smoothing-equation"] = "avro"; // other option is bossen-heckbert
}

} // avro
