#include "adaptation/parameters.h"

#include <math.h>

namespace avro
{

void
AdaptationParameters::standard()
{
  if (!boolParams_.has("parallel")) boolParams_["parallel"] = true;
  if (!intParams_.has("nb_iter")) intParams_["nb_iter"] = 1;
  if (!intParams_.has("adapt_iter")) intParams_["adapt_iter"] = 1;
  if (!stringParams_.has("algorithm")) stringParams_["algorithm"] = "patch";
  if (!stringParams_.has("directory")) stringParams_["directory"] = "./";
  if (!boolParams_.has("prepared")) boolParams_["prepared"] = false;
  if (!boolParams_.has("curved")) boolParams_["curved"] = true;
  if (!stringParams_.has("boundarySubdirectory")) stringParams_["boundarySubdirectory"] = "./";
  if (!boolParams_.has("write_json")) boolParams_["write_json"] = true;
  if (!boolParams_.has("write_meshb")) boolParams_["write_meshb"] = true;
  if (!stringParams_.has("prefix")) stringParams_["prefix"] = "mesh";
  if (!realParams_.has("insertion_volume_factor")) realParams_["insertion_volume_factor"] = sqrt(2.);
  if (!boolParams_.has("write_mesh")) boolParams_["write_mesh"] = true;
  if (!boolParams_.has("write_conformity")) boolParams_["write_conformity"] = true;
  if (!boolParams_.has("has_interior_boundaries")) boolParams_["has_interior_boundaries"] = false;
  if (!boolParams_.has("limit_insertion_length")) boolParams_["limit_insertion_length"] = true;
  if (!boolParams_.has("swapout")) boolParams_["swapout"] = true;
  if (!realParams_.has("lt_min")) realParams_["lt_min"] = sqrt(2.0);
  if (!realParams_.has("lt_max")) realParams_["lt_max"] = 2.0;
  if (!boolParams_.has("use_smoothing")) boolParams_["use_smoothing"] = true;
  if (!boolParams_.has("fefloa")) boolParams_["fefloa"] = false;
  if (!stringParams_.has("output_redirect")) stringParams_["output_redirect"] = std::string();
  if (!boolParams_.has("has_uv")) boolParams_["has_uv"] = false;
  if (!boolParams_.has("debug")) boolParams_["debug"] = true;
  if (!boolParams_.has("limit_metric")) boolParams_["limit_metric"] = false;
}

} // avro
