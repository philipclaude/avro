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
    names_ = {"parallel","nb_iter","algorithm",
      "directory","adapt_iter","prepared","curved","prefix","boundarySubdirectory",
      "write_json","write_meshb","insertion_volume_factor",
      "write_mesh","write_conformity","has_interior_boundaries",
      "limit_insertion_length","swapout","lt_min","lt_max","use_smoothing","fefloa",
      "output_redirect","has_uv"};
  }

  bool& parallel() { return boolParams_["parallel"]; }
  int& nb_iter() { return intParams_["nb_iter"]; }
  std::string& algorithm() { return stringParams_["algorithm"]; }
  std::string& directory() { return stringParams_["directory"]; }
  int& adapt_iter() { return intParams_["adapt_iter"]; }
  bool& prepared() { return boolParams_["prepared"]; }
  bool& curved() { return boolParams_["curved"]; }

  std::string& prefix() { return stringParams_["prefix"]; }
  std::string& boundarySubdirectory()
    { return stringParams_["boundarySubdirectory"]; }
  bool& write_json() { return boolParams_["write_json"]; }
  bool& write_meshb() { return boolParams_["write_meshb"]; }

  real_t& insertion_volume_factor() { return realParams_["insertion_volume_factor"]; }

  bool& write_mesh() { return boolParams_["write_mesh"]; }
  bool& write_conformity() { return boolParams_["write_conformity"]; }

  bool& has_interior_boundaries() { return boolParams_["has_interior_boundaries"]; }

  bool& limit_insertion_length() { return boolParams_["limit_insertion_length"]; }
  bool& swapout() { return boolParams_["swapout"]; }
  real_t& lt_min() { return realParams_["lt_min"]; }
  real_t& lt_max() { return realParams_["lt_max"]; }
  bool& use_smoothing() { return boolParams_["use_smoothing"]; }
  bool& fefloa() { return boolParams_["fefloa"]; }

  bool& limit_metric() { return boolParams_["limit_metric"]; }

  std::string& output_redirect() { return stringParams_["output_redirect"]; }

  bool& has_uv() { return boolParams_["has_uv"]; }

  bool& debug() { return boolParams_["debug"]; }

  void standard();
};

} // avro

#endif
