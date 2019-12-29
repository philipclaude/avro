#ifndef luma_LIB_ADAPTATION_PARAMETERS_H_
#define luma_LIB_ADAPTATION_PARAMETERS_H_

#include "common/types.h"

#include <map>
#include <string>
#include <vector>

namespace luma
{

template<typename type>
class ParameterContainer : public std::map<std::string,type>
{
public:
  bool has( const std::string& name ) const
  {
    typename std::map<std::string,type>::const_iterator it =
                                      std::map<std::string,type>::find(name);
    if (it==std::map<std::string,type>::end()) return false;
    return true;
  }
};

class Parameters
{
protected:
  ParameterContainer<std::string> stringParams_;
  ParameterContainer<int>         intParams_;
  ParameterContainer<real_t>      realParams_;
  ParameterContainer<bool>        boolParams_;

  std::vector<std::string> names_;
};

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

  std::string& output_redirect() { return stringParams_["output_redirect"]; }

  bool& has_uv() { return boolParams_["has_uv"]; }

  bool& debug() { return boolParams_["debug"]; }

  void standard();

  void print()
  {
    printf("Parameters:\n");
    for (index_t k=0;k<names_.size();k++)
    {
      if (realParams_.has(names_[k]))
        printf("\t%s = %g\n",names_[k].c_str(),realParams_[names_[k]]);
      if (intParams_.has(names_[k]))
        printf("\t%s = %d\n",names_[k].c_str(),intParams_[names_[k]]);
      if (boolParams_.has(names_[k]))
        printf("\t%s = %s\n",names_[k].c_str(),(boolParams_[names_[k]])?"true":"false");
      if (stringParams_.has(names_[k]))
        printf("\t%s = %s\n",names_[k].c_str(),stringParams_[names_[k]].c_str());
    }
  }
};

} // luma

#endif
