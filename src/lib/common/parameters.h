#ifndef AVRO_LIB_COMMON_PARAMETERS_H_
#define AVRO_LIB_COMMON_PARAMETERS_H_

#include "common/types.h"

#include <map>
#include <string>
#include <vector>

namespace avro
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
public:
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

protected:
  ParameterContainer<std::string> stringParams_;
  ParameterContainer<int>         intParams_;
  ParameterContainer<real_t>      realParams_;
  ParameterContainer<bool>        boolParams_;

  std::vector<std::string> names_;
};


} // avro

#endif
