#include "common/error.h"

#include "avro_params.h"
#include "avro_types.h"

namespace avro
{

ParameterSet::ParameterSet() {
  set_defaults();
}

ParameterSet::ParameterSet( const ParameterSet& params ) {
  parameters_ = params.parameters();
}

void
ParameterSet::set_defaults() {
  register_parameter( "output redirect" , std::string() , "defines where avro output will be redirected (if supported)" );
  register_parameter( "directory" , std::string("./") , "defines where output files will be written" );
  register_parameter( "prefix" , std::string("avro_") , "output prefix" );
  register_parameter( "write mesh" , true , "controls whether the mesh is written upon completion of the algortihm" );
  register_parameter( "limit metric" , false , "controls whether the metric is limited from the current mesh" );
  register_parameter( "metric limiting factor" , 2.0 , "factor by which entries of the step matrix are limited when limiting the metric" );
  register_parameter( "cvt iterations" , index_t(10) , "controls the number of iterations performed when computing a Centroidal Voronoi Tesselation" );
  register_parameter( "sdot iterations" , index_t(10) , "controls the number of iterations performed when optimizing the weights of a power diagram" );
  register_parameter( "max parallel passes" , index_t(3) , "controls the maximum number of passes used in a parallel adaptation run" );
  register_parameter( "curved" , true , "controls whether an input geometry is curved (which will make for a faster adaptation if false)" );
}

template<typename T>
void
ParameterSet::set_param( const std::string& name , const T& value ) {
  param_itr it = parameters_.find(name);
  if (it == parameters_.end()) {
    printf("parameter \"%s\" is not valid\n",name.c_str());
    avro_assert_not_reached;
  }
  else {
    Parameter& param = *it->second.get();
    static_cast<ParameterType<T>&>(param).set(value);
  }
}

template<typename T>
T
ParameterSet::get_param( const std::string& name ) const {
  param_const_itr it = parameters_.find(name);
  if (it == parameters_.end()) {
    printf("parameter \"%s\" is not valid\n",name.c_str());
    T x(0);
    return x;
  }
  const Parameter& param = *it->second.get();
  return static_cast<const ParameterType<T>&>(param);
}

#define INSTANTIATE(T) template void ParameterSet::set_param( const std::string& , const T& ); \
                       template T ParameterSet::get_param( const std::string& ) const;

INSTANTIATE(std::string)
INSTANTIATE(bool)
INSTANTIATE(index_t)
INSTANTIATE(real_t)

#undef INSTANTIATE

} // avro
