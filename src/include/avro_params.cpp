#include "common/error.h"

#include "avro_params.h"
#include "avro_types.h"

namespace avro
{

Parameter::operator index_t() const {
  avro_assert( type_name_ == typeid(index_t).name());
  return static_cast<const ParameterType<index_t>&>(*this);
}

Parameter::operator real_t() const {
  avro_assert( type_name_ == typeid(real_t).name());
  return static_cast<const ParameterType<real_t>&>(*this);
}

Parameter::operator std::string() const {
  avro_assert( type_name_ == typeid(std::string).name());
  return static_cast<const ParameterType<std::string>&>(*this);
}

Parameter::operator bool() const {
  avro_assert( type_name_ == typeid(bool).name());
  return static_cast<const ParameterType<bool>&>(*this);
}

ParameterSet::ParameterSet() {
  set_defaults();
}

ParameterSet::ParameterSet( const ParameterSet& params ) {
  parameters_ = params.parameters();
}

void
ParameterSet::set_defaults() {
  register_parameter( "output redirect" , "" , "defines where avro output will be redirected (if supported)" );
  register_parameter( "directory" , "./" , "defines where output files will be written" );
  register_parameter( "prefix" , "avro_" , "output prefix" );
  register_parameter( "write mesh" , true , "controls whether the mesh is written upon completion of the algortihm" );
  register_parameter( "limit metric" , false , "controls whether the metric is limited from the current mesh" );
  register_parameter( "metric limiting factor" , 2.0 , "factor by which entries of the step matrix are limited when limiting the metric" );
  register_parameter( "cvt iterations" , index_t(10) , "controls the number of iterations performed when computing a Centroidal Voronoi Tesselation" );
  register_parameter( "sdot iterations" , index_t(10) , "controls the number of iterations performed when optimizing the weights of a power diagram" );
  register_parameter( "max parallel passes" , index_t(3) , "controls the maximum number of passes used in a parallel adaptation run" );
  register_parameter( "curved" , true , "controls whether an input geometry is curved (which will make for a faster adaptation if false)" );
  register_parameter( "adapt iter" , index_t(1) , "adaptation iteration" );
  register_parameter( "domain type" , "polytope" , "representation of the domain when calculating Voronoi diagrams options: polytope, simplex, mesh, sphere");
  register_parameter( "debug level" , index_t(1) , "how much debugging info should be generated (used for writing input files in parallel)" );
  register_parameter( "geometry" , "" , "name of the geometry used by avro, can either be an internal library name or a file name (supported by EGADS)" );
}

template<typename T>
void
ParameterSet::set( const std::string& name , const T& value ) {
  param_itr it = parameters_.find(name);
  if (it == parameters_.end()) {
    printf("parameter \"%s\" is not valid\n",name.c_str());
    avro_assert_not_reached;
  }
  else {
    Parameter& param = *it->second.get();
    avro_assert( param.type_name() == typeid(T).name() );
    static_cast<ParameterType<T>&>(param).set(value);
  }
}

const Parameter&
ParameterSet::operator[] ( const std::string& name ) const {
  param_const_itr it = parameters_.find(name);
  if (it == parameters_.end()) {
    printf("parameter \"%s\" is not valid\n",name.c_str());
    avro_assert_not_reached;
  }
  return *it->second.get();
}

#define INSTANTIATE(T) template void ParameterSet::set( const std::string& , const T& );

INSTANTIATE(std::string)
INSTANTIATE(bool)
INSTANTIATE(index_t)
INSTANTIATE(real_t)

#undef INSTANTIATE

} // avro
