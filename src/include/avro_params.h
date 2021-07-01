#ifndef AVRO_API_PARAMETERS_H_
#define AVRO_API_PARAMETERS_H_

#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "avro_types.h"

namespace avro
{

template<typename T> class ParameterType;

class Parameter {
public:

  virtual ~Parameter() {}
  virtual void print() const = 0;
  virtual void help() const = 0;

  const std::string& type_name() const { return type_name_; }

protected:
  Parameter( const std::string& name , const std::string& description , const std::string& type_name ) :
    name_(name),
    description_(description),
    type_name_(type_name)
  {}

public:
  operator index_t() const;
  operator real_t() const;
  operator std::string() const;
  operator bool() const;

protected:
  std::string name_;
  std::string description_;
  std::string type_name_;
};

template<typename T>
class ParameterType : public Parameter {

public:
  ParameterType( const std::string& name , const T& def , const std::string& description ) :
    Parameter(name,description,typeid(T).name()),
    default_(def),
    value_(def)
  {}

  void set( const T& value ) { value_ = value; }

  operator T()       { return value_; }
  operator T() const { return value_; }

  void print() const {
    std::cout << std::boolalpha << "\tparameter[\"" << name_ << "\"]: " << value_ << std::endl;
  }

  void help() const {
    std::cout << std::boolalpha << "\tparameter[\"" << name_ << "\"]: " << description_ << "(default = " << default_ << ")" << std::endl;
  }

private:
  const T default_;
  T value_;
};

class ParameterSet
{
public:
  ParameterSet();
  ParameterSet( const ParameterSet& params );
  virtual ~ParameterSet() {}

  virtual void set_defaults();

  template<typename T> void set( const std::string& name , const T& value );

  void set( const std::string& name , const char* value ) {
    set(name,std::string(value));
  }

  void print() const {
    printf("ParameterSet:\n");
    for (param_const_itr it = parameters_.begin(); it != parameters_.end(); ++it)
      it->second->print();
  }

  void help() const {
    printf("ParameterSet:\n");
    for (param_const_itr it = parameters_.begin(); it != parameters_.end(); ++it)
      it->second->help();
  }

  const Parameter& operator[] ( const std::string& name ) const;


protected:
  template <typename T> void register_parameter( const std::string& name , const T& def , const std::string& description ) {
    if (parameters_.find(name) != parameters_.end()) return;
    std::shared_ptr<Parameter> param = std::make_shared<ParameterType<T>>(name,def,description);
    parameters_.insert( {name,param} );
  }

  void register_parameter( const std::string& name , const char* def , const std::string& description ) {
    register_parameter(name,std::string(def),description);
  }

  const std::map<std::string,std::shared_ptr<Parameter>>& parameters() const { return parameters_; }

private:
  typedef std::map<std::string,std::shared_ptr<Parameter>>::iterator param_itr;
  typedef std::map<std::string,std::shared_ptr<Parameter>>::const_iterator param_const_itr;
  std::map<std::string,std::shared_ptr<Parameter>> parameters_;
};

// these should be avoided in parameters because they are confusing
template<> void ParameterSet::set( const std::string& name , const int& value ) = delete;
template<> void ParameterSet::register_parameter( const std::string& , const int& , const std::string& description ) = delete;
template<> void ParameterSet::set( const std::string& name , const unsigned int& value ) = delete;
template<> void ParameterSet::register_parameter( const std::string& , const unsigned int& , const std::string& description ) = delete;
template<> void ParameterSet::set( const std::string& name , const short& value ) = delete;
template<> void ParameterSet::register_parameter( const std::string& , const short& , const std::string& description ) = delete;
template<> void ParameterSet::set( const std::string& name , const unsigned short& value ) = delete;
template<> void ParameterSet::register_parameter( const std::string& , const unsigned short& , const std::string& description ) = delete;

} // avro

#endif
