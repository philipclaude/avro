#ifndef AVRO_API_PARAMETERS_H_
#define AVRO_API_PARAMETERS_H_

#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace avro
{

class Parameter {
public:

  virtual ~Parameter() {}
  virtual void print() const = 0;
  virtual void help() const = 0;

protected:
  Parameter( const std::string& name , const std::string& description ) :
    name_(name),
    description_(description)
  {}

protected:
  std::string name_;
  std::string description_;
};

template<typename T>
class ParameterType : public Parameter {

public:
  ParameterType( const std::string& name , const T& def , const std::string& description ) :
    Parameter(name,description),
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
  template<typename T> void set_param( const std::string& name , const T& value );
  template<typename T> T get_param( const std::string& name ) const;

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

protected:
  template <typename T> void register_parameter( const std::string& name , const T& def , const std::string& description ) {
    std::shared_ptr<Parameter> param = std::make_shared<ParameterType<T>>(name,def,description);
    parameters_.insert( {name,param} );
  }

  template <typename T> void register_parameter( const std::string& name , const char* def , const std::string& description ) {
    register_parameter(name,std::string(def),description);
  }

  const std::map<std::string,std::shared_ptr<Parameter>>& parameters() const { return parameters_; }

private:
  typedef std::map<std::string,std::shared_ptr<Parameter>>::iterator param_itr;
  typedef std::map<std::string,std::shared_ptr<Parameter>>::const_iterator param_const_itr;
  std::map<std::string,std::shared_ptr<Parameter>> parameters_;
};

// integers should be avoided in parameters because they are confusing
template<> void ParameterSet::set_param( const std::string& name , const int& value ) = delete;
template<> int  ParameterSet::get_param( const std::string& name ) const = delete;

} // avro

#endif
