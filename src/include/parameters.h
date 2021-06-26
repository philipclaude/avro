#ifndef AVRO_API_PARAMETERS_H_
#define AVRO_API_PARAMETERS_H_

#include <iostream>
#include <map>
#include <string>

namespace avro
{

class Parameter {
public:

  virtual ~Parameter() {}
  virtual void print() const = 0;

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
    std::cout << std::boolalpha << "parameter[" << name_ << "]: " << value_ << std::endl;
  }

  void help() const {
    std::cout << "parameter[" << name_ << "]: " << description_ << "(default = " << default_ << ")" << std::endl;
  }

private:
  const T default_;
  T value_;
};

class ParameterSet
{
public:
  ParameterSet();

  void copy( const ParameterSet& params ) { assert(false); }

  template<typename T> void set_param( const std::string& name , const T& value ) {
    param_itr it = parameters_.find(name);
    if (it == parameters_.end()) {
      printf("parameter \"%s\" is not valid\n",name.c_str());
    }
    else {
      Parameter& param = *it->second.get();
      static_cast<ParameterType<T>&>(param).set(value);
    }
  }

  template<typename T> T get_param( const std::string& name ) const {
    param_const_itr it = parameters_.find(name);
    if (it == parameters_.end()) {
      printf("parameter \"%s\" is not valid\n",name.c_str());
      T x(0);
      return x;
    }
    const Parameter& param = *it->second.get();
    return static_cast<const ParameterType<T>&>(param);
  }

  void print() const {
    for (param_const_itr it = parameters_.begin(); it != parameters_.end(); ++it)
      it->second->print();
  }

protected:
  template <typename T> void register_parameter( const std::string& name , const T& def , const std::string& description ) {
    std::shared_ptr<Parameter> param = std::make_shared<ParameterType<T>>(name,def,description);
    parameters_.insert( {name,param} );
  }

private:
  typedef std::map<std::string,std::shared_ptr<Parameter>>::iterator param_itr;
  typedef std::map<std::string,std::shared_ptr<Parameter>>::const_iterator param_const_itr;
  std::map<std::string,std::shared_ptr<Parameter>> parameters_;
};

} // avro

#endif
