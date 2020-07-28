//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MESH_FIELD_H_
#define avro_LIB_MESH_FIELD_H_

#include "common/table.h"
#include "common/json.h"
#include "common/types.h"

#include "element/element.h"
#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/dof.h"

#include <map>
#include <string>
#include <vector>

namespace avro
{

template<typename T> class Table;
template<typename type> class Topology;

class Parameter;
class Coordinate;

class FieldHolder
{
public:
  virtual void evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const = 0;
  virtual ~FieldHolder() {}

  virtual real_t min( index_t rank ) const = 0;
  virtual real_t max( index_t rank ) const = 0;

  virtual std::string get_name( index_t j ) const = 0;
  virtual index_t nb_rank() const = 0;

protected:
  std::string name_;
};

enum FieldType
{
  CONTINUOUS,DISCONTINUOUS
};

template<typename T>
class FieldBase : public FieldHolder, public Table<index_t>
{
public:
  virtual ~FieldBase() {}

  T& eval( index_t elem , const Parameter& u ) const;
  T& eval( index_t elem , const Coordinate& x ) const;
  T& eval( const Coordinate& x ) const;

  virtual std::string get_name( index_t j ) const { return name_+std::to_string(j); }
  virtual index_t nb_rank() const { return 0; }

  void allocate( index_t n )
    { data_.allocate(n); }

  T& operator()( index_t i , index_t j )
    { return *data_[Table<index_t>::operator()(i,j)]; }

  const T& operator()( index_t i , index_t j ) const
    { return *data_[Table<index_t>::operator()(i,j)]; }

  FieldType& type() { return type_; }
  const FieldType& type() const { return type_; }

  index_t nb_elem() const { return Table<index_t>::nb(); }
  index_t nb_data() const { return data_.nb(); }

  T& value( index_t k ) { return *data_[k]; }
  const T& value( index_t k ) const { return *data_[k]; }

  const DOF<T>& dof() const { return data_; }

  real_t min( index_t rank ) const;
  real_t max( index_t rank ) const;

  virtual void evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const = 0;

  const ElementBase& element() const { return element_; }

protected:
  FieldBase( FieldType type , ElementBase& element , TableLayoutCategory category=TableLayout_Jagged );
  DOF<T> data_;

private:
  FieldType type_;
  ElementBase& element_;
};

template<typename type,typename T> class Field;

template<typename T>
class Field<Simplex,T> : public FieldBase<T>
{
public:
  Field( const Topology<Simplex>& topology , coord_t order , FieldType type );
  void build();

  const Simplex& element() const { return element_; }
  Simplex& element() { return element_; }

  template<typename Function>
  void evaluate( const Function& function );

  const Topology<Simplex>& topology() const { return topology_; }

  void evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const;

private:
  const Topology<Simplex>& topology_;
  Simplex element_;
};

template<typename T>
class Field<Polytope,T> : public FieldBase<T>
{
public:
  Field( Topology<Polytope>& topology , coord_t order , FieldType type );
  void build();

  void evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const;

private:
  const Topology<Polytope>& topology_;
  Polytope element_;
};

class TopologyBase;

template<typename derived_t>
class FieldClass : public FieldHolder
{
public:
  FieldClass( const std::string& name , derived_t& base ) :
    name_(name),
    base_(base)
  {}

  derived_t& base() { return base_; }
  const derived_t& base() const { return base_; }

  void evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const
    { base_.evaluate(rank,parents,alpha,result); }

  real_t min( index_t rank ) const
    { return base_.min(rank); }
  real_t max( index_t rank ) const
    { return base_.max(rank); }

  std::string get_name( index_t j ) const { return base_.get_name(j); }
  index_t nb_rank() const { return base_.nb_rank(); }

private:
  std::string name_;
  derived_t&  base_;
};

class Fields
{

public:
    Fields( const TopologyBase& topology );
    Fields( const json& J );

    bool has( const std::string& name ) const
    {
      if (fields_.find(name)==fields_.end())
        return false;
      return true;
    }

    template<typename type>
    void make( const std::string& name , std::shared_ptr<type>& f )
    {
      fields_.insert( { name , std::make_shared< FieldClass<type> >(name,*f.get()) } );
      identifiers_.insert( {name,pointer_string(f.get())} );
    }

    template<typename type>
    type* get( const std::string& name )
    {
      avro_assert_msg( fields_.find(name)!=fields_.end() , "could not find field %s" , name.c_str() );
      FieldHolder* f = fields_[name].get();
      return static_cast<FieldClass<type>*>(f)->base();
    }

    template<typename type>
    const type* get( const std::string& name ) const
    {
      avro_assert_msg( fields_.find(name)!=fields_.end() , "could not find field %s" , name.c_str() );
      const FieldHolder* f = fields_.at(name).get();
      return static_cast<const FieldClass<type>*>(f)->base();
    }

    const FieldHolder& operator[] ( const std::string& name ) const
    {
      avro_assert_msg( fields_.find(name)!=fields_.end() , "could not find field %s" , name.c_str() );
      return *fields_.at(name).get();
    }

    void remove( const std::string& name )
    {
      fields_.erase(name);
      identifiers_.erase(name);
    }

    index_t nb() const { return fields_.size(); }

    void from_json( const json& J );
    void get_names( std::vector<std::string>& names , std::vector<std::string>& ids ) const;

    std::string id2name( const std::string& id ) const;

    void print() const;

private:
  std::map<std::string,std::shared_ptr<FieldHolder>> fields_;
  std::map<std::string,std::string> identifiers_;

  const TopologyBase& topology_;
};

} // avro

#endif
