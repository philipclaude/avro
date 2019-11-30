#ifndef LUNA_LIB_MESH_FIELD_H_
#define LUNA_LIB_MESH_FIELD_H_

#include "common/table.h"
#include "common/json.h"
#include "common/types.h"

#include "master/master.h"
#include "master/polytope.h"
#include "master/simplex.h"

#include <map>
#include <string>
#include <vector>

namespace luna
{

template<typename T> class Table;
template<typename type> class Topology;

class Parameter;
class Coordinate;

class FieldHolder
{
private:
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
  T& eval( index_t elem , const Parameter& u ) const;
  T& eval( index_t elem , const Coordinate& x ) const;
  T& eval( const Coordinate& x ) const;

  T& operator()( index_t i , index_t j )
    { return data_[Table<index_t>::operator()(i,j)]; }

  const T& operator()( index_t i , index_t j ) const
    { return data_[Table<index_t>::operator()(i,j)]; }

  FieldType& type() { return type_; }
  const FieldType& type() const { return type_; }

  index_t nb_elem() const { return Table<index_t>::nb(); }
  index_t nb_data() const { return data_.size(); }

  const std::vector<T>& data() const { return data_; }

protected:
  FieldBase( FieldType type );
  std::vector<T> data_;

private:
  FieldType type_;
};

template<typename Shape,typename T> class Field;

template<typename T>
class Field<Simplex,T> : public FieldBase<T>
{
public:
  Field( const Topology<Simplex>& topology , coord_t order , FieldType type );
  void build();

  const Simplex& master() const { return master_; }

private:
  const Topology<Simplex>& topology_;
  const Simplex master_;
};

template<typename T>
class Field<Polytope,T> : public FieldBase<T>
{
public:
  Field( Topology<Polytope>& topology , coord_t order , FieldType type );
  void build();

private:
  const Topology<Polytope>& topology_;
  const Polytope master_;
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
    }

    template<typename type>
    type* get( const std::string& name )
    {
      FieldHolder* f = fields_[name].get();
      return static_cast<FieldClass<type>*>(f)->base();
    }

    template<typename type>
    const type* get( const std::string& name ) const
    {
      const FieldHolder* f = fields_.at(name).get();
      return static_cast<const FieldClass<type>*>(f)->base();
    }

    void remove( const std::string& name )
    {
      fields_.erase(name);
    }

    void fromJSON( const json& J );

private:
  std::map<std::string,std::shared_ptr<FieldHolder>> fields_;

  const TopologyBase& topology_;
};

} // luna

#endif
