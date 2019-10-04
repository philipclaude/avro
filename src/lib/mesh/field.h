#ifndef URSA_LIB_MESH_FIELD_H_
#define URSA_LIB_MESH_FIELD_H_

#include "common/json.h"
#include "common/types.h"

#include "master/master.h"
#include "master/polytope.h"
#include "master/simplex.h"

#include <map>
#include <string>
#include <vector>

namespace ursa
{

template<typename T> class Data;
template<typename type> class Topology;

class Parameter;
class Coordinate;

class FieldHolder
{
private:
//  const TopologyHolder& topology_;
std::string name_;
};

enum FieldType
{
  CONTINUOUS,DISCONTINUOUS
};

template<typename T>
class FieldBase : public FieldHolder, public Data<index_t>
{
public:
  FieldBase();

  T& eval( index_t elem , const Parameter& u ) const;
  T& eval( index_t elem , const Coordinate& x ) const;
  T& eval( const Coordinate& x ) const;

  T& operator()( index_t i , index_t j )
    { return data_[Data<index_t>::operator()(i,j)]; }

  const T& operator()( index_t i , index_t j ) const
    { return data_[Data<index_t>::operator()(i,j)]; }

  FieldType& type() { return type_; }
  const FieldType& type() const { return type_; }

  index_t nb_elem() const { return Data<index_t>::nb(); }
  index_t nb_data() const { return data_.size(); }

protected:
  template<typename Shape_t,typename Master_t>
  void build( const Topology<Shape_t>& topology , const Master_t& master  );
  std::vector<T> data_;

private:
  FieldType type_;
};

template<typename Master_t,typename T> class Field;

template<typename Basis_t,typename T>
class Field<Simplex<Basis_t>,T> : public FieldBase<T>
{
  typedef Simplex<Basis_t>  Master_t;
  typedef Simplex<Lagrange> Shape_t;

public:
  Field( const Topology<Shape_t>& topology , coord_t order );
  void build();

private:
  const Topology<Shape_t>& topology_;
  const Master_t master_;
};

template<typename T>
class Field<Polytope,T> : public FieldBase<T>
{
  typedef Polytope Shape_t;
  typedef Polytope Master_t;

public:
  Field( Topology<Shape_t>& topology , coord_t order );

  void build();

private:
  const Topology<Shape_t>& topology_;
  const Master_t master_;
};


} // ursa

#endif
