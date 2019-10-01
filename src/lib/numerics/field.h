#ifndef URSA_LIB_NUMERICS_FIELD_H_
#define URSA_LIB_NUMERICS_FIELD_H_

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

template<typename T>
class FieldBase : public FieldHolder
{
public:
  T& eval( index_t elem , const Parameter& u ) const;
  T& eval( index_t elem , const Coordinate& x ) const;
  T& eval( const Coordinate& x ) const;

  void add( const T& x );
  void remove( index_t k );

  index_t nb() const { return data_.nb(); }

protected:
  Data<T> data_;
};

template<typename Master_t,typename T> class Field;

template<typename Basis_t,typename T>
class Field<Simplex<Basis_t>,T> : public FieldBase<T>
{
  typedef Simplex<Basis_t>  Master_t;
  typedef Simplex<Lagrange> Shape_t;

public:
  Field( const Topology<Shape_t>& topology , coord_t order );

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

private:
  const Topology<Shape_t>& topology_;
  const Master_t master_;
};

class Fields
{

public:
    Fields() {}
    Fields( const json& J );

    void fromJSON( const json& J );
private:
  std::map<std::string,std::shared_ptr<FieldHolder>> fields_;

};


} // ursa

#endif
