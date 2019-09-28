#ifndef URSA_LIB_NUMERICS_FIELD_H_
#define URSA_LIB_NUMERICS_FIELD_H_

#include "common/types.h"

#include "master/master.h"

namespace ursa
{

template<typename T> class Data;
template<typename type> class Topology;

class Parameter;
class Coordinate;

class FieldBase
{};

// a field of T's defined on a mesh with elements that have a master element Master_t
template<typename T>
class _Field : public FieldBase
{
public:
  T& eval( index_t elem , const Parameter& u ) const;
  T& eval( index_t elem , const Coordinate& x ) const;
  T& eval( const Coordinate& x ) const;

protected:
  Data<T> data_;
};

template<typename Master_t,typename T> class Field;

template<typename Basis,typename T>
class Field<Simplex<Basis>,T> : public _Field<T>
{
  typedef Simplex<Basis>    Master_t;
  typedef Simplex<Lagrange> Shape_t;

public:
  Field( const Topology<Shape_t>& , coord_t order );

private:
  const Master_t master_;
  const Topology<Shape_t>& topology_;
};

template<typename T>
class Field<Polytope,T> : public _Field<T>
{
  typedef Polytope Shape_t;
  typedef Polytope Master_t;

public:
  Field( Topology<Shape_t>& topology , coord_t order );

private:
  const Topology<Shape_t>& topology_;
  const Master_t master_;
};


} // ursa

#endif
