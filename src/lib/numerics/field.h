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

// a field of T's defined on a mesh with elements that have a master element M
template<typename M,typename T>
class Field
{

public:
  Field( const Topology<M>& topology , int order );


  T& eval( index_t elem , const Parameter& u ) const;
  T& eval( index_t elem , const Coordinate& x ) const;
  T& eval( const Coordinate& x ) const;


protected:
  M master_;
  Data<T> data_;
  const Topology<M>& topology_;

};

} // ursa

#endif
