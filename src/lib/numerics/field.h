#ifndef URSA_LIB_NUMERICS_FIELD_H_
#define URSA_LIB_NUMERICS_FIELD_H_

namespace ursa
{

template<typename T> class Data;
template<typename type> class Topology;

// a field of T's defined on a mesh with elements that have a master element M
template<typename M,typename T>
class Field
{

public:
  Field( Topology<M>& topology );


  T& eval();
  const T& eval() const;

protected:
  Data<T> data_;
  Topology<M>& topology_;

};

} // ursa

#endif
