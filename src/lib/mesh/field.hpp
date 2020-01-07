#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/field.h"
#include "mesh/points.h"
#include "mesh/topology.h"

namespace avro
{

template<typename T>
template<typename Function>
void
Field<Simplex,T>::evaluate( const Function& function )
{
  std::vector<real_t> x( topology_.points().dim() );
  std::vector<const real_t*> xk;
  std::vector<real_t> phi( topology_.master().nb_basis() );

  for (index_t k=0;k<topology_.nb();k++)
  {
    xk.resize( topology_.nv(k) );
    for (index_t j=0;j<topology_.nv(k);j++)
      xk[j] = topology_.points()[topology_(k,j)];

    // loop through the reference points of the master simplex
    for (index_t j=0;j<master_.nb_basis();j++)
    {
      const real_t* xref = master_.reference().get_reference_coordinate(j);

      // evaluate the basis functions at the quadrature point
      topology_.master().basis().evaluate( xref , phi.data() );

      // evaluate the physical coordinates
      topology_.points().interpolate( xk , phi , x.data() );

      (*this)(k,j) = function(x.data());
    }
  }
}

template<typename T>
void
Field<Polytope,T>::evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const
{
  avro_assert( master_.order() == 0 );
  avro_assert( this->type() == DISCONTINUOUS );

  result.resize( parents.size() );
  for (index_t k=0;k<parents.size();k++)
  {
    // compute a linear combination of the vertex dof
    result[k] = this->data_[parents[k]][0];
  }
}

template<typename T>
real_t
FieldBase<T>::min( index_t rank ) const
{
  avro_implement;
}

template<typename T>
real_t
FieldBase<T>::max( index_t rank ) const
{
  avro_implement;
}

template<>
inline real_t
FieldBase<real_t>::min( index_t rank ) const
{
  const std::vector<real_t>& u = data_.template DOF<real_t>::data();
  return *std::min_element( u.begin() , u.end() );
}

template<>
inline real_t
FieldBase<real_t>::max( index_t rank ) const
{
  const std::vector<real_t>& u = data_.template DOF<real_t>::data();
  return *std::max_element( u.begin() , u.end() );
}

} // avro
