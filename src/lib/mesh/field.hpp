//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/field.h"
#include "mesh/points.h"
#include "mesh/topology.h"

namespace avro
{

template<typename T> real_t __at_rank__( const T& x , index_t r );

template<typename T>
struct __comparator__
{
  __comparator__( index_t rank ) :
    rank(rank)
  {}

  bool operator() ( const T& x , const T& y ) const
  {
    return __at_rank__(x,rank) < __at_rank__(y,rank);
  }

  index_t rank;
};

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
Field<Simplex,T>::evaluate( index_t rank , const std::vector<index_t>& parents , const Table<real_t>& alpha , std::vector<real_t>& result ) const
{
  avro_assert( parents.size() == alpha.nb() );
  std::vector<real_t> phi( master_.nb_basis() );

  result.resize( parents.size() );
  for (index_t k=0;k<parents.size();k++)
  {
    // evaluate the basis functions for this reference coordinate
    master_.basis().evaluate( alpha(k) , phi.data() );

    // evaluate the field
    result[k] = 0.0;
    for (index_t j=0;j<phi.size();j++)
      result[k] += __at_rank__( (*this)(parents[k],j) , rank )*phi[j];
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
  const std::vector<T>& u = data_.DOF<T>::data();
  typename std::vector<T>::const_iterator it = std::min_element( u.begin() , u.end() , __comparator__<T>(rank) );
  return __at_rank__<T>(*it,rank);
}

template<typename T>
real_t
FieldBase<T>::max( index_t rank ) const
{
  const std::vector<T>& u = data_.DOF<T>::data();
  typename std::vector<T>::const_iterator it = std::max_element( u.begin() , u.end() , __comparator__<T>(rank) );
  return __at_rank__<T>(*it,rank);
}

} // avro
