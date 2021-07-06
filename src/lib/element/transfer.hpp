//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "simplex.h"

namespace avro
{

template<typename T>
void
Simplex::transfer( const Simplex& element , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  avro_assert( Y.size()==nb_basis() );
  avro_assert( X.size()==element.nb_basis() );

  std::vector<real_t> phi( element.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // get the reference coordinate for Y
    const real_t* y = reference_.get_reference_coordinate(k);

    // evaluate the basis function of the original  at this coordinate
    element.reference().basis().evaluate( y , phi.data() );

    // evaluate the basis functions in the original  element
    for (coord_t d=0;d<dim;d++)
    {
      Y[k][d] = 0;
      for (index_t j=0;j<element.nb_basis();j++)
        Y[k][d] += phi[j]*X[j][d];
    }
  }
}

template<typename T>
void
Simplex::transfer( const Simplex& element , const std::vector<const T>& X , std::vector<T>& Y ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  avro_assert( Y.size()==nb_basis() );
  avro_assert( X.size()==element.nb_basis() );

  std::vector<real_t> phi( element.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // get the reference coordinate for Y
    const real_t* y = reference_.get_reference_coordinate(k);

    // evaluate the basis function of the original  at this coordinate
    element.reference().basis().evaluate( y , phi.data() );

    // evaluate the basis functions in the original  element
    for (index_t j=0;j<element.nb_basis();j++)
      Y[k] += phi[j]*X[j];
  }
}

#if 0
template<typename BasisFrom_t,typename T>
void
Simplex<Bezier>::transfer( const Simplex<BasisFrom_t>& element , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  avro_assert( Y.size()==nb_basis() );
  avro_assert( X.size()==element.nb_basis() );

  std::vector<real_t> phi( element.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // evaluate all basis functions associated with this point
    avro_implement;

    // evaluate the basis functions in the original  element
    for (coord_t d=0;d<dim;d++)
    {
      Y[k][d] = 0;
      for (index_t j=0;j<element.nb_basis();j++)
        Y[k][d] += phi[j]*X[j][d];
    }
  }
}

template<typename BasisFrom_t,typename T>
void
Simplex<Bezier>::transfer( const Simplex<BasisFrom_t>&  , const std::vector<const T>& X , std::vector<T>& Y ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  avro_assert( Y.size()==nb_basis() );
  avro_assert( X.size()==.nb_basis() );

  std::vector<real_t> phi( .nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // evaluate all basis functions associated with this point
    avro_implement;

    // evaluate the basis functions in the original  element
    for (index_t j=0;j<.nb_basis();j++)
      Y[k] += phi[j]*X[j];
  }
}

#endif

} // avro
