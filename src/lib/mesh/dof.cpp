//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/metric.h"

#include "mesh/dof.h"

#include "numerics/linear_algebra.h"

namespace avro
{

template<typename type>
bool
DOF<type>::interpolate( const index_t* idx , index_t nv , const std::vector<real_t>& phi , type* u ) const
{
  avro_assert_msg( nv == phi.size() , "nv = %lu, |phi| = %lu" , nv , phi.size() );
  for (index_t j=0;j<rank();j++)
    u[j] = type(0);
  for (index_t j=0;j<rank();j++)
  {
    for (index_t i=0;i<phi.size();i++)
      u[j] = u[j] + (*this)( idx[i] , 0 )*phi[i]; // should it be j instead of 0?
  }
  return true;
}

template<>
bool
DOF<Metric>::interpolate( const index_t* idx , index_t nv , const std::vector<real_t>& phi , Metric* pT ) const
{
  avro_assert_msg( nv == phi.size() , "nv = %lu, |phi| = %lu" , nv , phi.size() );

  symd<real_t> T(*pT);

  T.zero();
  for (index_t k=0;k<nv;k++)
	{
		T = T + numerics::logm( (*this)( idx[k] , 0 ) )*phi[k];
	}
	T = numerics::expm(T);

  #if 0
  // check the eigenvalues
  std::pair< vecd<real_t> , matd<real_t> > decomp = numerics::eig(T);
  real_t hmin = 1e-8;
  for (index_t d = 0; d < T.n(); d++) {
    real_t h = 1./sqrt(decomp.first(d));
    avro_assert( h > hmin );
  }
  #endif

  real_t d = numerics::det(T);
  if (d<=0 || std::isnan(d))
  {
    // find all nonzero coefficients
    #if 0
    std::vector<real_t> alpha;
    std::vector<index_t> idx_alpha;
    for (index_t k = 0; k < nv; k++) {
      if (phi[k] < 1e-8) continue;
      alpha.push_back( phi[k] );
      idx_alpha.push_back( idx[k] );
    }

    bool result = this->interpolate( idx_alpha.data() , idx_alpha.size() , alpha, pT );
    avro_assert( result );
    return true;
    #endif

    for (index_t k=0;k<nv;k++)
    {
      //std::cout << (*this)( idx[k],0 ) << std::endl;
      //std::cout << numerics::logm( (*this)(idx[k],0)) << std::endl;
    }
    //print_inline(phi);
    //std::cout << T << std::endl;
    return false;
  }
  pT->set(T);
  return true;
}

template class DOF<double>;

} // avro
