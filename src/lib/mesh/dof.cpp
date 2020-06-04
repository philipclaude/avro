#include "adaptation/metric.h"

#include "mesh/dof.h"

#include "numerics/linear_algebra.h"

namespace avro
{

template<>
bool
DOF<Metric>::interpolate( const index_t* idx , index_t nv , const std::vector<real_t>& phi , Metric* pT ) const
{
  avro_assert_msg( nv == phi.size() , "nv = %lu, |phi| = %lu" , nv , phi.size() );

  numerics::SymMatrixD<real_t> T(*pT);

  T = 0;
  for (index_t k=0;k<nv;k++)
	{
		T = T + numerics::logm( (*this)( idx[k] , 0 ) )*phi[k];
	}
	T = numerics::expm(T);

  real_t d = numerics::determinant(T);
  if (d<=0 || std::isnan(d))
  {
    for (index_t k=0;k<nv;k++)
    {
      std::cout << (*this)( idx[k],0 ) << std::endl;
      std::cout << numerics::logm( (*this)(idx[k],0)) << std::endl;
    }
    print_inline(phi);
    std::cout << T << std::endl;
    return false;
  }
  pT->set(T);
  return true;
}

} // avro
