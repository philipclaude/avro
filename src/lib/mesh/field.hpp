#include "mesh/field.h"

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

} // avro
