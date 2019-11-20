#include "simplex.h"

namespace luna
{

template<typename T>
void
Simplex::transfer( const Simplex& master , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  luna_assert( Y.size()==nb_basis() );
  luna_assert( X.size()==master.nb_basis() );

  std::vector<real_t> phi( master.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // get the reference coordinate for Y
    const real_t* y = reference_.get_reference_coordinate(k);

    // evaluate the basis function of the original master at this coordinate
    master.basis().evaluate( y , phi.data() );

    // evaluate the basis functions in the original master element
    for (coord_t d=0;d<dim;d++)
    {
      Y[k][d] = 0;
      for (index_t j=0;j<master.nb_basis();j++)
        Y[k][d] += phi[j]*X[j][d];
    }
  }
}

template<typename T>
void
Simplex::transfer( const Simplex& master , const std::vector<const T>& X , std::vector<T>& Y ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  luna_assert( Y.size()==nb_basis() );
  luna_assert( X.size()==master.nb_basis() );

  std::vector<real_t> phi( master.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // get the reference coordinate for Y
    const real_t* y = reference_.get_reference_coordinate(k);

    // evaluate the basis function of the original master at this coordinate
    master.basis().evaluate( y , phi.data() );

    // evaluate the basis functions in the original master element
    for (index_t j=0;j<master.nb_basis();j++)
      Y[k] += phi[j]*X[j];
  }
}

#if 0
template<typename BasisFrom_t,typename T>
void
Simplex<Bezier>::transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  luna_assert( Y.size()==nb_basis() );
  luna_assert( X.size()==master.nb_basis() );

  std::vector<real_t> phi( master.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // evaluate all basis functions associated with this point
    luna_implement;

    // evaluate the basis functions in the original master element
    for (coord_t d=0;d<dim;d++)
    {
      Y[k][d] = 0;
      for (index_t j=0;j<master.nb_basis();j++)
        Y[k][d] += phi[j]*X[j][d];
    }
  }
}

template<typename BasisFrom_t,typename T>
void
Simplex<Bezier>::transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T>& X , std::vector<T>& Y ) const
{
  // assume the coordinates in Y have correspond to the local coordinates
  luna_assert( Y.size()==nb_basis() );
  luna_assert( X.size()==master.nb_basis() );

  std::vector<real_t> phi( master.nb_basis() , 0. );
  for (index_t k=0;k<Y.size();k++)
  {
    // evaluate all basis functions associated with this point
    luna_implement;

    // evaluate the basis functions in the original master element
    for (index_t j=0;j<master.nb_basis();j++)
      Y[k] += phi[j]*X[j];
  }
}

#endif

} // luna
