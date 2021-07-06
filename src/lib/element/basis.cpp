#include "common/error.h"
#include "element/basis.h"

namespace avro
{

// yes, a p = 0 Lagrange element doesn't make sense, but it's okay
template<int N>
void constant_basis_func( const real_t* x , real_t* phi ) {
  phi[0] = 1;
}

template<int N>
void constant_basis_grad( const real_t* x , real_t* gphi ) {
  for (int i = 0; i < N; i++)
    gphi[i] = 0.0;
}

template<int N>
void constant_basis_hess( const real_t* x , real_t* hphi ) {
  for (int i = 0; i < N*N; i++)
    hphi[i] = 0.0;
}

template<typename type>
typename Basis<type>::eval_func_ptr
Basis<type>::get_func( coord_t number, coord_t order, BasisFunctionCategory category )
{
  if (category == BasisFunctionCategory_Lagrange) {
    if (number == 1) {
      static const int N = 1;
      if (order == 0)      return &constant_basis_func<N>;
      else if (order == 1) return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 0)      return &constant_basis_func<N>;
      else if (order == 1) return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 0)      return &constant_basis_func<N>;
      else if (order == 1) return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 0)      return &constant_basis_func<N>;
      else if (order == 1) return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  if (category == BasisFunctionCategory_Bernstein) {
    if (number == 1) {
      static const int N = 1;
      if (order == 1)      return Bernstein<type,N,1>::eval;
      else if (order == 2) return Bernstein<type,N,2>::eval;
      else if (order == 3) return Bernstein<type,N,3>::eval;
      else if (order == 4) return Bernstein<type,N,4>::eval;
      else if (order == 5) return Bernstein<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 1)      return Bernstein<type,N,1>::eval;
      else if (order == 2) return Bernstein<type,N,2>::eval;
      else if (order == 3) return Bernstein<type,N,3>::eval;
      else if (order == 4) return Bernstein<type,N,4>::eval;
      else if (order == 5) return Bernstein<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 1)      return Bernstein<type,N,1>::eval;
      else if (order == 2) return Bernstein<type,N,2>::eval;
      else if (order == 3) return Bernstein<type,N,3>::eval;
      else if (order == 4) return Bernstein<type,N,4>::eval;
      else if (order == 5) return Bernstein<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 1)      return Bernstein<type,N,1>::eval;
      else if (order == 2) return Bernstein<type,N,2>::eval;
      else if (order == 3) return Bernstein<type,N,3>::eval;
      else if (order == 4) return Bernstein<type,N,4>::eval;
      else if (order == 5) return Bernstein<type,N,5>::eval;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  if (category == BasisFunctionCategory_Legendre) {
    if (number == 1) {
      static const int N = 1;
      if (order == 0)      return Legendre<type,N,0>::eval;
      else if (order == 1) return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 0)      return Legendre<type,N,0>::eval;
      else if (order == 1) return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 0)      return Legendre<type,N,0>::eval;
      else if (order == 1) return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 0)      return Legendre<type,N,0>::eval;
      else if (order == 1) return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  return nullptr;
}

template<typename type>
typename Basis<type>::eval_grad_ptr
Basis<type>::get_grad( coord_t number, coord_t order, BasisFunctionCategory category )
{
  if (category == BasisFunctionCategory_Lagrange) {
    if (number == 1) {
      static const int N = 1;
      if (order == 0)      return &constant_basis_grad<N>;
      else if (order == 1) return Lagrange<type,N,1>::grad;
      else if (order == 2) return Lagrange<type,N,2>::grad;
      else if (order == 3) return Lagrange<type,N,3>::grad;
      else if (order == 4) return Lagrange<type,N,4>::grad;
      else if (order == 5) return Lagrange<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 0)      return &constant_basis_grad<N>;
      else if (order == 1) return Lagrange<type,N,1>::grad;
      else if (order == 2) return Lagrange<type,N,2>::grad;
      else if (order == 3) return Lagrange<type,N,3>::grad;
      else if (order == 4) return Lagrange<type,N,4>::grad;
      else if (order == 5) return Lagrange<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 0)      return &constant_basis_grad<N>;
      else if (order == 1) return Lagrange<type,N,1>::grad;
      else if (order == 2) return Lagrange<type,N,2>::grad;
      else if (order == 3) return Lagrange<type,N,3>::grad;
      else if (order == 4) return Lagrange<type,N,4>::grad;
      else if (order == 5) return Lagrange<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 0)      return &constant_basis_grad<N>;
      else if (order == 1) return Lagrange<type,N,1>::grad;
      else if (order == 2) return Lagrange<type,N,2>::grad;
      else if (order == 3) return Lagrange<type,N,3>::grad;
      else if (order == 4) return Lagrange<type,N,4>::grad;
      else if (order == 5) return Lagrange<type,N,5>::grad;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  if (category == BasisFunctionCategory_Bernstein) {
    if (number == 1) {
      static const int N = 1;
      if (order == 1)      return Bernstein<type,N,1>::grad;
      else if (order == 2) return Bernstein<type,N,2>::grad;
      else if (order == 3) return Bernstein<type,N,3>::grad;
      else if (order == 4) return Bernstein<type,N,4>::grad;
      else if (order == 5) return Bernstein<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 1)      return Bernstein<type,N,1>::grad;
      else if (order == 2) return Bernstein<type,N,2>::grad;
      else if (order == 3) return Bernstein<type,N,3>::grad;
      else if (order == 4) return Bernstein<type,N,4>::grad;
      else if (order == 5) return Bernstein<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 1)      return Bernstein<type,N,1>::grad;
      else if (order == 2) return Bernstein<type,N,2>::grad;
      else if (order == 3) return Bernstein<type,N,3>::grad;
      else if (order == 4) return Bernstein<type,N,4>::grad;
      else if (order == 5) return Bernstein<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 1)      return Bernstein<type,N,1>::grad;
      else if (order == 2) return Bernstein<type,N,2>::grad;
      else if (order == 3) return Bernstein<type,N,3>::grad;
      else if (order == 4) return Bernstein<type,N,4>::grad;
      else if (order == 5) return Bernstein<type,N,5>::grad;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  if (category == BasisFunctionCategory_Legendre) {
    if (number == 1) {
      static const int N = 1;
      if (order == 0)      return Legendre<type,N,0>::grad;
      else if (order == 1) return Legendre<type,N,1>::grad;
      else if (order == 2) return Legendre<type,N,2>::grad;
      else if (order == 3) return Legendre<type,N,3>::grad;
      else if (order == 4) return Legendre<type,N,4>::grad;
      else if (order == 5) return Legendre<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 0)      return Legendre<type,N,0>::grad;
      else if (order == 1) return Legendre<type,N,1>::grad;
      else if (order == 2) return Legendre<type,N,2>::grad;
      else if (order == 3) return Legendre<type,N,3>::grad;
      else if (order == 4) return Legendre<type,N,4>::grad;
      else if (order == 5) return Legendre<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 0)      return Legendre<type,N,0>::grad;
      else if (order == 1) return Legendre<type,N,1>::grad;
      else if (order == 2) return Legendre<type,N,2>::grad;
      else if (order == 3) return Legendre<type,N,3>::grad;
      else if (order == 4) return Legendre<type,N,4>::grad;
      else if (order == 5) return Legendre<type,N,5>::grad;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 0)      return Legendre<type,N,0>::grad;
      else if (order == 1) return Legendre<type,N,1>::grad;
      else if (order == 2) return Legendre<type,N,2>::grad;
      else if (order == 3) return Legendre<type,N,3>::grad;
      else if (order == 4) return Legendre<type,N,4>::grad;
      else if (order == 5) return Legendre<type,N,5>::grad;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  return nullptr;
}

template<typename type>
typename Basis<type>::eval_hess_ptr
Basis<type>::get_hess( coord_t number, coord_t order, BasisFunctionCategory category )
{
  if (category == BasisFunctionCategory_Lagrange) {
    if (number == 1) {
      static const int N = 1;
      if (order == 0)      return &constant_basis_hess<N>;
      else if (order == 1) return Lagrange<type,N,1>::hess;
      else if (order == 2) return Lagrange<type,N,2>::hess;
      else if (order == 3) return Lagrange<type,N,3>::hess;
      else if (order == 4) return Lagrange<type,N,4>::hess;
      else if (order == 5) return Lagrange<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 0)      return &constant_basis_hess<N>;
      else if (order == 1) return Lagrange<type,N,1>::hess;
      else if (order == 2) return Lagrange<type,N,2>::hess;
      else if (order == 3) return Lagrange<type,N,3>::hess;
      else if (order == 4) return Lagrange<type,N,4>::hess;
      else if (order == 5) return Lagrange<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 0)      return &constant_basis_hess<N>;
      else if (order == 1) return Lagrange<type,N,1>::hess;
      else if (order == 2) return Lagrange<type,N,2>::hess;
      else if (order == 3) return Lagrange<type,N,3>::hess;
      else if (order == 4) return Lagrange<type,N,4>::hess;
      else if (order == 5) return Lagrange<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 0)      return &constant_basis_hess<N>;
      else if (order == 1) return Lagrange<type,N,1>::hess;
      else if (order == 2) return Lagrange<type,N,2>::hess;
      else if (order == 3) return Lagrange<type,N,3>::hess;
      else if (order == 4) return Lagrange<type,N,4>::hess;
      else if (order == 5) return Lagrange<type,N,5>::hess;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  if (category == BasisFunctionCategory_Bernstein) {
    if (number == 1) {
      static const int N = 1;
      if (order == 1)      return Bernstein<type,N,1>::hess;
      else if (order == 2) return Bernstein<type,N,2>::hess;
      else if (order == 3) return Bernstein<type,N,3>::hess;
      else if (order == 4) return Bernstein<type,N,4>::hess;
      else if (order == 5) return Bernstein<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 1)      return Bernstein<type,N,1>::hess;
      else if (order == 2) return Bernstein<type,N,2>::hess;
      else if (order == 3) return Bernstein<type,N,3>::hess;
      else if (order == 4) return Bernstein<type,N,4>::hess;
      else if (order == 5) return Bernstein<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 1)      return Bernstein<type,N,1>::hess;
      else if (order == 2) return Bernstein<type,N,2>::hess;
      else if (order == 3) return Bernstein<type,N,3>::hess;
      else if (order == 4) return Bernstein<type,N,4>::hess;
      else if (order == 5) return Bernstein<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 1)      return Bernstein<type,N,1>::hess;
      else if (order == 2) return Bernstein<type,N,2>::hess;
      else if (order == 3) return Bernstein<type,N,3>::hess;
      else if (order == 4) return Bernstein<type,N,4>::hess;
      else if (order == 5) return Bernstein<type,N,5>::hess;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  if (category == BasisFunctionCategory_Legendre) {
    if (number == 1) {
      static const int N = 1;
      if (order == 0)      return Legendre<type,N,0>::hess;
      else if (order == 1) return Legendre<type,N,1>::hess;
      else if (order == 2) return Legendre<type,N,2>::hess;
      else if (order == 3) return Legendre<type,N,3>::hess;
      else if (order == 4) return Legendre<type,N,4>::hess;
      else if (order == 5) return Legendre<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 0)      return Legendre<type,N,0>::hess;
      else if (order == 1) return Legendre<type,N,1>::hess;
      else if (order == 2) return Legendre<type,N,2>::hess;
      else if (order == 3) return Legendre<type,N,3>::hess;
      else if (order == 4) return Legendre<type,N,4>::hess;
      else if (order == 5) return Legendre<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 0)      return Legendre<type,N,0>::hess;
      else if (order == 1) return Legendre<type,N,1>::hess;
      else if (order == 2) return Legendre<type,N,2>::hess;
      else if (order == 3) return Legendre<type,N,3>::hess;
      else if (order == 4) return Legendre<type,N,4>::hess;
      else if (order == 5) return Legendre<type,N,5>::hess;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 0)      return Legendre<type,N,0>::hess;
      else if (order == 1) return Legendre<type,N,1>::hess;
      else if (order == 2) return Legendre<type,N,2>::hess;
      else if (order == 3) return Legendre<type,N,3>::hess;
      else if (order == 4) return Legendre<type,N,4>::hess;
      else if (order == 5) return Legendre<type,N,5>::hess;
      else avro_implement;
    }
    else {
      printf("unsupported shape dimension %u\n",number);
      avro_implement;
    }
  }
  return nullptr;
}

template class Basis<Simplex>;

} // avro
