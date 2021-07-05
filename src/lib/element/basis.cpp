#include "common/error.h"
#include "element/basis.h"

namespace avro
{

template<typename type>
typename Basis<type>::eval_func_ptr
Basis<type>::get_func( coord_t number, coord_t order, BasisFunctionCategory category )
{
  if (category == BasisFunctionCategory_Lagrange) {
    if (number == 1) {
      static const int N = 1;
      if (order == 1)      return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 1)      return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 1)      return Lagrange<type,N,1>::eval;
      else if (order == 2) return Lagrange<type,N,2>::eval;
      else if (order == 3) return Lagrange<type,N,3>::eval;
      else if (order == 4) return Lagrange<type,N,4>::eval;
      else if (order == 5) return Lagrange<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 1)      return Lagrange<type,N,1>::eval;
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
  if (category == BasisFunctionCategory_Bezier) {
    if (number == 1) {
      static const int N = 1;
      if (order == 1)      return Bezier<type,N,1>::eval;
      else if (order == 2) return Bezier<type,N,2>::eval;
      else if (order == 3) return Bezier<type,N,3>::eval;
      else if (order == 4) return Bezier<type,N,4>::eval;
      else if (order == 5) return Bezier<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 1)      return Bezier<type,N,1>::eval;
      else if (order == 2) return Bezier<type,N,2>::eval;
      else if (order == 3) return Bezier<type,N,3>::eval;
      else if (order == 4) return Bezier<type,N,4>::eval;
      else if (order == 5) return Bezier<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 1)      return Bezier<type,N,1>::eval;
      else if (order == 2) return Bezier<type,N,2>::eval;
      else if (order == 3) return Bezier<type,N,3>::eval;
      else if (order == 4) return Bezier<type,N,4>::eval;
      else if (order == 5) return Bezier<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 1)      return Bezier<type,N,1>::eval;
      else if (order == 2) return Bezier<type,N,2>::eval;
      else if (order == 3) return Bezier<type,N,3>::eval;
      else if (order == 4) return Bezier<type,N,4>::eval;
      else if (order == 5) return Bezier<type,N,5>::eval;
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
      if (order == 1)      return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 2) {
      static const int N = 2;
      if (order == 1)      return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 3) {
      static const int N = 3;
      if (order == 1)      return Legendre<type,N,1>::eval;
      else if (order == 2) return Legendre<type,N,2>::eval;
      else if (order == 3) return Legendre<type,N,3>::eval;
      else if (order == 4) return Legendre<type,N,4>::eval;
      else if (order == 5) return Legendre<type,N,5>::eval;
      else avro_implement;
    }
    else if (number == 4) {
      static const int N = 4;
      if (order == 1)      return Legendre<type,N,1>::eval;
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

template class Basis<Simplex>;

} // avro
