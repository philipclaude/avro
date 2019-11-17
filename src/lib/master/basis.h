#ifndef LUNA_LIB_MASTER_BASIS_H_
#define LUNA_LIB_MASTER_BASIS_H_

namespace luna
{

#if 1

class Lagrange;
class Bezier;
class Legendre;

#else

template<typename Shape> class Lagrange;
template<typename Shape> class Bezier;
template<typename Shape> class Legendre;

enum BasisFunctionCategory
{
  BasisFunctionCategory_Lagrange,
  BasisFunctionCategory_Legendre,
  BasisFunctionCategory_Bezier,
  BasisFunctionCategory_None
};

typedef void (*eval_func_ptr)(coord_t,coord_t,double*,double*);
typedef void (*eval_grad_ptr)(coord_t,coord_t,double*,double*);
typedef void (*eval_hess_ptr)(coord_t,coord_t,double*,double*);

template<typename Shape>
const eval_func_ptr
get_func( BasisFunctionCategory category )
{
  if (category==BasisFunctionCategory_Lagrange)
    return Lagrange<Shape>::eval;
  if (category==BasisFunctionCategory_Bezier)
    return Bezier<Shape>::eval;
  return NULL;
}

template<typename Shape>
const eval_grad_ptr
get_grad( BasisFunctionCategory category )
{
  if (category==BasisFunctionCategory_Lagrange)
    return Lagrange<Shape>::grad;
  if (category==BasisFunctionCategory_Bezier)
    return Bezier<Shape>::grad;
  return NULL;
}

template<typename Shape>
const eval_hess_ptr
get_hess( BasisFunctionCategory category )
{
  if (category==BasisFunctionCategory_Lagrange)
    return Lagrange<Shape>::hess;
  if (category==BasisFunctionCategory_Bezier)
    return Bezier<Shape>::hess;
  return NULL;
}

template<typename Shape>
class BasisFunction
{
public:
  BasisFunction( coord_t number , coord_t order , const BasisFunctionCategory category ) :
    number_(number),
    order_(order),
    category_(category),
    fptr_( get_func<Shape>(category_) ),
    gptr_( get_grad<Shape>(category_) ),
    hptr_( get_hess<Shape>(category_) )
  {}

  void eval( double* x , double* phi) const
  {
    fptr_(number_,order_,x,phi);
  }

  void grad( double* x , double* gphi ) const
  {
    gptr_(number_,order_,x,gphi);
  }

  void hess( double* x , double* hphi ) const
  {
    hptr_(number_,order_,x,hphi);
  }

private:
  const coord_t number_;
  const coord_t order_;
  BasisFunctionCategory category_;

  const eval_func_ptr fptr_;
  const eval_grad_ptr gptr_;
  const eval_hess_ptr hptr_;
};

#endif

} // luna

#endif
