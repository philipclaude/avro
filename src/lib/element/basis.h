//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MASTER_BASIS_H_
#define avro_LIB_MASTER_BASIS_H_

#include "avro_types.h"

#include "element/reference.h"

namespace avro
{

template<typename type> class Bezier;
template<typename type> class Lagrange;
class Simplex;

enum BasisFunctionCategory
{
  BasisFunctionCategory_Lagrange,
  BasisFunctionCategory_Legendre,
  BasisFunctionCategory_Bezier,
  BasisFunctionCategory_None
};

template<>
class Lagrange<Simplex>
{
public:
  static void eval( const ReferenceElement<Simplex>& ref , const double* x , double* phi );
  static void grad( const ReferenceElement<Simplex>& ref , const double* x , double* gphi );
  static void hess( const ReferenceElement<Simplex>& ref , const double* x , double* hphi );
};

template<>
class Bezier<Simplex>
{
public:
  static void eval(const ReferenceElement<Simplex>& ref , const double* x , double* phi );
  static void grad( const ReferenceElement<Simplex>& ref , const double* x , double* gphi );
  static void hess(const ReferenceElement<Simplex>& ref , const double* x , double* hphi );
};

template<typename type>
class Basis
{
private:

  typedef void (*eval_func_ptr)(const ReferenceElement<type>&,const double*,double*);
  typedef void (*eval_grad_ptr)(const ReferenceElement<type>&,const double*,double*);
  typedef void (*eval_hess_ptr)(const ReferenceElement<type>&,const double*,double*);

  eval_func_ptr
  get_func( BasisFunctionCategory category )
  {
    if (category==BasisFunctionCategory_Lagrange)
      return Lagrange<type>::eval;
    if (category==BasisFunctionCategory_Bezier)
      return Bezier<type>::eval;
    return NULL;
  }

  eval_grad_ptr
  get_grad( BasisFunctionCategory category )
  {
    if (category==BasisFunctionCategory_Lagrange)
      return Lagrange<type>::grad;
    if (category==BasisFunctionCategory_Bezier)
      return Bezier<type>::grad;
    return NULL;
  }

  eval_hess_ptr
  get_hess( BasisFunctionCategory category )
  {
    if (category==BasisFunctionCategory_Lagrange)
      return Lagrange<type>::hess;
    if (category==BasisFunctionCategory_Bezier)
      return Bezier<type>::hess;
    return NULL;
  }

public:

  Basis( const ReferenceElement<type>& reference , const BasisFunctionCategory category ) :
    reference_(reference),
    category_(category),
    fptr_( get_func(category_) ),
    gptr_( get_grad(category_) ),
    hptr_( get_hess(category_) )
  {}

  void evaluate( const double* x , double* phi) const
  {
    fptr_(reference_,x,phi);
  }

  void gradient( const double* x , double* gphi ) const
  {
    gptr_(reference_,x,gphi);
  }

  void hessian( const double* x , double* hphi ) const
  {
    hptr_(reference_,x,hphi);
  }

private:
  const ReferenceElement<type>& reference_;
  BasisFunctionCategory category_;

  const eval_func_ptr fptr_;
  const eval_grad_ptr gptr_;
  const eval_hess_ptr hptr_;
};

} // avro

#endif
