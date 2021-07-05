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

template<typename type,int number,int order> class Bezier;
template<typename type,int number,int order> class Lagrange;
template<typename type,int number,int order> class Legendre;

class Simplex;

enum BasisFunctionCategory
{
  BasisFunctionCategory_Lagrange,
  BasisFunctionCategory_Legendre,
  BasisFunctionCategory_Bezier,
  BasisFunctionCategory_None
};


template<int N,int P>
class Lagrange<Simplex,N,P>
{
public:
  static void eval( const double* x , double* phi );
  static void grad( const double* x , double* gphi );
  static void hess( const double* x , double* hphi );
};

template<int N,int P>
class Bezier<Simplex,N,P>
{
public:
  static void eval( const double* x , double* phi );
  static void grad( const double* x , double* gphi );
  static void hess( const double* x , double* hphi );
};

template<int N,int P>
class Legendre<Simplex,N,P>
{
public:
  static void eval( const double* x , double* phi );
  static void grad( const double* x , double* gphi );
  static void hess( const double* x , double* hphi );
};

template<typename type>
class Basis
{
private:

  typedef void (*eval_func_ptr)(const double*,double*);
  typedef void (*eval_grad_ptr)(const double*,double*);
  typedef void (*eval_hess_ptr)(const double*,double*);

  /*
  eval_func_ptr
  get_func( coord_t number, coord_t order, BasisFunctionCategory category )
  {
    if (category == BasisFunctionCategory_Lagrange) {
      if (number == 1) {
        static int N = 1;
        if (order == 1)      return Lagrange<type,N,1>::eval;
        else if (order == 2) return Lagrange<type,N,2>::eval;
        else if (order == 3) return Lagrange<type,N,3>::eval;
        else if (order == 4) return Lagrange<type,N,4>::eval;
        else if (order == 5) return Lagrange<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 2) {
        static int N = 2;
        if (order == 1)      return Lagrange<type,N,1>::eval;
        else if (order == 2) return Lagrange<type,N,2>::eval;
        else if (order == 3) return Lagrange<type,N,3>::eval;
        else if (order == 4) return Lagrange<type,N,4>::eval;
        else if (order == 5) return Lagrange<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 3) {
        static int N = 3;
        if (order == 1)      return Lagrange<type,N,1>::eval;
        else if (order == 2) return Lagrange<type,N,2>::eval;
        else if (order == 3) return Lagrange<type,N,3>::eval;
        else if (order == 4) return Lagrange<type,N,4>::eval;
        else if (order == 5) return Lagrange<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 4) {
        static int N = 4;
        if (order == 1)      return Lagrange<type,N,1>::eval;
        else if (order == 2) return Lagrange<type,N,2>::eval;
        else if (order == 3) return Lagrange<type,N,3>::eval;
        else if (order == 4) return Lagrange<type,N,4>::eval;
        else if (order == 5) return Lagrange<type,N,5>::eval;
        else avro_implement;
      }
      else {
        printf("unsupported shape dimension %lu\n",number);
        avro_implement;
      }
    }
    if (category == BasisFunctionCategory_Bezier) {
      if (number == 1) {
        static int N = 1;
        if (order == 1)      return Bezier<type,N,1>::eval;
        else if (order == 2) return Bezier<type,N,2>::eval;
        else if (order == 3) return Bezier<type,N,3>::eval;
        else if (order == 4) return Bezier<type,N,4>::eval;
        else if (order == 5) return Bezier<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 2) {
        static int N = 2;
        if (order == 1)      return Bezier<type,N,1>::eval;
        else if (order == 2) return Bezier<type,N,2>::eval;
        else if (order == 3) return Bezier<type,N,3>::eval;
        else if (order == 4) return Bezier<type,N,4>::eval;
        else if (order == 5) return Bezier<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 3) {
        static int N = 3;
        if (order == 1)      return Bezier<type,N,1>::eval;
        else if (order == 2) return Bezier<type,N,2>::eval;
        else if (order == 3) return Bezier<type,N,3>::eval;
        else if (order == 4) return Bezier<type,N,4>::eval;
        else if (order == 5) return Bezier<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 4) {
        static int N = 4;
        if (order == 1)      return Bezier<type,N,1>::eval;
        else if (order == 2) return Bezier<type,N,2>::eval;
        else if (order == 3) return Bezier<type,N,3>::eval;
        else if (order == 4) return Bezier<type,N,4>::eval;
        else if (order == 5) return Bezier<type,N,5>::eval;
        else avro_implement;
      }
      else {
        printf("unsupported shape dimension %lu\n",number);
        avro_implement;
      }
    }
    if (category == BasisFunctionCategory_Legendre) {
      if (number == 1) {
        static int N = 1;
        if (order == 1)      return Legendre<type,N,1>::eval;
        else if (order == 2) return Legendre<type,N,2>::eval;
        else if (order == 3) return Legendre<type,N,3>::eval;
        else if (order == 4) return Legendre<type,N,4>::eval;
        else if (order == 5) return Legendre<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 2) {
        static int N = 2;
        if (order == 1)      return Legendre<type,N,1>::eval;
        else if (order == 2) return Legendre<type,N,2>::eval;
        else if (order == 3) return Legendre<type,N,3>::eval;
        else if (order == 4) return Legendre<type,N,4>::eval;
        else if (order == 5) return Legendre<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 3) {
        static int N = 3;
        if (order == 1)      return Legendre<type,N,1>::eval;
        else if (order == 2) return Legendre<type,N,2>::eval;
        else if (order == 3) return Legendre<type,N,3>::eval;
        else if (order == 4) return Legendre<type,N,4>::eval;
        else if (order == 5) return Legendre<type,N,5>::eval;
        else avro_implement;
      }
      else if (number == 4) {
        static int N = 4;
        if (order == 1)      return Legendre<type,N,1>::eval;
        else if (order == 2) return Legendre<type,N,2>::eval;
        else if (order == 3) return Legendre<type,N,3>::eval;
        else if (order == 4) return Legendre<type,N,4>::eval;
        else if (order == 5) return Legendre<type,N,5>::eval;
        else avro_implement;
      }
      else {
        printf("unsupported shape dimension %lu\n",number);
        avro_implement;
      }
    }
    return nullptr;
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
      return Lagrange<type,>::hess;
    if (category==BasisFunctionCategory_Bezier)
      return Bezier<type>::hess;
    return NULL;
  }
  */

  eval_func_ptr get_func( coord_t number, coord_t order, BasisFunctionCategory );
  eval_func_ptr get_grad( coord_t number, coord_t order, BasisFunctionCategory );
  eval_func_ptr get_hess( coord_t number, coord_t order, BasisFunctionCategory );

public:

  Basis( coord_t number , coord_t order , const BasisFunctionCategory category ) :
    category_(category),
    fptr_( get_func(number,order,category_) ),
    gptr_( get_grad(number,order,category_) ),
    hptr_( get_hess(number,order,category_) )
  {}

  void evaluate( const double* x , double* phi) const
  {
    fptr_(x,phi);
  }

  void gradient( const double* x , double* gphi ) const
  {
    gptr_(x,gphi);
  }

  void hessian( const double* x , double* hphi ) const
  {
    hptr_(x,hphi);
  }

private:

  BasisFunctionCategory category_;

  const eval_func_ptr fptr_;
  const eval_grad_ptr gptr_;
  const eval_hess_ptr hptr_;
};

} // avro

#endif
