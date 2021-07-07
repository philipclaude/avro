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

namespace avro
{

template<typename type,int number,int order> class Bernstein;
template<typename type,int number,int order> class Lagrange;
template<typename type,int number,int order> class Legendre;

class Simplex;

enum BasisFunctionCategory
{
  BasisFunctionCategory_Bernstein,
  BasisFunctionCategory_Lagrange,
  BasisFunctionCategory_Legendre,
  BasisFunctionCategory_None
};

template<int N,int P>
class Bernstein<Simplex,N,P> {
public:
  static void eval( const double* x , double* phi );
  static void grad( const double* x , double* gphi );
  static void hess( const double* x , double* hphi );
};

template<int N,int P>
class Lagrange<Simplex,N,P> {
public:
  static void eval( const double* x , double* phi );
  static void grad( const double* x , double* gphi );
  static void hess( const double* x , double* hphi );
};

template<int N,int P>
class Legendre<Simplex,N,P> {
public:
  static void eval( const double* x , double* phi );
  static void grad( const double* x , double* gphi );
  static void hess( const double* x , double* hphi );
};

template<typename type>
class Basis {

private:

  typedef void (*eval_func_ptr)(const double*,double*);
  typedef void (*eval_grad_ptr)(const double*,double*);
  typedef void (*eval_hess_ptr)(const double*,double*);

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

  void evaluate( const double* x , double* phi) const {
    fptr_(x,phi);
  }

  void gradient( const double* x , double* gphi ) const {
    gptr_(x,gphi);
  }

  void hessian( const double* x , double* hphi ) const {
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
