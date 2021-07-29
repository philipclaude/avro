//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_LINEAR_ALGEBRA_H_
#define avro_LIB_LINEAR_ALGEBRA_H_

#include <cmath>
#include "avro_types.h"

namespace avro
{

// forward declarations
template<index_t M, index_t N,typename T> class mats;
template<typename type> class matd;
template <typename T> class symd;
template <typename T> class vecd;

namespace numerics
{

/*\
 * =============================================================================
 *
 * mats
 *
 * =============================================================================
\*/
template<index_t M, index_t N,typename T> mats<N,M,T> transpose( const mats<M,N,T>& A );
template<index_t M,typename T> T det( const mats<M,M,T>& A );
template<index_t M,typename T> T trace( const mats<M,M,T>& A );
template<index_t M,typename T> mats<M,M,T> inverse( const mats<M,M,T>& A );

/*\
 * =============================================================================
 *
 * matd
 *
 * =============================================================================
\*/
template<typename T> matd<T> transpose( const matd<T>& A );
template<typename type> int  kernel( const matd<type>& A , matd<type>& K );
template<typename type> int  range( const matd<type>& A , matd<type>& U );
template<typename T> matd<T> diag( const vecd<T>& d );
template<typename T> void    solveLUP( const matd<T>& A , const vecd<T>& b , vecd<T>& x );
template<typename T> void    inverseLUP( const matd<T>& A , matd<T>& Ainv );
template<typename T> matd<T> inverse( const matd<T>& M );
template<typename T> T det(const matd<T>& A);
template<typename T> void eign( const matd<T>& A , vecd<T>& L , matd<T>& Q );
template<typename T> void eign( const matd<T>& A , vecd<T>& L );

/*\
 * =============================================================================
 *
 * symd
 *
 * =============================================================================
\*/
template<typename type> symd<type> expm( const symd<type>& M );
template<typename type> symd<type> logm( const symd<type>& M );
template<typename type> symd<type> powm( const symd<type>& M , real_t p );
template<typename type> symd<type> sqrtm( const symd<type>& M );
template<typename type> type       det( const symd<type>& M );
template<typename T> symd<T>       inverse( const symd<T>& M );
template<typename T> void eig( const symd<T>& m , vecd<T>& L , matd<T>& Q );
template<typename T> std::pair< vecd<T> , matd<T> > eig( const symd<T>& m );
template<typename T> std::pair< vecd<T> , matd<T> > eign( const symd<T>& m);
template<typename T> symd<T> interp( const std::vector<real_t>& alpha , const std::vector<symd<T>>& tensors);
template<typename type> type quadratic_form( const symd<type>& M , const vecd<real_t>& e );
template<typename type> type quadratic_form( const symd<type>& M , const type* e );

} // numerics

} // avro

#endif
