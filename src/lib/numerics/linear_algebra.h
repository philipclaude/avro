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

#include "common/error.h"

#include "numerics/mat.h"

#include "numerics/determinant.h"

#include <cmath>

namespace avro
{

namespace numerics
{

template<index_t M,typename T> mats<M,M,T> inverse( const mats<M,M,T>& A );
template<index_t M,typename T> T det( const mats<M,M,T>& A );
template<index_t M,typename T> T trace( const mats<M,M,T>& A );
template<typename T> matd<T> transpose( const matd<T>& A );

template<typename type> symd<type> expm( const symd<type>& M ) { return M.exp(); }
template<typename type> symd<type> logm( const symd<type>& M ) { return M.log(); }
template<typename type> symd<type> powm( const symd<type>& M , real_t p ) { return M.pow(p); }
template<typename type> symd<type> sqrtm( const symd<type>& M ) { return M.sqrt(); }

template<typename type> type det( const symd<type>& M );

template<typename T> void solveLUP( const matd<T>& A , const vecd<T>& b , vecd<T>& x );
template<typename T> void inverseLUP( const matd<T>& A , matd<T>& Ainv );

template<typename type> int kernel( const matd<type>& A , matd<type>& K );

template<typename type> int range( const matd<type>& A , matd<type>& U );

template<typename T> matd<T> diag( const vecd<T>& d );

template<typename T>
inline matd<T>
inverse( matd<T>& M )
{
	avro_assert( M.m() == M.n() );
	T idetM = 1./det(M);

  matd<T> Minv(M.m(),M.n());

	if (M.n()==1)
	{
		Minv(0,0) = idetM;
	}
	else if (M.n()==2)
	{
		Minv(0,0) =  M(1,1)*idetM;
		Minv(0,1) = -M(0,1)*idetM;
		Minv(1,0) = -M(1,0)*idetM;
		Minv(1,1) =  M(0,0)*idetM;
	}
	else if (M.n()==3)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2);
		Minv(0,0) = (a2_2*a3_3 -a2_3*a3_2)*idetM;
		Minv(0,1) = (a1_3*a3_2 -a1_2*a3_3)*idetM;
		Minv(0,2) = (a1_2*a2_3 -a1_3*a2_2)*idetM;
		Minv(1,0) = (a2_3*a3_1 -a2_1*a3_3)*idetM;
		Minv(1,1) = (a1_1*a3_3 -a1_3*a3_1)*idetM;
		Minv(1,2) = (a1_3*a2_1 -a1_1*a2_3)*idetM;
		Minv(2,0) = (a2_1*a3_2 -a2_2*a3_1)*idetM;
		Minv(2,1) = (a1_2*a3_1 -a1_1*a3_2)*idetM;
		Minv(2,2) = (a1_1*a2_2 -a1_2*a2_1)*idetM;
	}
	else if (M.n()==4)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2); T a1_4 = M(0,3);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2); T a2_4 = M(1,3);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2); T a3_4 = M(2,3);
		T a4_1 = M(3,0); T a4_2 = M(3,1); T a4_3 = M(3,2); T a4_4 = M(3,3);

		Minv(0,0) = (a2_2*a3_3*a4_4-a2_2*a3_4*a4_3-a2_3*a3_2*a4_4+a2_3*a3_4*a4_2+a2_4*a3_2*a4_3-a2_4*a3_3*a4_2)*idetM;
		Minv(0,1) = (-a1_2*a3_3*a4_4+a1_2*a3_4*a4_3+a1_3*a3_2*a4_4-a1_3*a3_4*a4_2-a1_4*a3_2*a4_3+a1_4*a3_3*a4_2)*idetM;
		Minv(0,2) = (a1_2*a2_3*a4_4-a1_2*a2_4*a4_3-a1_3*a2_2*a4_4+a1_3*a2_4*a4_2+a1_4*a2_2*a4_3-a1_4*a2_3*a4_2)*idetM;
		Minv(0,3) = (-a1_2*a2_3*a3_4+a1_2*a2_4*a3_3+a1_3*a2_2*a3_4-a1_3*a2_4*a3_2-a1_4*a2_2*a3_3+a1_4*a2_3*a3_2)*idetM;
		Minv(1,0) = (-a2_1*a3_3*a4_4+a2_1*a3_4*a4_3+a2_3*a3_1*a4_4-a2_3*a3_4*a4_1-a2_4*a3_1*a4_3+a2_4*a3_3*a4_1)*idetM;
		Minv(1,1) = (a1_1*a3_3*a4_4-a1_1*a3_4*a4_3-a1_3*a3_1*a4_4+a1_3*a3_4*a4_1+a1_4*a3_1*a4_3-a1_4*a3_3*a4_1)*idetM;
		Minv(1,2) = (-a1_1*a2_3*a4_4+a1_1*a2_4*a4_3+a1_3*a2_1*a4_4-a1_3*a2_4*a4_1-a1_4*a2_1*a4_3+a1_4*a2_3*a4_1)*idetM;
		Minv(1,3) = (a1_1*a2_3*a3_4-a1_1*a2_4*a3_3-a1_3*a2_1*a3_4+a1_3*a2_4*a3_1+a1_4*a2_1*a3_3-a1_4*a2_3*a3_1)*idetM;
		Minv(2,0) = (a2_1*a3_2*a4_4-a2_1*a3_4*a4_2-a2_2*a3_1*a4_4+a2_2*a3_4*a4_1+a2_4*a3_1*a4_2-a2_4*a3_2*a4_1)*idetM;
		Minv(2,1) = (-a1_1*a3_2*a4_4+a1_1*a3_4*a4_2+a1_2*a3_1*a4_4-a1_2*a3_4*a4_1-a1_4*a3_1*a4_2+a1_4*a3_2*a4_1)*idetM;
		Minv(2,2) = (a1_1*a2_2*a4_4-a1_1*a2_4*a4_2-a1_2*a2_1*a4_4+a1_2*a2_4*a4_1+a1_4*a2_1*a4_2-a1_4*a2_2*a4_1)*idetM;
		Minv(2,3) = (-a1_1*a2_2*a3_4+a1_1*a2_4*a3_2+a1_2*a2_1*a3_4-a1_2*a2_4*a3_1-a1_4*a2_1*a3_2+a1_4*a2_2*a3_1)*idetM;
		Minv(3,0) = (-a2_1*a3_2*a4_3+a2_1*a3_3*a4_2+a2_2*a3_1*a4_3-a2_2*a3_3*a4_1-a2_3*a3_1*a4_2+a2_3*a3_2*a4_1)*idetM;
		Minv(3,1) = (a1_1*a3_2*a4_3-a1_1*a3_3*a4_2-a1_2*a3_1*a4_3+a1_2*a3_3*a4_1+a1_3*a3_1*a4_2-a1_3*a3_2*a4_1)*idetM;
		Minv(3,2) = (-a1_1*a2_2*a4_3+a1_1*a2_3*a4_2+a1_2*a2_1*a4_3-a1_2*a2_3*a4_1-a1_3*a2_1*a4_2+a1_3*a2_2*a4_1)*idetM;
		Minv(3,3) = (a1_1*a2_2*a3_3-a1_1*a2_3*a3_2-a1_2*a2_1*a3_3+a1_2*a2_3*a3_1+a1_3*a2_1*a3_2-a1_3*a2_2*a3_1)*idetM;
	}
	else
		avro_implement;
  return Minv;
}

template<typename T>
inline symd<T>
inverse( const symd<T>& M )
{
	T idetM = 1./det(M);

  symd<T> Minv(M.m(),M.n());

	if (M.n()==1)
	{
		Minv(0,0) = idetM;
	}
	else if (M.n()==2)
	{
		Minv(0,0) =  M(1,1)*idetM;
		Minv(0,1) = -M(0,1)*idetM;
		Minv(1,0) = -M(1,0)*idetM;
		Minv(1,1) =  M(0,0)*idetM;
	}
	else if (M.n()==3)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2);
		Minv(0,0) = (a2_2*a3_3 -a2_3*a3_2)*idetM;
		Minv(0,1) = (a1_3*a3_2 -a1_2*a3_3)*idetM;
		Minv(0,2) = (a1_2*a2_3 -a1_3*a2_2)*idetM;
		Minv(1,0) = (a2_3*a3_1 -a2_1*a3_3)*idetM;
		Minv(1,1) = (a1_1*a3_3 -a1_3*a3_1)*idetM;
		Minv(1,2) = (a1_3*a2_1 -a1_1*a2_3)*idetM;
		Minv(2,0) = (a2_1*a3_2 -a2_2*a3_1)*idetM;
		Minv(2,1) = (a1_2*a3_1 -a1_1*a3_2)*idetM;
		Minv(2,2) = (a1_1*a2_2 -a1_2*a2_1)*idetM;
	}
	else if (M.n()==4)
	{
		T a1_1 = M(0,0); T a1_2 = M(0,1); T a1_3 = M(0,2); T a1_4 = M(0,3);
		T a2_1 = M(1,0); T a2_2 = M(1,1); T a2_3 = M(1,2); T a2_4 = M(1,3);
		T a3_1 = M(2,0); T a3_2 = M(2,1); T a3_3 = M(2,2); T a3_4 = M(2,3);
		T a4_1 = M(3,0); T a4_2 = M(3,1); T a4_3 = M(3,2); T a4_4 = M(3,3);

		Minv(0,0) = (a2_2*a3_3*a4_4-a2_2*a3_4*a4_3-a2_3*a3_2*a4_4+a2_3*a3_4*a4_2+a2_4*a3_2*a4_3-a2_4*a3_3*a4_2)*idetM;
		Minv(0,1) = (-a1_2*a3_3*a4_4+a1_2*a3_4*a4_3+a1_3*a3_2*a4_4-a1_3*a3_4*a4_2-a1_4*a3_2*a4_3+a1_4*a3_3*a4_2)*idetM;
		Minv(0,2) = (a1_2*a2_3*a4_4-a1_2*a2_4*a4_3-a1_3*a2_2*a4_4+a1_3*a2_4*a4_2+a1_4*a2_2*a4_3-a1_4*a2_3*a4_2)*idetM;
		Minv(0,3) = (-a1_2*a2_3*a3_4+a1_2*a2_4*a3_3+a1_3*a2_2*a3_4-a1_3*a2_4*a3_2-a1_4*a2_2*a3_3+a1_4*a2_3*a3_2)*idetM;
		Minv(1,0) = (-a2_1*a3_3*a4_4+a2_1*a3_4*a4_3+a2_3*a3_1*a4_4-a2_3*a3_4*a4_1-a2_4*a3_1*a4_3+a2_4*a3_3*a4_1)*idetM;
		Minv(1,1) = (a1_1*a3_3*a4_4-a1_1*a3_4*a4_3-a1_3*a3_1*a4_4+a1_3*a3_4*a4_1+a1_4*a3_1*a4_3-a1_4*a3_3*a4_1)*idetM;
		Minv(1,2) = (-a1_1*a2_3*a4_4+a1_1*a2_4*a4_3+a1_3*a2_1*a4_4-a1_3*a2_4*a4_1-a1_4*a2_1*a4_3+a1_4*a2_3*a4_1)*idetM;
		Minv(1,3) = (a1_1*a2_3*a3_4-a1_1*a2_4*a3_3-a1_3*a2_1*a3_4+a1_3*a2_4*a3_1+a1_4*a2_1*a3_3-a1_4*a2_3*a3_1)*idetM;
		Minv(2,0) = (a2_1*a3_2*a4_4-a2_1*a3_4*a4_2-a2_2*a3_1*a4_4+a2_2*a3_4*a4_1+a2_4*a3_1*a4_2-a2_4*a3_2*a4_1)*idetM;
		Minv(2,1) = (-a1_1*a3_2*a4_4+a1_1*a3_4*a4_2+a1_2*a3_1*a4_4-a1_2*a3_4*a4_1-a1_4*a3_1*a4_2+a1_4*a3_2*a4_1)*idetM;
		Minv(2,2) = (a1_1*a2_2*a4_4-a1_1*a2_4*a4_2-a1_2*a2_1*a4_4+a1_2*a2_4*a4_1+a1_4*a2_1*a4_2-a1_4*a2_2*a4_1)*idetM;
		Minv(2,3) = (-a1_1*a2_2*a3_4+a1_1*a2_4*a3_2+a1_2*a2_1*a3_4-a1_2*a2_4*a3_1-a1_4*a2_1*a3_2+a1_4*a2_2*a3_1)*idetM;
		Minv(3,0) = (-a2_1*a3_2*a4_3+a2_1*a3_3*a4_2+a2_2*a3_1*a4_3-a2_2*a3_3*a4_1-a2_3*a3_1*a4_2+a2_3*a3_2*a4_1)*idetM;
		Minv(3,1) = (a1_1*a3_2*a4_3-a1_1*a3_3*a4_2-a1_2*a3_1*a4_3+a1_2*a3_3*a4_1+a1_3*a3_1*a4_2-a1_3*a3_2*a4_1)*idetM;
		Minv(3,2) = (-a1_1*a2_2*a4_3+a1_1*a2_3*a4_2+a1_2*a2_1*a4_3-a1_2*a2_3*a4_1-a1_3*a2_1*a4_2+a1_3*a2_2*a4_1)*idetM;
		Minv(3,3) = (a1_1*a2_2*a3_3-a1_1*a2_3*a3_2-a1_2*a2_1*a3_3+a1_2*a2_3*a3_1+a1_3*a2_1*a3_2-a1_3*a2_2*a3_1)*idetM;
	}
	else
		avro_implement;
  return Minv;
}

template<typename T>
inline
std::pair< vecd<T> , matd<T> >
eig( const symd<T>& m ) {
  return m.eig();
}

template<typename T>
inline
void
eig( const symd<T>& m , vecd<T>& L , matd<T>& Q ) {
  std::pair< vecd<T> , matd<T> > decomp = m.eig();
  L.set( decomp.first );
  Q.set( decomp.second );
}

} // numerics

} // avro

#endif
