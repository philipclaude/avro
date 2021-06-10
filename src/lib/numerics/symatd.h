// avro: Adaptive Voronoi Remesher
// Copyright 2017-2021, Philip Claude Caplan
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_NUMERICS_symd_H_
#define AVRO_NUMERICS_symd_H_

#include "common/error.h"
#include "common/types.h"

#include "numerics/mat.h"
#include "numerics/vec.h"

#include <vector>
#include <string>
#include <utility>

namespace avro
{

namespace numerics
{

template<typename type>
class symd
{

public:
	symd() {}
	symd( const coord_t _n );
	symd( const std::vector<type>& _data );
	symd( const vecd<type>& lambda , const matd<type>& q );
	symd( const std::pair< vecd<type>,matd<type> >& decomp );
  symd( const matd<type>& M );
	template <typename typeR> symd( const symd<typeR>& M );

	void allocate( const coord_t n ) {
		n_ = n;
		data_.resize( nb() );
	}

  void fromeig( const vecd<type>& lambda , const matd<type>& q );

	index_t nb() const { return n_*(n_ +1)/2; }
	coord_t n() const { return n_; }
	coord_t m() const { return n_; }

	void copy( const symd& T )
	{
		avro_assert( n_==T.n() );
		for (index_t k=0;k<nb();k++)
			data_[k] = T.data(k);
	}

	void zero() {
		for (index_t i=0;i<nb();i++)
			data_[i] = 0.;
	}

	void identity() {
		zero();
		for (coord_t i=0;i<n_;i++)
			operator()(i,i) = 1.;
	}

	void scale( const type alpha ) {
		for (index_t k=0;k<nb();k++)
			data_[k] *= alpha;
	}

	type& data( const index_t k ) { return data_[k]; }

	type data( const index_t k ) const { return data_[k]; }

	type& operator() ( const index_t i, const index_t j )
	{
		return (i>j) ? data_[i*(i+1)/2 +j] : data_[j*(j +1)/2 +i];
	}

	type operator() ( const index_t i, const index_t j ) const
	{
		return (i>j) ? data_[i*(i+1)/2 +j] : data_[j*(j +1)/2 +i];
	}

	symd operator+() const { return *this; }
	symd operator+( const symd& U )
	{
		symd V(U.n());
		for (coord_t i=0;i<U.nb();i++)
			V.data(i) = data_[i] +U.data(i);
		return V;
	}

	symd operator*( const type a ) const
	{
		symd C(n());
		for (index_t i=0;i<nb();i++)
			C.data(i) = a*data_[i];
		return C;
	}

	symd& operator=( const matd<type>& A ) {
		avro_implement;
		return *this;
	}

	symd& operator=( const real_t& a ) {
		std::fill( data_.begin() , data_.end() , a );
		return *this;
	}

	template<typename Method>
	void interpolate( const std::vector<type>& alpha ,
										const std::vector<symd>& tensors )
	{
		Method::interpolate(alpha,tensors,*this);
	}

	symd exp() const;
	symd log() const;
	symd inv() const;
	symd pow( real_t p ) const;
	symd sqrt() const;

	// eigenvalues and eigenvectors
	std::pair< vecd<type>,matd<type> > eig() const;

	type quadratic_form( const type* x ) const;
  //type quadraticFormReal( real_t* x ) const;
	type determinant() const;
	void display( const std::string& title=std::string() ) const;
	void matlabize( const std::string& title ) const;
	void forunit( const std::string& title ) const;

	void dump() const;

  bool check();

private:
	coord_t n_;
	std::vector<type> data_;

	std::pair< vecd<type>,matd<type> > __eig2__() const;
	std::pair< vecd<type>,matd<type> > __eigivens__() const;

};

template<typename type>
class LogEuclidean
{
public:
	static void interpolate( const std::vector<type>& alpha ,
													 const std::vector<symd<type>>& tensors , symd<type>& T );
};

template<typename type>
class AffineInvariant
{
	static void interpolate( const std::vector<type>& alpha ,
													 const std::vector<symd<type>>& tensors , symd<type>& T );
};

template<typename type> symd<type> expm( const symd<type>& M );
template<typename type> symd<type> logm( const symd<type>& M );
template<typename type> symd<type> powm( const symd<type>& M , real_t p );
template<typename type> symd<type> sqrtm( const symd<type>& M );

template<typename type> type det( const symd<type>& M );

template<typename R, typename S>
symd< typename result_of<R,S>::type> operator*( const symd<R>& A , const symd<S>& B );

template<typename type> symd<type> inverse2( const symd<type>& M );


} // numerics

} // avro

#endif
