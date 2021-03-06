//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_NUMERICS_symd_H_
#define AVRO_NUMERICS_symd_H_

#include "common/error.h"
#include "avro_types.h"

#include "numerics/mat.h"
#include "numerics/vec.h"

#include <vector>
#include <string>
#include <utility>

namespace avro
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

  void from_eig( const vecd<type>& lambda , const matd<type>& q );

	index_t nb() const { return n_*(n_ +1)/2; }
	coord_t n() const { return n_; }
	coord_t m() const { return n_; }

	void copy( const symd& T ) {
		avro_assert( n_ == T.n() );
		for (index_t k = 0; k < nb(); k++)
			data_[k] = T.data(k);
	}

	void set( const matd<type>& M );

	void zero() {
		for (index_t i=0;i<nb();i++)
			data_[i] = 0.;
	}

	void identity() {
		zero();
		for (coord_t i=0;i<n_;i++)
			operator()(i,i) = 1.;
	}

	type& data( const index_t k ) { return data_[k]; }
	type data( const index_t k ) const { return data_[k]; }

	inline type& operator() ( const index_t i, const index_t j ) {
		return (i>j) ? data_[i*(i+1)/2 +j] : data_[j*(j +1)/2 +i];
	}

	inline type operator() ( const index_t i, const index_t j ) const {
		return (i>j) ? data_[i*(i+1)/2 +j] : data_[j*(j +1)/2 +i];
	}

	symd operator+() const { return *this; }
	symd operator+( const symd& U ) {
		symd V(U.n());
		for (coord_t i=0;i<U.nb();i++)
			V.data(i) = data_[i] +U.data(i);
		return V;
	}

	symd operator*( const type a ) const {
		symd C(n());
		for (index_t i=0;i<nb();i++)
			C.data(i) = a*data_[i];
		return C;
	}

	symd& operator=( const matd<type>& A ) {
		// assumes A is already symmetric
		avro_assert( A.m() == A.n() );
		allocate( A.m() );
		for (index_t i = 0; i < n_; i++)
		for (index_t j = 0; j < n_; j++)
			(*this)(i,j) = A(i,j);
		return *this;
	}

	symd& operator=( const real_t& a ) {
		std::fill( data_.begin() , data_.end() , a );
		return *this;
	}

	symd<type> sandwich( const symd<type>& B ) const;

	// eigenvalues and eigenvectors
	std::pair< vecd<type>,matd<type> > eig() const;

	void display( const std::string& title=std::string() ) const;
	void for_matlab( const std::string& title ) const;
	void for_unit( const std::string& title ) const;
	void dump() const { display(); }

private:
	coord_t n_;
	std::vector<type> data_;

	std::pair< vecd<type>,matd<type> > __eig2__() const;
	std::pair< vecd<type>,matd<type> > __eigivens__() const;
};

template<typename R, typename S> symd< typename result_of<R,S>::type > operator*( const symd<R>& A , const symd<S>& B );

} // avro

#endif
