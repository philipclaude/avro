// avro: Adaptive Voronoi Remesher
// Copyright 2017-2021, Philip Claude Caplan
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_NUMERICS_symatd_H_
#define AVRO_NUMERICS_symatd_H_

#include "common/error.h"
#include "common/types.h"

#include "numerics/matd.h"

#include <vector>
#include <string>
#include <utility>

namespace avro
{

namespace numerics
{

template<typename type>
class symatd
{

public:
	symatd() {}
	symatd( const coord_t _n );
	symatd( const std::vector<type>& _data );
	symatd( const std::vector<type>& lambda , matd<type>& q );
	symatd( std::pair< std::vector<type>,matd<type> >& decomp );
  symatd( const matd<type>& M );
	template <typename typeR> symatd( const symatd<typeR>& M );

  void fromeig( const std::vector<type>& lambda , matd<type>& q );

	index_t nb() const
	{
		return n_*(n_ +1)/2;
	}

	coord_t n() const
	{
		return n_;
	}

	void copy( const symatd& T )
	{
		avro_assert( n_==T.n() );
		for (index_t k=0;k<nb();k++)
			data_[k] = T.data(k);
	}

	void zero()
	{
		for (index_t i=0;i<nb();i++)
			data_[i] = 0.;
	}

	void identity()
	{
		zero();
		for (coord_t i=0;i<n_;i++)
			operator()(i,i) = 1.;
	}

	void scale( const type alpha )
	{
		for (index_t k=0;k<nb();k++)
			data_[k] *= alpha;
	}

	type& data( const index_t k )
	{
		return data_[k];
	}

	type data( const index_t k ) const
	{
		return data_[k];
	}

	type& operator() ( const index_t i, const index_t j )
	{
		return (i>j) ? data_[i*(i+1)/2 +j] : data_[j*(j +1)/2 +i];
	}

	type operator() ( const index_t i, const index_t j ) const
	{
		return (i>j) ? data_[i*(i+1)/2 +j] : data_[j*(j +1)/2 +i];
	}

	symatd operator+() const { return *this; }
	symatd operator+( const symatd& U )
	{
		symatd V(U.n());
		for (coord_t i=0;i<U.nb();i++)
			V.data(i) = data_[i] +U.data(i);
		return V;
	}

	symatd operator*( const type a ) const
	{
		symatd C(n());
		for (index_t i=0;i<nb();i++)
			C.data(i) = a*data_[i];
		return C;
	}

	symatd<type> sandwich( const symatd<type>& B ) const;

	template<typename Method>
	void interpolate( const std::vector<type>& alpha ,
										const std::vector<symatd>& tensors )
	{
		Method::interpolate(alpha,tensors,*this);
	}

  void intersect( const symatd& M2 , symatd& T );

	symatd exp() const;
	symatd log() const;
	symatd inv() const;
	symatd pow( real_t p ) const;
	symatd sqrt() const;

	// eigenvalues and eigenvectors
	std::pair< std::vector<type>,matd<type> > eig() const;

	type quadratic_form( const type* x ) const;
  //type quadraticFormReal( real_t* x ) const;
	type determinant() const;
	void display( const std::string& title=std::string() ) const;
	void matlabize( const std::string& title ) const;
	void forunit( const std::string& title ) const;

  bool check();

private:
	coord_t n_;
	std::vector<type> data_;

	std::pair< std::vector<type>,matd<type> > __eig2__() const;
	std::pair< std::vector<type>,matd<type> > __eigivens__() const;

};

template<typename type>
class LogEuclidean
{
public:
	static void interpolate( const std::vector<type>& alpha ,
													 const std::vector<symatd<type>>& tensors , symatd<type>& T );
};

template<typename type>
class AffineInvariant
{
	static void interpolate( const std::vector<type>& alpha ,
													 const std::vector<symatd<type>>& tensors , symatd<type>& T );
};

} // numerics

} // avro

#endif
