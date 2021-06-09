// avro: Adaptive Voronoi Remesher
// Copyright 2017-2021, Philip Claude Caplan
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef AVRO_NUMERICS_DENSMAT_H_
#define AVRO_NUMERICS_DENSMAT_H_

#include "common/types.h"

#include <vector>
#include <string>

namespace avro
{

namespace numerics
{

template<typename type> class symatd;

#define MATRIX_LIMIT 10000

template <typename type>
class matd
{

	public:
		matd();
		matd(const int,const int);
    matd( const symatd<type>& T );
		~matd();

		/* Variables */
		int  m,n;     /* m rows by n columns */
		std::vector<type> data;

		void copy( const matd<type>& A , bool transpose=false );

		/* Functions */
		void allocate();
		void display() const;
		void zeros();
		void eye();
		void solveLU(type*);
		void solveLU(matd<type>&);
		std::vector<double> eig( matd<type>& q );
		void multiply(const matd<type>&,matd<type>&) const;
		void multiply(const type*,type*) const;
		void print();
    int solve(type*);
		void inv( matd<type>& Minv ) const;
    type det();

    void rankone( std::vector<type>& x );
    void scale( const type& alpha );
    void add( const matd<type>& X );

		matd<type> transpose() const;
		matd<type> inverse() const;

		/*\
		 *
		 * Call LAPACK to solve (for M) M*X = Y by solving X^t * M^t = Y^t
		 *
		 * in: X matrix, reverse flag
		 * out:
		 * in/out: Y matrix
		 *
		\*/
		int mrdivide(const matd<type>& Y , matd<type>& M );

		int kernel( matd<type>& K );

    void column( index_t j , std::vector<type>& col ) const;

		void range( matd<type>& U ) const;

		/* Operator overloading */
		type& operator() (const unsigned i,const unsigned j);
		type  operator() (const unsigned i,const unsigned j) const;
		type& val(const unsigned i,const unsigned j);

    void tomatlab( const std::string& var ) const;

		matd<type> operator* (const matd<type>& B) const;
		matd<type> operator* (const symatd<type>& B) const;

		//matd<type> operator= (const matd<real>& B);
		//matd<type> operator= (const symatd<real>& B);

	private:
		void assertInRange(const int i,const int j);
};

} // numerics

} // avro

#endif
