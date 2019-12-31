#include "numerics/lapack.h"
#include "numerics/linear_algebra.h"

#include <algorithm>
#include <vector>

namespace avro
{

namespace numerics
{

real_t
eps( const real_t& x )
{
  // this is like matlab's eps function
  return ::nextafter( x , x +1.0f ) -x;
}

template<>
int
kernel( const MatrixD<real_t>& A , MatrixD<real_t>& K )
{
  int info;

	char jobu = 'A';
	char jobvt = 'A';

  int m = A.m();
  int n = A.n();

	int M = m;
	int N = n;

	int lda = m;
	int ldu = m;
	int ldvt = n;
	std::vector<double> S( m*n , 0. );
	std::vector<double> U( m*m , 0. );
	std::vector<double> VT( n*n , 0. );

	int lwork = std::max( 3*std::min(M,N)+std::max(M,N),5*std::min(M,N)-4 ) +10;
	std::vector<double> work( lwork );

	std::vector<double> tdata( m*n );
	for (int i=0;i<m;i++)
	for (int j=0;j<n;j++)
		tdata[j*m+i] = A(i,j);

  // perform the svd
	dgesvd_( &jobu , &jobvt , &M , &N , tdata.data() , &lda , S.data() , U.data() , &ldu , VT.data() , &ldvt , work.data() , &lwork , &info );

	// now we analyze the result to get the nullspace
	std::vector<double> s;
	if (m==0) s.push_back( S[0] );
	else if (m>0)
	{
		for (int i=0;i<m;i++)
			s.push_back( S[i] );
	}
	else s.push_back(0.);


	double maxS = *std::max_element( s.begin() , s.end() );
	double tol = std::max(m,n)*eps(maxS);
	int r = 0;
	for (int i=0;i<int(s.size());i++)
	{
		if (s[i]>tol) r++;
	}

	MatrixD<real_t> v(N,N);
	for (int i=0;i<n;i++)
	for (int j=0;j<n;j++)
		v(i,j) = VT[j*n+i];

	// save the result
	//K.m = n;
	//K.n = n -r;
	//K.allocate();
  K.resize(n,n-r);
	for (int i=0;i<n;i++)
	for (int j=r;j<n;j++)
		K(i,j-r) = v(j,i);

	return info;
}

} // numerics

} // avro
