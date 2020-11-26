// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "../MatrixS_SVD.h"
#include "../VectorS.h"
#include "../MatrixSymS.h"

#include "tools/SANSnumerics.h"

#include "tinymat/types/SurrealS.h"

#include <cmath>

namespace tinymat 
{
namespace DLA
{

//
// original source from https://github.com/ericjang/svd3
// also in LinearAlgebra/dense/static/Documentation/SVD_3x3
//

#define _gamma 5.828427124746190 // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532511287 // cos(pi/8)
#define _sstar 0.382683432365090 // sin(p/8)
#define EPSILON 1e-6

template< class T >
inline void condSwap(bool c, T &X, T &Y)
{
    // used in step 2
    T Z = X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

template< class T >
inline void condNegSwap(bool c, T &X, T &Y)
{
    // used in step 2 and 3
    T Z = -X;
    X = c ? Y : X;
    Y = c ? Z : Y;
}

template< int N, class T >
inline void quatToMat3(const T * qV, MatrixS<N,N,T>& M)
{
  T w = qV[3];
  T x = qV[0];
  T y = qV[1];
  T z = qV[2];

  T qxx = x*x;
  T qyy = y*y;
  T qzz = z*z;
  T qxz = x*z;
  T qxy = x*y;
  T qyz = y*z;
  T qwx = w*x;
  T qwy = w*y;
  T qwz = w*z;

  M(0,0) = 1 - 2*(qyy + qzz); M(0,1) = 2*(qxy - qwz);     M(0,2) = 2*(qxz + qwy);
  M(1,0) = 2*(qxy + qwz);     M(1,1) = 1 - 2*(qxx + qzz); M(1,2) = 2*(qyz - qwx);
  M(2,0) = 2*(qxz - qwy);     M(2,1) = 2*(qyz + qwx);     M(2,2) = 1 - 2*(qxx + qyy);
}

template< class T >
inline void approximateGivensQuaternion(T a11, T a12, T a22, T &ch, T &sh)
{
  /*
   * Given the givens angle computed by approximateGivensAngles,
   * compute the corresponding rotation quaternion.
   */
  ch = 2.0*(a11-a22);
  sh = a12;
  bool b = _gamma*sh*sh < ch*ch;

  T w = 1.0/sqrt(ch*ch+sh*sh);
  ch = b ? w*ch : (T)_cstar;
  sh = b ? w*sh : (T)_sstar;
}

template< int M, class T >
inline void jacobiConjugation( const int x, const int y, const int z,
                               MatrixSymS<M,T>& S, T * qV)
{
  T ch, sh;
  approximateGivensQuaternion( S(0,0), S(1,0), S(1,1), ch, sh );

  T scale = ch*ch + sh*sh;
  T a = (ch*ch - sh*sh)/scale;
  T b = (2*sh*ch)/scale;

  // make temporary copy of S
  MatrixSymS<M,T> Stmp = S;

  // perform conjugation S = Q'*S*Q
  // Q already implicitly solved from a, b
  S(0,0) =  a*( a*Stmp(0,0) + b*Stmp(1,0)) + b*( a*Stmp(1,0) + b*Stmp(1,1));
  S(1,0) =  a*(-b*Stmp(0,0) + a*Stmp(1,0)) + b*(-b*Stmp(1,0) + a*Stmp(1,1));
  S(1,1) = -b*(-b*Stmp(0,0) + a*Stmp(1,0)) + a*(-b*Stmp(1,0) + a*Stmp(1,1));
  S(2,0) =  a*Stmp(2,0) + b*Stmp(2,1);
  S(2,1) = -b*Stmp(2,0) + a*Stmp(2,1);
  S(2,2) =  Stmp(2,2);

  // update cumulative rotation qV
  T tmp[3];
  tmp[0] = qV[0]*sh;
  tmp[1] = qV[1]*sh;
  tmp[2] = qV[2]*sh;
  sh *= qV[3];

  qV[0] *= ch;
  qV[1] *= ch;
  qV[2] *= ch;
  qV[3] *= ch;

  // (x,y,z) corresponds to ((0,1,2),(1,2,0),(2,0,1))
  // for (p,q) = ((0,1),(1,2),(0,2))
  qV[z] += sh;
  qV[3] -= tmp[z]; // w
  qV[x] += tmp[y];
  qV[y] -= tmp[x];

  // re-arrange matrix for next iteration
  Stmp(0,0) = S(1,1);
  Stmp(1,0) = S(2,1); Stmp(1,1) = S(2,2);
  Stmp(2,0) = S(1,0); Stmp(2,1) = S(2,0); Stmp(2,2) = S(0,0);

  S = Stmp;
}

template< class T >
inline T dist2(T x, T y, T z)
{
    return x*x + y*y + z*z;
}

// finds transformation that diagonalizes a symmetric matrix
template< int M, class T >
inline void jacobiEigenanlysis( MatrixSymS<M,T>& S, T * qV)
{
  qV[3]=1; qV[0]=0; qV[1]=0; qV[2]=0; // follow same indexing convention as GLM

  for (int i=0;i<4;i++)
  {
    // we wish to eliminate the maximum off-diagonal element
    // on every iteration, but cycling over all 3 possible rotations
    // in fixed order (p,q) = (1,2) , (2,3), (1,3) still retains
    // asymptotic convergence
    jacobiConjugation(0, 1, 2, S, qV); // p,q = 0,1
    jacobiConjugation(1, 2, 0, S, qV); // p,q = 1,2
    jacobiConjugation(2, 0, 1, S, qV); // p,q = 0,2
  }
}

template< int M, int N, class T >
inline void sortSingularValues(MatrixS<M,N,T>& B, MatrixS<N,N,T>& V)
{
  T rho1 = dist2( B(0,0), B(1,0), B(2,0) );
  T rho2 = dist2( B(0,1), B(1,1), B(1,2) );
  T rho3 = dist2( B(0,2), B(1,2), B(2,2) );

  bool c;
  c = rho1 < rho2;
  condNegSwap(c, B(0,0), B(0,1)); condNegSwap(c, V(0,0), V(0,1));
  condNegSwap(c, B(1,0), B(1,1)); condNegSwap(c, V(1,0), V(1,1));
  condNegSwap(c, B(2,0), B(2,1)); condNegSwap(c, V(2,0), V(2,1));
  condSwap(c, rho1, rho2);

  c = rho1 < rho3;
  condNegSwap(c, B(0,0), B(0,2)); condNegSwap(c, V(0,0), V(0,2));
  condNegSwap(c, B(1,0), B(1,2)); condNegSwap(c, V(1,0), V(1,2));
  condNegSwap(c, B(2,0), B(2,2)); condNegSwap(c, V(2,0), V(2,2));
  condSwap(c, rho1, rho3);

  c = rho2 < rho3;
  condNegSwap(c, B(0,1), B(0,2)); condNegSwap(c, V(0,1), V(0,2));
  condNegSwap(c, B(1,1), B(1,2)); condNegSwap(c, V(1,1), V(1,2));
  condNegSwap(c, B(2,1), B(2,2)); condNegSwap(c, V(2,1), V(2,2));
}

template< class T >
void QRGivensQuaternion(T a1, T a2, T &ch, T &sh)
{
  // a1 = pivot point on diagonal
  // a2 = lower triangular entry we want to annihilate
  T epsilon = (T)EPSILON;
  T rho = sqrt(a1*a1 + a2*a2);

  sh = rho > epsilon ? a2 : 0;
  ch = fabs(a1) + max(rho,epsilon);
  bool b = a1 < 0;
  condSwap(b,sh,ch);
  T w = 1./sqrt(ch*ch + sh*sh);
  ch *= w;
  sh *= w;
}

template< int M, int N, class T >
inline void QRDecomposition(MatrixS<M,N,T>& A, MatrixS<M,M,T>& Q, MatrixS<M,N,T>& R)
{
  T ch1,sh1,ch2,sh2,ch3,sh3;
  T a,b;

  // first givens rotation (ch,0,0,sh)
  QRGivensQuaternion(A(0,0), A(1,0), ch1, sh1);
  a = 1.0 - 2.0*sh1*sh1;
  b =       2.0*ch1*sh1;

  // apply A = Q' * A
  R(0,0) = a*A(0,0)+b*A(1,0); R(0,1) = a*A(0,1)+b*A(1,1);  R(0,2) = a*A(0,2)+b*A(1,2);
  R(1,0) =-b*A(0,0)+a*A(1,0); R(1,1) =-b*A(0,1)+a*A(1,1);  R(1,2) =-b*A(0,2)+a*A(1,2);
  R(2,0) =   A(2,0);          R(2,1) =   A(2,1);           R(2,2) =   A(2,2);

  // second givens rotation (ch,0,-sh,0)
  QRGivensQuaternion(R(0,0),R(2,0),ch2,sh2);
  a = 1.0 - 2.0*sh2*sh2;
  b =       2.0*ch2*sh2;

  // apply A = Q' * A;
  A(0,0) = a*R(0,0)+b*R(2,0);  A(0,1) = a*R(0,1)+b*R(2,1);  A(0,2) = a*R(0,2)+b*R(2,2);
  A(1,0) =   R(1,0);           A(1,1) =   R(1,1);           A(1,2) =   R(1,2);
  A(2,0) =-b*R(0,0)+a*R(2,0);  A(2,1) =-b*R(0,1)+a*R(2,1);  A(2,2) =-b*R(0,2)+a*R(2,2);

  // third givens rotation (ch,sh,0,0)
  QRGivensQuaternion(A(1,1),A(2,1),ch3,sh3);
  a = 1.0 - 2.0*sh3*sh3;
  b =       2.0*ch3*sh3;

  // R is now set to desired value
  R(0,0) =   A(0,0);           R(0,1) =   A(0,1);           R(0,2) =   A(0,2);
  R(1,0) = a*A(1,0)+b*A(2,0);  R(1,1) = a*A(1,1)+b*A(2,1);  R(1,2) = a*A(1,2)+b*A(2,2);
  R(2,0) =-b*A(1,0)+a*A(2,0);  R(2,1) =-b*A(1,1)+a*A(2,1);  R(2,2) =-b*A(1,2)+a*A(2,2);

  // construct the cumulative rotation Q=Q1 * Q2 * Q3
  // the number of floating point operations for three quaternion multiplications
  // is more or less comparable to the explicit form of the joined matrix.
  // certainly more memory-efficient!
  T sh12=sh1*sh1;
  T sh22=sh2*sh2;
  T sh32=sh3*sh3;

  Q(0,0) = (-1+2*sh12)*(-1+2*sh22);
  Q(0,1) = 4*ch2*ch3*(-1+2*sh12)*sh2*sh3+2*ch1*sh1*(-1+2*sh32);
  Q(0,2) = 4*ch1*ch3*sh1*sh3-2*ch2*(-1+2*sh12)*sh2*(-1+2*sh32);

  Q(1,0) = 2*ch1*sh1*(1-2*sh22);
  Q(1,1) = -8*ch1*ch2*ch3*sh1*sh2*sh3+(-1+2*sh12)*(-1+2*sh32);
  Q(1,2) = -2*ch3*sh3+4*sh1*(ch3*sh1*sh3+ch1*ch2*sh2*(-1+2*sh32));

  Q(2,0) = 2*ch2*sh2;
  Q(2,1) = 2*ch3*(1-2*sh22)*sh3;
  Q(2,2) = (-1+2*sh22)*(-1+2*sh32);
}

template< int M, int N, class T >
void
SVD( const MatrixS<M,N,T>& A, MatrixS<M,M,T>& U, VectorS<MIN(M,N),T>& S, MatrixS<N,N,T>& VT )
{
  // normal equations matrix
  MatrixSymS<M,T> ATA = Transpose(A)*A;

  MatrixS<N,N,T> V = Identity();

  T qV[4]; // quaternion representation of V

  // symmetric eigenalysis
  jacobiEigenanlysis(ATA, qV);

  quatToMat3(qV,V);

  MatrixS<M,N,T> B = A*V;

  // sort singular values and find V
  sortSingularValues(B, V);

  MatrixS<M,N,T> Smat;

  // QR decomposition
  QRDecomposition(B, U, Smat);

  //Copy singular values
  for (int i = 0; i < 3; i++)
  {
    if (Smat(i,i) < 0)
    {
      Smat(i,i) = -Smat(i,i);
      V(0,i) = -V(0,i);
      V(1,i) = -V(1,i);
      V(2,i) = -V(2,i);
    }
    S[i] = Smat(i,i);
  }

  //Need to return the transpose of V
  VT = Transpose(V);
}


#define INSTANTIATE_SVD(T) \
template void SVD<3,3,T>(const MatrixS<3,3,T>& A, MatrixS<3,3,T>& U, VectorS<3,T>& S, MatrixS<3,3,T>& VT);

INSTANTIATE_SVD(Real)
//INSTANTIATE_SVD(SurrealS<1>)

}
}
