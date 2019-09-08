// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(MATRIXSYMS_DECOMPOSE_CHOLESKY_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <cmath> //sqrt

#include "MatrixSymS_Decompose_Cholesky.h"
#include "numpack/DenseLinAlg/tools/SingularException.h"
#include "MatrixSymS.h"

//Perform a Cholesky decomposition on a symmetric positive definite matrix

namespace numpack 
{
namespace DLA
{

//-----------------------------------------------------------------------------
template< int M, class T >
void MatrixSymSCholesky<M, T>::Decompose( MatrixSymS<M, T>& A )
{
  //Perform the Cholesky decomposition
  for (int i = 0; i < M; i++)
  {
    for (int j = 0; j < (i+1); j++)
    {
      T s = 0;
      for (int k = 0; k < j; k++)
        s += A(i, k) * A(j, k);

      if (i == j)
      {
        SANS_ASSERT_MSG( A(i, i) - s >= 0, "Matrix is not SPD: A(%d,%d) = %f, s = %f", i, i, A(i, i), s);
        A(i, j) = sqrt(A(i, i) - s);
      }
      else
      {
        T denom = A(j, j);
        T numer = (A(i, j) - s);

        SANS_ASSERT_NONSINGULAR(denom);

        A(i, j) = numer/denom;
      }
    }
  }

}

#if DOCUMENTATION && 0

// http://rosettacode.org/wiki/Cholesky_decomposition

#include <stdio.h>
#include <stdlib.h>

double *cholesky(double *A, int n)
{
    double *L = (double*)calloc(n * n, sizeof(double));
    if (L == NULL)
        exit(EXIT_FAILURE);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < (i+1); j++)
        {
            double s = 0;
            for (int k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
            L[i * n + j] = (i == j) ?
                           sqrt(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }

    return L;
}

void show_matrix(double *A, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            printf("%2.5f ", A[i * n + j]);
        printf("\n");
    }
}

int main()
{
    int n = 3;
    double m1[] = {25, 15, -5,
                   15, 18,  0,
                   -5,  0, 11};
    double *c1 = cholesky(m1, n);
    show_matrix(c1, n);
    printf("\n");
    free(c1);

    n = 4;
    double m2[] = {18, 22,  54,  42,
                   22, 70,  86,  62,
                   54, 86, 174, 134,
                   42, 62, 134, 106};
    double *c2 = cholesky(m2, n);
    show_matrix(c2, n);
    free(c2);

    return 0;
}

/*
Output:
 5.00000 0.00000 0.00000
 3.00000 3.00000 0.00000
-1.00000 1.00000 3.00000

 4.24264 0.00000 0.00000 0.00000
 5.18545 6.56591 0.00000 0.00000
12.72792 3.04604 1.64974 0.00000
 9.89949 1.62455 1.84971 1.39262
*/
#endif

} //namespace DLA
} //namespace numpack 
