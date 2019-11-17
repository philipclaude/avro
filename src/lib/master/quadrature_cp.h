// no include guard, should only be included once

#include "common/error.h"

#include <algorithm>
#include <cmath> // sqrt
#include <string>
#include <iostream>
#include <limits>

void cgqf(int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[]);
void cdgqf(int nt, int kind, double alpha, double beta, double t[], double wts[]);
void parchk(int kind, int m, double alpha, double beta);
void scqf(int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[],
      double st[], int kind, double alpha, double beta, double a, double b);
void imtqlx(int n, double d[], double e[], double z[]);
double class_matrix(int kind, int m, double alpha, double beta, double aj[], double bj[]);
void sgqf(int nt, double aj[], double bj[], double zemu, double t[], double wts[]);

static int
factorial( int n )
{
  int m;
  if (n <= 1)
    m = 1;
  else
    m = n * factorial(n - 1);

  return m;
}

/**
 * calculateJacobiQuadrature calculates a Gauss-Jacobi quadrature rule.
 *
 * $\int_{0}^{1} (1 - x)^\alpha x^\beta f(x) \mathrm{~d} x \approx \sum_i w_i x_i$;
 *
 * This function exists to wrap the code of John Burkardt.
 * @see http://people.sc.fsu.edu/~jburkardt/cpp_src/jacobi_rule/jacobi_rule.html
 *
 * @param n[in]       the number of points (i.e. quadrature order) requested
 * @param alpha[in]   $\alpha$
 * @param beta[in]    $\beta$
 * @param t[out]      the $n$ quadrature points, $x_i$
 * @param wts[out]    the $n$ quadrature weights, $w_i$
 */
void
calculateJacobiQuadrature(int order, double alpha, double beta, double *t, double *wts)
{
  const int kind = 4;

  // get the weights on the [0, 1] interval
  cgqf(order, kind, alpha, beta, 0.0, 1.0, t, wts);
}

/**
 * CGQF computes knots and weights of a Gauss quadrature formula.
 *
 * Purpose:
 *   CGQF computes knots and weights of a Gauss quadrature formula.
 *
 * Discussion:
 *   The user may specify the interval (A,B).
 *
 *   Only simple knots are produced.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   16 February 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 * Various quadrature rules can be selected by this program, indicated by
 * the parameter kind, which is an integer. Options include:
 *   - 1, Legendre,             (a,b)       1.0
 *   - 2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
 *   - 3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
 *   - 4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
 *   - 5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
 *   - 6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
 *   - 7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
 *   - 8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
 *
 * Converted for SANS by Cory Frontin in 2019
 * imported into luna by Philip Caplan in 2019
 *
 * @param nt      the number of knots
 * @param kind    the type of quadrature used
 * @param alpha   floating point number for alpha, possibly required for a given quadrature
 * @param beta    floating point number for beta, possibly required for a given quadrature
 * @param a       left endpoint for interval domain
 * @param b       right endpoint for interval domain
 * @param t       the knots of the quadrature rule
 * @param wts     the weights of the quadrature rule
 *
 */
void
cgqf(int nt, int kind, double alpha, double beta, double a, double b, double t[], double wts[])
{
  int i;
  int *mlt;
  int *ndx;

  // compute the Gauss quadrature formula for default values of A and B.

  cdgqf(nt, kind, alpha, beta, t, wts);

  // prepare to scale the quadrature formula to other weight function with valid A and B.

  mlt= new int[nt];
  for (i = 0; i < nt; i++)
    mlt[i] = 1;
  ndx= new int[nt];
  for (i = 0; i < nt; i++)
    ndx[i] = i + 1;
  scqf(nt, t, mlt, wts, nt, ndx, wts, t, kind, alpha, beta, a, b);

  delete [] mlt;
  delete [] ndx;

  return;
}

/**
 * sign calculates the sign of a double.
 *
 * @param x         input
 * @return value    +1.0 if x is positive or -1.0 if x is negative
 */
double inline
sign(double x)
{
  double value;
  if (x < 0.0)
    value = -1.0;
  else
    value = 1.0;
  return value;
}

/**
 * SCQF scales a quadrature formula to a nonstandard interval.
 *
 * Purpose:
 *   SCQF scales a quadrature formula to a nonstandard interval.
 *
 * Discussion:
 *   The arrays WTS and SWTS may coincide.
 *
 *   The arrays T and ST may coincide.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   16 February 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 * @param nt[in]      the number of knots
 * @param t[in]       the original knots of the quadrature rule
 * @param mlt[in]     the multiplicity of the knots
 * @param wts[in]     the weights
 * @param nwts[in]    the number of weights
 * @param ndx[in]     used to index the array wts
 * @param swts[out]   the scaled output weights
 * @param st[out]     the scaled knows
 * @param kind[in]    the type of quadrature used
 * @param alpha[in]   the value of alpha, if needed
 * @param beta[in]    the value of beta, if needed
 * @param a[in]       the left endpoint of the interval
 * @param b[in]       the right endpoint of the interval
 *
 */
void
scqf(int nt, double t[], int mlt[], double wts[], int nwts, int ndx[], double swts[],
      double st[], int kind, double alpha, double beta, double a, double b)
{
  double al = 0;
  double be = 0;
  int i = 0;
  int k = 0;
  int lz = 0;
  double p = 0;
  double shft = 0;
  double slp = 0;
  double temp = 0;
  double tmp = 0;

  temp = std::numeric_limits<double>::epsilon();

  parchk(kind, 1, alpha, beta);

  if (kind == 1)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    al = 0.0;
    be = 0.0;
    if (fabs(b - a) <= temp)
      luna_assert_msg(false,"SCQF - Fatal error!\n|B - A| too small.");
    shft = (a + b)/2.0;
    slp = (b - a)/2.0;
  }
  else if (kind == 2)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    al = -0.5;
    be = -0.5;
    if (fabs(b - a) <= temp)
      luna_assert_msg(false,"SCQF - Fatal error!\n|B - A| too small.");
    shft = (a + b)/2.0;
    slp = (b - a)/2.0;
  }
  else if (kind == 3)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    al = alpha;
    be = alpha;
    if (fabs(b - a) <= temp)
      luna_assert_msg(false,"SCQF - Fatal error!\n|B - A| too small.");
    shft = (a + b)/2.0;
    slp = (b - a)/2.0;
  }
  else if (kind == 4)
  {
    al = alpha;
    be = beta;

    if (fabs(b - a) <= temp)
      luna_assert_msg(false,"SCQF - Fatal error!\n|B - A| too small.");
    shft = (a + b)/2.0;
    slp = (b - a)/2.0;
  }
  else if ( kind == 5 )
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    if (b <= 0.0)
      luna_assert_msg(false,"SCQF - Fatal error!\nB <= 0.");
    shft = a;
    slp = 1.0/b;
    al = alpha;
    be = 0.0;
  }
  else if (kind == 6)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    if (b <= 0.0)
      luna_assert_msg(false,"SCQF - Fatal error!\nB <= 0.");
    shft = a;
    slp = 1.0/sqrt(b);
    al = alpha;
    be = 0.0;
  }
  else if (kind == 7)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    al = alpha;
    be = 0.0;
    if (fabs(b - a) <= temp)
      luna_assert_msg(false,"SCQF - Fatal error!\n|B - A| too small.");
    shft = (a + b)/2.0;
    slp = (b - a)/2.0;
  }
  else if (kind == 8)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    if (a + b <= 0.0)
      luna_assert_msg(false,"SCQF - Fatal error!\nA + B <= 0.");
    shft = a;
    slp = a + b;
    al = alpha;
    be = beta;
  }
  else if (kind == 9)
  {
    luna_assert_msg(false,"KIND != 4 NOT UNIT TESTED.");

    al = 0.5;
    be = 0.5;
    if (fabs(b - a) <= temp)
      luna_assert_msg(false,"SCQF - Fatal error!\n|B - A| too small.");
    shft = (a + b)/2.0;
    slp = (b - a)/2.0;
  }
  else
    luna_assert_msg(false,"KIND specified is not implemented for shifts.");

  p = pow(slp, al + be + 1.0);

  for (k = 0; k < nt; k++)
  {
    st[k] = shft + slp*t[k];
    lz = abs(ndx[k]);

    if (lz != 0)
    {
      tmp = p;
      for (i = lz - 1; i <= lz - 1 + mlt[k] - 1; i++)
      {
        swts[i] = wts[i]*tmp;
        tmp = tmp*slp;
      }
    }
  }
  return;
}

/**
 * CDGQF computes a Gauss quadrature formula on the unit interval
 *
 * Purpose:
 *   CDGQF computes a Gauss quadrature formula with default A, B and simple knots.
 *
 * Discussion:
 *   This routine computes all the knots and weights of a Gauss quadrature
 *   formula with a classical weight function with default values for A and B,
 *   and only simple knots.
 *
 *   There are no moments checks and no printing is done.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   08 January 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 * Various quadrature rules can be selected by this program, indicated by
 * the parameter kind, which is an integer. Options include:
 *   - 1, Legendre,             (a,b)       1.0
 *   - 2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
 *   - 3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
 *   - 4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
 *   - 5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
 *   - 6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
 *   - 7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
 *   - 8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
 *
 * Converted for SANS by Cory Frontin in 2019
 * imported into luna by Philip Caplan in 2019
 *
 * @param nt      the number of knots
 * @param kind    the type of quadrature used
 * @param alpha   floating point number for alpha, possibly required for a given quadrature
 * @param beta    floating point number for beta, possibly required for a given quadrature
 * @param t       the knots of the quadrature rule
 * @param wts     the weights of the quadrature rule
 *
 */
void
cdgqf(int nt, int kind, double alpha, double beta, double t[], double wts[])
{
  double *aj;
  double *bj;
  double zemu;

  parchk( kind, 2*nt, alpha, beta);

  // get the Jacobi matrix and zero-th moment.

  aj = new double[nt];
  bj = new double[nt];

  zemu = class_matrix(kind, nt, alpha, beta, aj, bj);

  // compute the knots and weights.

  sgqf(nt, aj, bj, zemu, t, wts);

  delete [] aj;
  delete [] bj;

  return;
}

/**
 * SGQF computes knots and weights of a Gauss Quadrature formula.
 *
 * Purpose:
 *   SGQF computes knots and weights of a Gauss Quadrature formula.
 *
 * Discussion:
 *   This routine computes all the knots and weights of a Gauss quadrature
 *   formula with simple knots from the Jacobi matrix and the zero-th
 *   moment of the weight function, using the Golub-Welsch technique.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   08 January 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 * @param nt              the number of knots
 * @param aj[in]          the diagonal of the Jacobi matrix
 * @param bj[in]          the subdiagonal of the Jacobi matrix
 * @param bj[out]         the overwritten subdiagonal of the Jacobi matrix
 * @param zemu[in]        the zeroth moment of the weight function
 * @oaram t[out]          the knots of the quadrature rule
 * @oaram wts[out]        the weights of the quadrature rule
 *
 */
void
sgqf(int nt, double aj[], double bj[], double zemu, double t[], double wts[])
{
  int i;

  // exit if the zero-th moment is not positive.

  if (zemu <= 0.0)
    luna_assert_msg(false,"SGQF error: zemu is not positive.");

  // set up vectors for IMTQLX.

  for (i = 0; i < nt; i++)
    t[i] = aj[i];
  wts[0] = sqrt(zemu);
  for (i = 1; i < nt; i++)
    wts[i] = 0.0;

  imtqlx(nt, t, bj, wts);

  for (i = 0; i < nt; i++)
    wts[i] = wts[i]*wts[i];

  return;
}

/**
 * CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
 *
 * Purpose:
 *   CLASS_MATRIX computes the Jacobi matrix for a quadrature rule.
 *
 * Discussion:
 *   This routine computes the diagonal AJ and sub-diagonal BJ
 *   elements of the order M tridiagonal symmetric Jacobi matrix
 *   associated with the polynomials orthogonal with respect to
 *   the weight function specified by KIND.
 *
 *   For weight functions 1-7, M elements are defined in BJ even
 *   though only M-1 are needed.  For weight function 8, BJ(M) is
 *   set to zero.
 *
 *   The zero-th moment of the weight function is returned in ZEMU.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   08 January 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 * Various quadrature rules can be selected by this program, indicated by
 * the parameter kind, which is an integer. Options include:
 *   - 1, Legendre,             (a,b)       1.0
 *   - 2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
 *   - 3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
 *   - 4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
 *   - 5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
 *   - 6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
 *   - 7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
 *   - 8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
 *
 * Converted for SANS by Cory Frontin in 2019
 * imported into luna by Philip Caplan in 2019
 *
 * @param kind[in]    the type of quadrature used
 * @param m[in]       the order of the Jacobi matrix
 * @param alpha[in]   floating point number for alpha, possibly required for a given quadrature
 * @param beta[in]    floating point number for beta, possibly required for a given quadrature
 * @param aj[in]      the diagonal of the Jacobi matrix
 * @param bj[in]      the subdiagonal of the Jacobi matrix
 * @return zemu       the zeroth moment
 *
 */
double
class_matrix(int kind, int m, double alpha, double beta, double aj[], double bj[])
{
  double a2b2;
  double ab;
  double aba;
  double abi;
  double abj;
  double abti;
  double apone;
  int i;
  double pi = 3.14159265358979323846264338327950;
  double temp;
  double temp2;
  double zemu = 0;

  temp = std::numeric_limits<double>::epsilon();

  parchk(kind, 2*m - 1, alpha, beta);

  temp2 = 0.5;

  if (500.0*temp < fabs(pow(tgamma(temp2), 2) - pi))
    luna_assert_msg(false,"CLASS_MATRIX: Gamma function does not match machine parameters.");

  if (kind == 1)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    ab = 0.0;

    zemu = 2.0/(ab + 1.0);

    for (i = 0; i < m; i++)
      aj[i] = 0.0;

    for (i = 1; i <= m; i++)
    {
      abi = i + ab*(i % 2);
      abj = 2*i + ab;
      bj[i - 1]= sqrt(abi*abi/(abj*abj - 1.0));
    }
  }
  else if (kind == 2)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    zemu = pi;

    for (i = 0; i < m; i++)
      aj[i] = 0.0;

    bj[0] = sqrt(0.5);
    for (i = 1; i < m; i++)
      bj[i] = 0.5;
  }
  else if (kind == 3)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    ab = alpha*2.0;
    zemu = pow(2.0, ab + 1.0)*pow(tgamma(alpha + 1.0), 2)/tgamma(ab + 2.0);

    for (i = 0; i < m; i++)
      aj[i] = 0.0;

    bj[0] = sqrt(1.0/(2.0*alpha + 3.0));
    for (i = 2; i <= m; i++)
      bj[i - 1] = sqrt(i*(i + ab)/(4.0*pow(i + alpha, 2) - 1.0));
  }
  else if (kind == 4)
  {
    ab = alpha + beta;
    abi = 2.0 + ab;
    zemu = pow(2.0, ab + 1.0)*tgamma(alpha + 1.0)*tgamma(beta + 1.0)/tgamma(abi);
    aj[0] = (beta - alpha)/abi;
    bj[0] = sqrt(4.0*(1.0 + alpha)*(1.0 + beta)/((abi + 1.0)*abi*abi));
    a2b2 = beta*beta - alpha*alpha;

    for (i = 2; i <= m; i++)
    {
      abi = 2.0*i + ab;
      aj[i - 1]= a2b2/((abi - 2.0)*abi);
      abi = abi*abi;
      bj[i - 1] = sqrt(4.0*i*(i + alpha)*(i + beta)*(i + ab)/((abi - 1.0)*abi));
    }
  }
  else if (kind == 5)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    zemu = tgamma(alpha + 1.0);

    for (i = 1; i <= m; i++)
    {
      aj[i - 1] = 2.0*i - 1.0 + alpha;
      bj[i - 1] = sqrt(i*(i + alpha));
    }
  }
  else if (kind == 6)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    zemu = tgamma((alpha + 1.0)/2.0);

    for (i = 0; i < m; i++)
      aj[i] = 0.0;

    for (i = 1; i <= m; i++)
      bj[i - 1] = sqrt((i + alpha*(i % 2))/2.0);
  }
  else if (kind == 7)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    ab = alpha;
    zemu = 2.0/(ab + 1.0);

    for (i = 0; i < m; i++)
      aj[i] = 0.0;

    for (i = 1; i <= m; i++)
    {
      abi = i + ab*(i % 2);
      abj = 2*i + ab;
      bj[i - 1] = sqrt(abi*abi/(abj*abj - 1.0));
    }
  }
  else if (kind == 8)
  {
    luna_assert_msg(false,"KIND != 4 HAS NOT BEEN UNIT TESTED.\n");

    ab = alpha + beta;
    zemu = tgamma(alpha + 1.0)*tgamma(-(ab + 1.0))/tgamma(-beta);
    apone = alpha + 1.0;
    aba = ab*apone;
    aj[0] = -apone/(ab + 2.0);
    bj[0] = -aj[0]*(beta + 1.0)/(ab + 2.0)/(ab + 3.0);
    for (i = 2; i <= m; i++)
    {
      abti = ab + 2.0*i;
      aj[i - 1] = aba + 2.0*(ab + i)*(i - 1);
      aj[i - 1] = -aj[i - 1]/abti/(abti - 2.0);
    }

    for (i = 2; i <= m - 1; i++)
    {
      abti = ab + 2.0*i;
      bj[i - 1] = i*(alpha + i)/(abti - 1.0)*(beta + i)/(abti*abti)*(ab + i)/(abti + 1.0);
    }
    bj[m - 1] = 0.0;
    for (i = 0; i < m; i++)
      bj[i] = sqrt(bj[i]);
  }
  else
    luna_assert_msg(false,"KIND NOT DETECTED.\n");


  return zemu;
}

/**
 * IMTQLX diagonalizes a symmetric tridiagonal matrix.
 *
 * Purpose:
 *   IMTQLX diagonalizes a symmetric tridiagonal matrix.
 *
 * Discussion:
 *   This routine is a slightly modified version of the EISPACK routine to
 *   perform the implicit QL algorithm on a symmetric tridiagonal matrix.
 *
 *   The authors thank the authors of EISPACK for permission to use this
 *   routine.
 *
 *   It has been modified to produce the product Q' * Z, where Z is an input
 *   vector and Q is the orthogonal matrix diagonalizing the input matrix.
 *   The changes consist (essentialy) of applying the orthogonal transformations
 *   directly to Z as they are generated.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   08 January 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 *   Roger Martin, James Wilkinson,
 *   The Implicit QL Algorithm,
 *   Numerische Mathematik,
 *   Volume 12, Number 5, December 1968, pages 377-383.
 *
 *
 * @param n       the order of the matrix
 *
 * @param d[in]   the n diagonal entries of the input matrix
 * @param d[out]  the n diagonal entries of the output matrix
 *
 * @param e[in]   the n - 1 subdiagonal entries of the input matrix; extra entry on end
 * @param e[out]  the n - 1 subdiagonal entries of the output matrix; extra entry on end
 *
 * @param z[in]   an input vector
 * @param z[out]  the value of Q'*Z where Q is the matrix that diagonalizes the input matrix
 *
 */
void
imtqlx(int n, double d[], double e[], double z[])
{
  double b;
  double c;
  double f;
  double g;
  int i;
  int ii;
  int itn= 30;
  int j;
  int k;
  int lz;
  int m= 0;
  int mml;
  double p;
  double prec;
  double r;
  double s;

  prec = std::numeric_limits<double>::epsilon();

  if (n == 1)
  {
    return;
  }

  e[n - 1] = 0.0;

  for (lz = 1; lz <= n; lz++)
  {
    j = 0;
    for ( ; ; )
    {
      for (m = lz; m <= n; m++)
      {
        if (m == n)
        {
          break;
        }

        if (fabs (e[m - 1]) <= prec*(fabs(d[m - 1]) + fabs(d[m])))
        {
          break;
        }
      }
      p = d[lz - 1];
      if (m == lz)
      {
        break;
      }
      if (itn <= j)
      {
        luna_assert_msg(false,"Iteration limit exceeded in IMTQLX.");
      }
      j = j + 1;
      g = (d[lz] - p)/(2.0*e[lz - 1]);
      r = sqrt(g*g + 1.0);
      g = d[m - 1] - p + e[lz - 1]/(g + fabs(r)*sign(g));
      s = 1.0;
      c = 1.0;
      p = 0.0;
      mml = m - lz;

      for (ii = 1; ii <= mml; ii++)
      {
        i = m - ii;
        f = s*e[i - 1];
        b = c*e[i - 1];

        if (fabs(g) <= fabs(f))
        {
          c = g/f;
          r = sqrt(c*c + 1.0);
          e[i] = f*r;
          s = 1.0/r;
          c = c*s;
        }
        else
        {
          s = f/g;
          r = sqrt(s*s + 1.0);
          e[i] = g*r;
          c = 1.0/r;
          s = s*c;
        }
        g = d[i] - p;
        r = (d[i - 1] - g)*s + 2.0*c*b;
        p = s*r;
        d[i] = g + p;
        g = c*r - b;
        f = z[i];
        z[i] = s*z[i - 1] + c*f;
        z[i - 1] = c*z[i - 1] - s*f;
      }
      d[lz - 1] = d[lz - 1] - p;
      e[lz - 1] = g;
      e[m - 1] = 0.0;
    }
  }

  //  Sorting.

  for (ii = 2; ii <= m; ii++)
  {
    i = ii - 1;
    k = i;
    p = d[i - 1];

    for (j = ii; j <= n; j++)
    {
      if (d[j - 1] < p)
      {
        k = j;
        p = d[j - 1];
      }
    }

    if (k != i)
    {
      d[k - 1] = d[i - 1];
      d[i - 1] = p;
      p = z[i - 1];
      z[i - 1] = z[k - 1];
      z[k - 1] = p;
    }
  }
  return;
}

/**
 * PARCHK checks parameters ALPHA and BETA for classical weight functions.
 *
 * Purpose:
 *   PARCHK checks parameters ALPHA and BETA for classical weight functions.
 *
 * Licensing:
 *   This code is distributed under the GNU LGPL license.
 *
 * Modified:
 *   07 January 2010
 *
 * Author:
 *   Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
 *   C++ version by John Burkardt.
 *
 * Reference:
 *   Sylvan Elhay, Jaroslav Kautsky,
 *   Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
 *   Interpolatory Quadrature,
 *   ACM Transactions on Mathematical Software,
 *   Volume 13, Number 4, December 1987, pages 399-415.
 *
 * Various quadrature rules can be selected by this program, indicated by
 * the parameter kind, which is an integer. Options include:
 *   - 1, Legendre,             (a,b)       1.0
 *   - 2, Chebyshev,            (a,b)       ((b-x)*(x-a))^(-0.5)
 *   - 3, Gegenbauer,           (a,b)       ((b-x)*(x-a))^alpha
 *   - 4, Jacobi,               (a,b)       (b-x)^alpha*(x-a)^beta
 *   - 5, Generalized Laguerre, (a,inf)     (x-a)^alpha*exp(-b*(x-a))
 *   - 6, Generalized Hermite,  (-inf,inf)  |x-a|^alpha*exp(-b*(x-a)^2)
 *   - 7, Exponential,          (a,b)       |x-(a+b)/2.0|^alpha
 *   - 8, Rational,             (a,inf)     (x-a)^alpha*(x+b)^beta
 *
 * Converted for SANS by Cory Frontin in 2019
 * imported into luna by Philip Caplan in 2019
 *
 * @param kind    the type of quadrature rule used
 * @param m       the order of the highest moment to be calculated; only for kind = 8
 * @param alpha   floating point number for alpha possibly required for a given quadrature
 * @param beta    floating point number for beta possibly required for a given quadrature
 *
 */
void
parchk(int kind, int m, double alpha, double beta)
{
  double tmp;

  if (kind <= 0)
    luna_assert_msg(false,"PARCHK - Fatal error!\nKIND <= 0.\n");

  // check ALPHA for Gegenbauer, Jacobi, Laguerre, Hermite, Exponential.
  if (3 <= kind && alpha <= -1.0)
    luna_assert_msg(false,"PARCHK - Fatal error!\n3 <= KIND and ALPHA <= -1.\n");

  // check BETA for Jacobi.
  if (kind == 4 && beta <= -1.0)
    luna_assert_msg(false,"PARCHK - Fatal error!\nKIND == 4 and BETA <= -1.0.\n");

  // check ALPHA and BETA for rational.
  if (kind == 8)
  {
    tmp = alpha + beta + m + 1.0;
    if (0.0 <= tmp || tmp <= beta)
      luna_assert_msg(false,"PARCHK - Fatal error!\nKIND == 8 but condition on ALPHA and BETA fails.\n");
  }

  return;
}

/**
 * nPointsStroudQuadrature calculates the number of points necessary for a quadrature
 * rule of the conical method of Stroud.
 *
 * Returns order^nDim.
 *
 * @param nDim[in]    number of dimensions of simplex of interest
 * @param order[in]   order of integration requested
 * @return            number of points for a quadrature rule
 */
int
nPointsStroudQuadrature(int nDim, int order)
{
  if (order==0) return 1;
  return (int) pow((double) order, (double) nDim);
}

/**
 * calculateStroudQuadrature calculates a quadrature rule with the conical method of Stroud.
 *
 * Resulting numerical integration of polynomials of up to a given order in any given
 * direction will be integrated exactly on a simplex of arbitrary dimension using the
 * resulting rule. Weight function is normalized in order to integrate to one on the
 * simplex.
 *
 * @param nDim[in]    the number of dimensions of the simplex of interest
 * @param order[in]   the order of integration requested
 * @param x[out]      the nDim*order points requested, each point of nDim coordinates
 * @param w[out]      the nDim*order weights requested
 *
 */
void
calculateStroudQuadrature(int nDim, int order, double x[], double w[])
{
  if (order==0)
  {
    w[0] = 1.0;
    luna_assert( nDim==4 );
    for (int i=0;i<nDim;i++)
      x[i] = 1./5.;
    return;
  }

  // there are $M^n$ points
  int nPoints = (int) pow((double) order, (double) nDim);
  double *y = new double[nPoints*nDim];

  //luna_assert_msg(order>0,"order = %d",order);
  luna_assert(order>0);

//  printf("nDim= %d\n", nDim);
//  printf("order= %d\n", order);
//  printf("nPoints= %d\n\n", nPoints);

  // if we want a 1D, we actually just want the one-dimensional Legendre quadrature
  if (nDim == 1)
    calculateJacobiQuadrature(order, 0.0, 0.0, x, w);
  else if (nDim == 2)
  {

    // otherwise we carry on with the conical product
    // first we need the jacobi rules for each selection possible
    double *XJacobi = new double[nDim*order];
    double *WJacobi = new double[nDim*order];

    // get the correct Jacobi-Gauss rules
    for (int j=0; j < nDim; j++)
      calculateJacobiQuadrature(order, (double) nDim - (j + 1), 0.0, XJacobi + j*order, WJacobi + j*order);

    // do the tensor product rule to get the box
    for (int j=0; j < order; j++)
      for (int i=0; i < order; i++)
      {
        int idx= order*j + i;

        y[nDim*idx + 0] = 1.0 - XJacobi[0*order + i];
        y[nDim*idx + 1] = 1.0 - XJacobi[1*order + j];

        w[idx] = WJacobi[0*order + i]*WJacobi[1*order + j];

      }

    for (int i = 0; i < nPoints; i++)
    {
      x[nDim*i + 0] = y[nDim*i + 0]*(1.0 - y[nDim*i + 1]);
      x[nDim*i + 1] = (1.0 - y[nDim*i + 0]);

      // normalize weights
      w[i] *= factorial(nDim);
    }

    // free memory
    delete [] XJacobi;
    delete [] WJacobi;

  }
  else if (nDim == 3)
  {

    // otherwise we carry on with the conical product
    // first we need the jacobi rules for each selection possible
    double *XJacobi = new double[nDim*order];
    double *WJacobi = new double[nDim*order];

    // get the correct Jacobi-Gauss rules
    for (int j = 0; j < nDim; j++)
      calculateJacobiQuadrature(order, (double) nDim - (j + 1), 0.0, XJacobi + j*order, WJacobi + j*order);

    // do the tensor product rule to get the box
      for (int k = 0; k < order; k++)
        for (int j = 0; j < order; j++)
          for (int i = 0; i < order; i++)
          {
            int idx = order*order*k + order*j + i;

            y[nDim*idx + 0] = 1.0 - XJacobi[0*order + i];
            y[nDim*idx + 1] = 1.0 - XJacobi[1*order + j];
            y[nDim*idx + 2] = 1.0 - XJacobi[2*order + k];

            w[idx] = WJacobi[0*order + i]*WJacobi[1*order + j]*WJacobi[2*order + k];

          }

    for (int i = 0; i < nPoints; i++)
    {
      x[nDim*i + 0] = y[nDim*i + 0]*y[nDim*i + 1]*(1.0 - y[nDim*i + 2]);
      x[nDim*i + 1] = (1.0 - y[nDim*i + 0]);
      x[nDim*i + 2] = y[nDim*i + 0]*(1.0 - y[nDim*i + 1]);

      // normalize weights
      w[i] *= factorial(nDim);
    }

    // free memory
    delete [] XJacobi;
    delete [] WJacobi;

  }
  else if (nDim == 4)
  {

    if (order==0)
    {
      x[0] = 1./5;
      x[1] = 1./5;
      x[2] = 1./5.;
      x[3] = 1./5;
      w[0] = 1.;
      return;
    }

    // otherwise we carry on with the conical product
    // first we need the jacobi rules for each selection possible
    double *XJacobi = new double[nDim*order];
    double *WJacobi = new double[nDim*order];

    // get the correct Jacobi-Gauss rules
    for (int j= 0; j < nDim; j++)
      calculateJacobiQuadrature(order, (double) nDim - (j + 1), 0.0, XJacobi + j*order, WJacobi + j*order);

    // do the tensor product rule to get the box
    for (int m = 0; m < order; m++)
      for (int k = 0; k < order; k++)
        for (int j = 0; j < order; j++)
          for (int i = 0; i < order; i++)
          {
            int idx= order*order*order*m + order*order*k + order*j + i;

            y[nDim*idx + 0] = 1.0 - XJacobi[0*order + i];
            y[nDim*idx + 1] = 1.0 - XJacobi[1*order + j];
            y[nDim*idx + 2] = 1.0 - XJacobi[2*order + k];
            y[nDim*idx + 3] = 1.0 - XJacobi[3*order + m];

            w[idx] = WJacobi[0*order + i]*WJacobi[1*order + j]*WJacobi[2*order + k]*WJacobi[3*order + m];

          }

    for (int i = 0; i < nPoints; i++)
    {
      x[nDim*i + 0] = y[nDim*i + 0]*y[nDim*i + 1]*y[nDim*i + 2]*(1.0 - y[nDim*i + 3]);
      x[nDim*i + 1] = (1.0 - y[nDim*i + 0]);
      x[nDim*i + 2] = y[nDim*i + 0]*(1.0 - y[nDim*i + 1]);
      x[nDim*i + 3] = y[nDim*i + 0]*y[nDim*i + 1]*(1.0 - y[nDim*i + 2]);

      // normalize weights
      w[i] *= factorial(nDim);
    }

    // free memory
    delete [] XJacobi;
    delete [] WJacobi;

  }
  else
    luna_assert_msg(false,"nDim > 4 STILL NEEDS IMPLEMENTATION.");

  // free memory
  delete [] y;

}
