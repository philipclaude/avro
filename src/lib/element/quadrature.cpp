//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "element/quadrature.h"

#include <stdio.h>

#include "element/quadrature_cp.h"
#include "element/simplex.h"

#include "numerics/linear_algebra.h"

namespace avro
{

QuadratureStore<Simplex> __store_simplex_lagrange__(BasisFunctionCategory_Lagrange);
QuadratureStore<Simplex> __store_simplex_legendre__(BasisFunctionCategory_Legendre);
QuadratureStore<Simplex> __store_simplex_bernstein__(BasisFunctionCategory_Bernstein);

template<typename type>
QuadratureStore<type>::QuadratureStore( BasisFunctionCategory category ) :
  category_(category) {

  nmax = 4;
  pmax = 4;
  qmax = 7;

  build();
}

template<typename type>
void
QuadratureStore<type>::build() {

  // create all the quadratures
  for (index_t n = 1; n <= nmax; n++) {
    quadrature_.insert( { n , std::map<coord_t,std::shared_ptr<Quadrature>>() } );
    for (index_t q = 1; q <= qmax; q++) {
      quadrature_[n].insert( {q,std::make_shared<ConicalProductQuadrature>(n,q)} );
      quadrature_[n][q]->define();
    }
  }

  memory_ = 0;
  for (index_t n = 1; n <= nmax; n++) {

    if (n == 4 && category_ == BasisFunctionCategory_Bernstein) {
      pmax = 1;
    }

    basis_.insert( { n , std::map< coord_t,std::map<coord_t,std::vector<matd<real_t>>>>() } );
    for (index_t p = 0; p <= pmax; p++)
    {
      basis_[n].insert( { p , std::map<coord_t,std::vector<matd<real_t>>>() } );
      const index_t nb_basis = type::nb_basis(n,p);
      std::vector<real_t> phi( nb_basis );
      std::vector<real_t> gphi( nb_basis*n );

      // initialize the basis function evaluator
      ReferenceElement<Simplex> reference(n,p);
      Basis<Simplex> basis(n,p,category_);

      // evaluate the basis at the Lagrange nodes
      matd<real_t> A( nb_basis , nb_basis );
      for (index_t i = 0; i < nb_basis; i++)
      {
        basis.evaluate( reference.get_reference_coordinate(i) , phi.data() );
        for (index_t j = 0; j < nb_basis; j++)
          A(j,i) = phi[j];
      }

      matd<real_t> Ainv(nb_basis,nb_basis);
      numerics::inverseLUP(A,Ainv);

      for (index_t q = 1; q <= qmax; q++)
      {
        // retrieve the quadrature rule
        const Quadrature& quad = quadrature(n,q);
        std::vector<matd<real_t>> B(n+1);
        for (index_t d = 0; d <= n; d++)
          B[d] = matd<real_t>( nb_basis , quad.nb() );

        // loop through the quadrature points and evaluate the basis
        for (index_t j = 0; j < quad.nb(); j++)
        {
          // evaluate the basis function at this quadrature points
          basis.evaluate( quad.x(j) , phi.data() );

          for (index_t i = 0; i < nb_basis; i++)
            B[0](i,j) = phi[i];

          basis.gradient( quad.x(j) , gphi.data() );
          for (coord_t d = 0; d < n; d++)
          {
            for (index_t i = 0; i < nb_basis; i++)
              B[d+1](i,j) = gphi[ d*nb_basis + i ];
          }
        }

        for (index_t d = 0; d < B.size(); d++)
        {
          B[d] = Ainv*B[d]; // convert modal basis to nodal basis
          memory_ += B[d].memory();
        }

        basis_[n][p].insert( { q , B } );
      }
    }
  }
}

Quadrature::Quadrature( coord_t dim , const int order) :
  dim_(dim),
  order_(order),
  nb_quad_(0),
  defined_(false)
{}

const real_t*
Quadrature::x( index_t k ) const
{
  avro_assert( k < nb_quad_ );
  return x_.data()+k*dim_;
}

real_t
Quadrature::w( index_t k ) const
{
  avro_assert( k < nb_quad_ );
  return w_[k];
}

void
ConicalProductQuadrature::define()
{
  if (order_<0) order_ = order_max_;

  nb_quad_ = nPointsStroudQuadrature(dim_,order_);
  x_.resize( nb_quad_*dim_ );
  w_.resize( nb_quad_ );

  calculateStroudQuadrature( dim_ , order_ , x_.data() , w_.data() );

  // volume of an orthogonal m-simplex
  volume_ = 1.0;
  for ( int i = 1; i <= dim_; i++ )
    volume_ = volume_ / ( ( real_t ) i );

  defined_ = true;
}

int
i4_choose( int n, int k )
{
  int i;
  int mn;
  int mx;
  int value;

  mn = std::min( k , n -k );

  if ( mn < 0 )
  {
    value = 0;
  }
  else if ( mn == 0 )
  {
    value = 1;
  }
  else
  {
    mx = std::max( k , n -k );
    value = mx + 1;

    for ( i = 2; i <= mn; i++ )
    {
      value = ( value * ( mx + i ) ) / i;
    }
  }
  return value;
}

void
comp_next( int n, int k, int a[], bool &more, int &h, int &t )
{
  int i;

  if ( !( more ) )
  {
    t = n;
    h = 0;
    a[0] = n;
    for ( i = 1; i < k; i++ )
    {
       a[i] = 0;
    }
  }
  else
  {
    if ( 1 < t )
    {
      h = 0;
    }
    h = h + 1;
    t = a[h-1];
    a[h-1] = 0;
    a[0] = t - 1;
    a[h] = a[h] + 1;
  }

  more = ( a[k-1] != n );

  return;
}

int
gm_rule_size( int rule, int m )
{
  int arg1;
  int n;

  arg1 = m + rule + 1;
  n = i4_choose( arg1, rule );
  return n;
}

void
GrundmannMoellerQuadrature::define()
{
  int *beta;
  int beta_sum;
  int d;
  int dim;
  int h;
  int i;
  int j;
  int j_hi;
  int k;
  bool more;
  double one_pm;
  int s;
  int t;
  double weight;

  if (order_<0) order_ = order_max_;
  if (order_%2==0)
    rule_ = order_/2;
  else
    rule_ = (order_-1)/2;

  nb_quad_ = gm_rule_size( rule_ , dim_ );

  int m = dim_;

  x_.resize( nb_quad_*dim_ );
  w_.resize( nb_quad_ );

  s = rule_;
  d = 2 * s + 1;
  k = 0;
  one_pm = 1.0;

  beta = new int[m+1];

  for ( i = 0; i <= s; i++ )
  {
    weight = ( double ) one_pm;

    j_hi = std::max( m , std::max( d , d +m -i ) );

    for ( j = 1; j <= j_hi; j++ )
    {
      if ( j <= m )
      {
        weight = weight * ( double ) ( j );
      }
      if ( j <= d )
      {
        weight = weight * ( double ) ( d + m - 2 * i );
      }
      if ( j <= 2 * s )
      {
        weight = weight / 2.0;
      }
      if ( j <= i )
      {
        weight = weight / ( double ) ( j );
      }
      if ( j <= d + m - i )
      {
        weight = weight / ( double ) ( j );
      }
    }

    one_pm = - one_pm;

    beta_sum = s - i;
    more = false;
    h = 0;
    t = 0;

    for (;;)
    {
      comp_next ( beta_sum, m + 1, beta, more, h, t );

      w_[k] = weight;
      for ( dim = 0; dim < m; dim++ )
      {
        x_[dim+k*m] = ( double ) ( 2 * beta[dim+1] + 1 )
                    / ( double ) ( d + m - 2 * i );
      }
      k = k + 1;

      if ( !more )
      {
        break;
      }
    }
  }

  // volume of an orthogonal m-simplex
  volume_ = 1.0;
  for ( int i = 1; i <= m; i++ )
    volume_ = volume_ / ( ( real_t ) i );

  // free memory.
  delete [] beta;

  defined_ = true;

  return;
}

} // avro
