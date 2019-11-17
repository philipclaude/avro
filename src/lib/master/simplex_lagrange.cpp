#include "common/tools.h"

#include "master/simplex.h"

namespace luna
{

Simplex<Lagrange>::Simplex( coord_t number , coord_t order ) :
  SimplexBase(number,order)
{
  precalculate();

  for (index_t k=0;k<number_+1;k++)
  for (index_t i=k+1;i<number_+1;i++)
  {
    this->edges_.push_back(k);
    this->edges_.push_back(i);
  }
}

static void
next_index( const int n , const int q , bool& more , std::vector<index_t>& x )
{
  int i,j;

  luna_assert( x.size()==index_t(n) );

  if (!more)
  {
    if (n<1)
      luna_assert(false);


    more = true;
    j = 1;

    x[0] = q;
    for (i=1;i<n;i++)
      x[i] = 0;

    if (n==1)
      more = false;
  }
  else
  {
    j = n -1;
    for (i=n-2;0<=i;i--)
    {
      if (0 < x[i])
      {
        j = i;
        break;
      }
    }

    x[j] = x[j] -1;
    x[j+1] = q;

    for (i=0;i<=j;i++)
      x[j+1] = x[j+1] -x[i];

    for (i=j+2;i<n;i++)
      x[i] = 0;

    if (x[n-1]==index_t(q))
      more = false;
  }
}

void
Simplex<Lagrange>::precalculate()
{
  // set the unit (equilateral) coordinates
  if (number_==0)
  {
    xunit_.push_back(1.); // not sure?
  }
  else if (number_==1)
  {
    xunit_.push_back(1.);
    xunit_.push_back(0.);
  }
  else if (number_==2)
  {
    xunit_.push_back(1.); xunit_.push_back(0.);
    xunit_.push_back(-.5); xunit_.push_back(  std::sqrt(.25*3.) );
    xunit_.push_back(-.5); xunit_.push_back( -std::sqrt(.25*3.) );
  }
  else if (number_==3)
  {
    xunit_.push_back(1.); xunit_.push_back(0.); xunit_.push_back(0.);
    xunit_.push_back(-1./3); xunit_.push_back(  std::sqrt(8./9.) ); xunit_.push_back(  0.);
    xunit_.push_back(-1./3); xunit_.push_back( -std::sqrt(2./9.) ); xunit_.push_back(  std::sqrt(2./3.));
    xunit_.push_back(-1./3); xunit_.push_back( -std::sqrt(2./9.) ); xunit_.push_back( -std::sqrt(2./3.));
  }
  else if (number_==4)
  {
    // this is not actually unit...lengths are sqrt(5/2)
    xunit_.push_back(1.); xunit_.push_back(0.); xunit_.push_back(0.); xunit_.push_back(0.);
    xunit_.push_back(-.25); xunit_.push_back(  std::sqrt(15./16.) ); xunit_.push_back(0.); xunit_.push_back(0.);
    xunit_.push_back(-.25); xunit_.push_back( -std::sqrt( 5./48.) ); xunit_.push_back(  std::sqrt(5./6. ) ); xunit_.push_back(0.);
    xunit_.push_back(-.25); xunit_.push_back( -std::sqrt( 5./48.) ); xunit_.push_back( -std::sqrt(5./24.) ); xunit_.push_back(  std::sqrt(5./8.) );
    xunit_.push_back(-.25); xunit_.push_back( -std::sqrt( 5./48.) ); xunit_.push_back( -std::sqrt(5./24.) ); xunit_.push_back( -std::sqrt(5./8.) );

    for (index_t k=0;k<xunit_.size();k++)
      xunit_[k] *= sqrt(2./5.);
  }
  else
  {
    printf("warning: this is untested\n");
    numerics::MatrixD<real_t> X(number_,number_+1);
    X = 0;

    for (coord_t i=0;i<number_;i++)
      X(i,i) = 1.;
    real_t a = (1. -std::sqrt(1.+number_))/number_;
    for (coord_t i=0;i<number_;i++)
      X(i,number_) = a;

    std::vector<real_t> c(number_,0.);
    for (coord_t i=0;i<number_;i++)
    for (coord_t j=0;j<number_+1;j++)
      c[i] += X(i,j)/(number_+1);

    for (coord_t j=0;j<number_+1;j++)
    for (coord_t i=0;i<number_;i++)
      X(i,j) = X(i,j) -c[i];

    real_t s = 0.;
    for (coord_t i=0;i<number_;i++)
      s += X(i,0)*X(i,0);
    s = std::sqrt(s);
    for (coord_t i=0;i<number_;i++)
    for (coord_t j=0;j<number_+1;j++)
      X(i,j) = X(i,j)/s;

    for (coord_t i=0;i<number_+1;i++)
    for (coord_t j=0;j<number_;j++)
      xunit_.push_back( X(j,i) );
  }

  // set the orthogonal coordinates
  xorth_.resize( number_*(number_+1) );
  std::fill( xorth_.begin() , xorth_.end() , 0. );
  for (index_t k=0;k<number_;k++)
    xorth_[ (k+1)*number_ +k ] = 1.;

  real_t length = order_;
  if (order_==0)
  {
    // it doesn't really make sense to have a zero-order lagrange element
    // but we can hack it by setting the coordinates to the centroid

    for (index_t j=0;j<number_;j++)
    {
      lref_.push_back( 1 );
      xref_.push_back( 1./(number_+1) );
    }
    return;
  }

  luna_assert( length>0 );

  bool more = false;
  std::vector<index_t> xp(number_+1);
  for (index_t i=0;;i++)
  {
    next_index( number_+1 , order_ , more , xp );

    for (index_t j=0;j<xp.size();j++)
    {
      lref_.push_back(xp[j]);
      xref_.push_back( real_t(xp[j])/real_t(length));
    }

    if (!more) break;
  }
}

void
Simplex<Lagrange>::get_edges( const index_t* v , const index_t nv , std::vector<index_t>& ek ) const
{
  ek.resize( edges_.size() );
  for (index_t j=0;j<edges_.size();j++)
    ek[j] = v[edges_[j]];
}

const real_t*
Simplex<Lagrange>::get_reference_coordinate( index_t k ) const
{
  return &xref_[k*(number_+1)];
}

const index_t*
Simplex<Lagrange>::get_lattice_coordinate( index_t k ) const
{
  return &lref_[k*(number_+1)];
}

template<typename type>
static type
vector_sum( index_t n, const type* a )
{
  type value(0);

  for (index_t i = 0; i < n; i++ )
  {
    value = value + a[i];
  }
  return value;
}

template<typename type>
static type
eval_lagrange_basis( index_t m, const index_t *alpha , const type* x )
{
  type c,l,w;

  index_t d = vector_sum<index_t>( m+1 , alpha );

  l = 1.0;
  c = 1.0;

  for (index_t q = 0; q < m; q++ )
  {
    for (index_t p = 0; p < alpha[q]; p++ )
    {
      l = l * ( type(d) * x[q] - type(p) );
      c = c * ( alpha[q]       - type(p) );
    }
  }

  w = 1.0 - vector_sum ( m, x );

  for (index_t p = 0; p < alpha[m]; p++ )
  {
    l = l * ( type(d) * w    - type(p) );
    c = c * ( type(alpha[m]) - type(p) );
  }

  l = l / c;

  return l;
}

void
Simplex<Lagrange>::evaluate( const real_t* x , std::vector<real_t>& phi ) const
{
  luna_assert( phi.size()==nb_basis() );
  for (index_t k=0;k<nb_basis();k++)
  {
    // get the lattice coordinates for this basis function
    const index_t* alpha = get_lattice_coordinate(k);

    // evaluate the basis function
    phi[k] = eval_lagrange_basis( number_ , alpha , x );
  }
}

void
Simplex<Lagrange>::evaluate( index_t k , std::vector<real_t>& phi ) const
{
  // evaluate the set of basis functions at the k'th interpolation point
  evaluate( get_reference_coordinate(k) , phi );
}

} // luna
