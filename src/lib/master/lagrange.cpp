#include "common/tools.h"

#include "master/basis.h"
#include "master/simplex.h"

namespace avro
{

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
Lagrange<Simplex>::eval( const ReferenceElement<Simplex>& reference , const real_t* x , real_t* phi )
{
  for (index_t k=0;k<reference.nb_basis();k++)
  {
    // get the lattice coordinates for this basis function
    const index_t* alpha = reference.get_lattice_coordinate(k);

    // evaluate the basis function
    phi[k] = eval_lagrange_basis( reference.number() , alpha , x );
  }
}

void
Lagrange<Simplex>::grad( const ReferenceElement<Simplex>& ref , const double* x , double* gphi )
{
  printf("grad in lagrange!\n");
}

void
Lagrange<Simplex>::hess( const ReferenceElement<Simplex>& ref , const double* x , double* hphi )
{
  printf("hess in lagrange!\n");
}

#if 0

void
Simplex::evaluate( index_t k , std::vector<real_t>& phi ) const
{
  // evaluate the set of basis functions at the k'th interpolation point
  evaluate( get_reference_coordinate(k) , phi );
}
#endif

} // avro
