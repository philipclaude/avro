#include "common/error.h"

#include "master/master.h"
#include "master/quadrature.h"

namespace luma
{

// needed to create a set/map of elements
bool
operator==( const Element& fx , const Element& fy )
{
  // assumes fx and fy have the same topological dimension
  // and that the indices are sorted
  luma_assert( fx.dim==fy.dim );
  for (index_t j=0;j<fx.dim;j++)
    if (fx.indices[j]!=fy.indices[j])
      return false;
  return true;
}

// needed to create a map of elements
bool
operator<( const Element& f , const Element& g )
{
  // first check the topological dimension
  if (f.dim < g.dim)
    return true;

  // lexicographically compare the indices
  return std::lexicographical_compare(f.indices.begin(), f.indices.end(),
                                      g.indices.begin(), g.indices.end());
}

template<typename Shape>
void
Master<Shape>::load_quadrature( Quadrature& quadrature )
{
  quadrature.retrieve(xquad_,wquad_);
}

template<typename Shape>
void
Master<Shape>::set_basis( BasisFunctionCategory category )
{
  basis_ = std::make_shared<Basis<Shape>>(reference_,category);
}

template class Master<Simplex>;

} // luma
