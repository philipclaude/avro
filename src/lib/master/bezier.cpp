#include "master/basis.h"
#include "master/simplex.h"

namespace luma
{

void
Bezier<Simplex>::eval(const ReferenceElement<Simplex>& ref , const double* x , double* phi )
{
  printf("eval in bezier!\n");
}

void
Bezier<Simplex>::grad( const ReferenceElement<Simplex>& ref , const double* x , double* gphi )
{
  printf("grad in bezier!\n");
}

void
Bezier<Simplex>::hess(const ReferenceElement<Simplex>& ref , const double* x , double* hphi )
{
  printf("hess in bezier!\n");
}

} // luma
