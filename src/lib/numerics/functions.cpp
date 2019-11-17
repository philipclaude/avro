#include "numerics/functions.h"

namespace luna
{

namespace numerics
{


index_t
factorial( const index_t i )
{
  if (i==0 || i==1) return 1;
  return i*factorial(i-1);
}

index_t
binomial( const index_t n , const index_t k )
{
	index_t num=1,den=1;
	for (index_t i=1;i<=k;i++)
	{
		num *= (n +1 -i);
		den *= i;
	}
	return num/den;;
}

} // numerics

} // avro
