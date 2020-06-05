
#include "numerics/expansion.h"
#include "numerics/predicates.h"

namespace GEO {

namespace PCK {

Sign avro_side1_nd_exact_pck(const double* p0, const double* p1, const double* q0 ,unsigned short dim
){

 expansion& l = expansion_sq_dist(p0,p1,dim);
 expansion& a = expansion_dot_at(p1,q0,p0,dim).scale_fast(2.0);
 expansion& r = expansion_diff(l,a);

 Sign result = r.sign();
 if(result != ZERO) return result;

 {
   const double* p_sort[2];
   p_sort[0] = p0;
   p_sort[1] = p1;
   std::sort(p_sort,p_sort+2);
   for (index_t i=0;i<2;++i)
   {
     if(p_sort[i] == p0) return POSITIVE;
     if(p_sort[i] == p1) return NEGATIVE;
   }
 }
 return result;
}

} // PCK

} // GEO
