#ifndef avro_LIB_NUMERICS_DETERMINANT_H_
#define avro_LIB_NUMERICS_DETERMINANT_H_

#include "common/types.h"

namespace avro
{

namespace numerics
{

template<typename type> class densMat;

template<typename type> type determinant(const densMat<type>& X);

} // numerics

} // avro

#endif
