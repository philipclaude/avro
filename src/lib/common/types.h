#ifndef avro_COMMON_TYPES_H_
#define avro_COMMON_TYPES_H_

#ifndef nil
#include <cstdlib>
#define nil NULL
#endif

namespace avro
{

typedef unsigned short coord_t;
typedef unsigned long  index_t;
typedef double         real_t;
typedef coord_t        coord_index_t;

enum Sign
{
  NEGATIVE = -1,
  ZERO = 0,
  POSITIVE = 1
};

} // avro

#endif
