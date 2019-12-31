#ifndef avro_LIB_LIBRARY_TESSERACT_H_
#define avro_LIB_LIBRARY_TESSERACT_H_

#include "geometry/psc/object.h"

namespace avro
{

namespace library
{

class Tesseract : public PSC::Body
{
public:
  Tesseract( const std::vector<real_t>& x0 , const std::vector<real_t>& length ) :
    Body(4,4),
    x0_(x0),
    length_(length)
  {
    build();
  }

  void build();

  void print() const { avro_implement; }

private:
  std::vector<real_t> x0_;
  std::vector<real_t> length_;
};

} // library

} // avro

#endif
