//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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

  ~Tesseract() {}

  void build();

  void print() const
  {
    printf("tesseract!\n");
    for (index_t k=0;k<nb_entities();k++)
    {
      printf("\t");
      entity_[k]->print(true);
    }
  }

  void tessellate( BodyTessellation& body_tess ) const { avro_implement; }

private:
  std::vector<real_t> x0_;
  std::vector<real_t> length_;
};

} // library

} // avro

#endif
