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

#include "geometry/psc/facet.h"
#include "geometry/psc/node.h"
#include "geometry/psc/object.h"
#include "mesh/points.h"

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


  static void expansion( const real_t* x , real_t* y ) {
    std::vector<real_t> X(x,x+4);

    // translate to the center of the domain (0.5)^4
    //for (coord_t d = 0; d < 4; d++)
    //  X[d] -= 0.5;

    // scale in time
    real_t t = x[3]; // time between 0 and 1
    real_t s = 1 + t;
    for (coord_t d = 0; d < 3; d++)
      X[d] *= s;
    X[3] *= 2;

    // translate back
    for (coord_t d = 0; d < 4; d++)
      y[d] = X[d];// + 0.5;
  }
  static void closingwall( const real_t* x , real_t* y ) {
    std::vector<real_t> X(x,x+4);

    // translate to the center of the domain (0.5)^4
    //for (coord_t d = 0; d < 4; d++)
    //  X[d] -= 0.5;

    // scale in time
    real_t t = x[3]; // time between 0 and 1
    real_t s = 1 - 0.5*t;
    X[2] *= s;
    X[3] *= 1;

    // translate back
    for (coord_t d = 0; d < 4; d++)
      y[d] = X[d];// + 0.5;
  }

  template<typename F>
  void
  map_to( const F& function , Points* points=nullptr ) {

    // temporary coordinates to do the evaluation
    std::vector<real_t> x(4);

    std::vector<Entity*> entities;
    get_entities(entities);

    // map the nodes
    for (index_t k = 0; k < entities.size(); k++) {

      if (entities[k]->number() != 0) continue;
      PSC::Node& node = *static_cast<PSC::Node*>(entities[k]);
      function( node.x() , x.data() );

      for (coord_t d = 0; d < x.size(); d++)
        node(d) = x[d];
    }

    // we need to rebuild the basis for all the entities
    for (index_t k = 0; k < entities.size(); k++) {

      if (entities[k]->number() == 0) continue;
      PSC::Facet& facet = *static_cast<PSC::Facet*>(entities[k]);
      facet.build_basis();
    }

    // map the provided points
    if (points == nullptr) return;
    for (index_t k = 0; k < points->nb(); k++) {
      function( (*points)[k] , x.data() );
      for (coord_t d = 0; d < x.size(); d++)
        (*points)[k][d] = x[d];
    }
  }

private:
  std::vector<real_t> x0_;
  std::vector<real_t> length_;
};

} // library

} // avro

#endif
