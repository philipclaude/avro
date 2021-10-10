//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_EGADS_CONTEXT_H_
#define avro_LIB_GEOMETRY_EGADS_CONTEXT_H_

#include "geometry/egads/egads_fwd.h"

namespace avro
{

namespace EGADS
{

class Context
{
public:
  Context();
  Context( ego context );
  ~Context();

  ego get();
  const ego get() const;

  void print() const;

private:
  ego  context_;
  bool mine_;
};

} // EGADS

} // avro

#endif
