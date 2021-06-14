//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_NUMERICS_QUATERNION_H_
#define AVRO_NUMERICS_QUATERNION_H_

#include "common/types.h"

#include "graphics/math.h"

#include "numerics/mat.h"

#include <math.h>

namespace avro
{

class Quaternion
{
public:
  Quaternion( const real_t theta , real_t* axis )
  {
    q_.push_back( cos(.5*theta) );
    for (coord_t d=0;d<3;d++)
      q_.push_back( axis[d]*sin(.5*theta) );
  }

  graphics::mat3
  rotation_matrix()
  {
    graphics::mat3 R;
    real_t s = magnitude();
    real_t qr = q_[0], qi = q_[1], qj = q_[2], qk = q_[3];
    R(0,0) = 1. -2.*s*(qj*qj +qk*qk);
    R(0,1) = 2.*s*(qi*qj -qk*qr);
    R(0,2) = 2.*s*(qi*qk +qj*qr);
    R(1,0) = 2.*s*(qi*qj +qk*qr);
    R(1,1) = 1. -2.*s*(qi*qi +qk*qk);
    R(1,2) = 2.*s*(qj*qk -qi*qr);
    R(2,0) = 2.*s*(qi*qk -qj*qr);
    R(2,1) = 2.*s*(qj*qk +qi*qr);
    R(2,2) = 1. -2.*s*(qi*qi +qj*qj);
    return R;
  }

  real_t
  magnitude()
  {
    // not exactly the magnitude but 's'
    real_t s = 0.;
    for (coord_t d=0;d<q_.size();d++)
      s += q_[d]*q_[d];
    return pow(s,-2.);
  }

private:
  std::vector<real_t> q_;
};

} // avro

#endif
