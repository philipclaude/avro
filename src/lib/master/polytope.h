#ifndef URSA_LIB_MASTER_POLYTOPE_H_
#define URSA_LIB_MASTER_POLYTOPE_H_

#include "common/data.h"

#include "master/master.h"
#include "master/simplex.h"

namespace ursa
{

class Polytope : public Master
{

public:
  Polytope( coord_t number , coord_t order , Data<int>& incidence );

  Data<int>& incidence() { return incidence_; }
  const Data<int>& incidence() const { return incidence_; }

private:
  Simplex<Lagrange> simplex_;
  Data<int>& incidence_;
};

} // ursa

#endif
