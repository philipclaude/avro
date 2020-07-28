#ifndef AVRO_LIB_SIMULATION_EQUATIONS_H_
#define AVRO_LIB_SIMULATION_EQUATIONS_H_

#include "numerics/sparse_matrix.h"

class Residual;
class Jacobian;

namespace avro
{

template<typename Residual_type>
class EquationSet
{

public:
  EquationSet( Residual_type& res );

private:
  SparseMatrix<real_t> matrix_;
  Residual_type& residual_;

};

} // avro

#endif
