#ifndef AVRO_LIB_PHYSICS_FLUID_H_
#define AVRO_LIB_PHYSICS_FLUID_H_

#include "common/parameters.h"

#include "voronoi/laguerre.h"

namespace avro
{

template<typename PDE_t>
class OptimalTransportSimulator
{

  typedef PDE_t::nb_equations nb_pde;

public:
  OptimalTransportSimulator();

  void

  void step();

private:
  PDE_t& pde_;
  LaguerreDiagram diagram_;

  std::vector<real_t> density_;
  std::vector<real_t> mass_;
  Points position_;
  DOF<real_t> velocity_;

};

class PDE
{

};

class PDE_Parameters : public Parameters
{
public:
  PDE_Parameters();

  void standard();

  real_t& reference_density() { return realParams_["reference_density"]; }
  real_t& reference_velocity() { return realParams_["reference_velocity"]; }

};

template<int dim>
class IncompressibleNavierStokes : public PDE
{
public:
  typedef dim+1 nb_equations; // mass first, then dim momentum equations

  IncompressibleNavierStokes( const PDE_Parameters& params );

  // computes the force for particle k
  real_t compute_force( index_t k , const LaguerreDiagram& laguerre ) const;

private:
  const PDE_Parameters& params_;

};

} // avro

#endif
