#ifndef URSA_LIB_MASTER_SIMPLEX_H_
#define URSA_LIB_MASTER_SIMPLEX_H_

#include "master/basis.h"
#include "master/master.h"

#include "numerics/matrix.h"
#include "numerics/types.h"

#include <vector>

namespace ursa
{

template<typename Shape_t> class Topology;

template<typename Basis> class Simplex;
class Quadrature;

template<typename Basis>
class SimplexBase : public Master
{
public:
  SimplexBase( const coord_t number , const coord_t order ) :
    Master(number,order,"simplex")
  {
    precalculate();
  }

  SimplexBase( const Topology<Simplex<Basis>>& topology , const coord_t order );

  void loadQuadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

  template<typename BasisFrom_t,typename dof_t> void convert( const BasisFrom_t& masterFrom , const std::vector<dof_t>& A , std::vector<dof_t>& B ) const;

  index_t nb_poly() const { return phi_.m(); }
  index_t nb_quad() const { return phi_.n(); }

private:

  void precalculate();

  numerics::MatrixD<real_t> phi_;
  std::vector< numerics::MatrixD<real_t> > dphi_;

  std::vector<real_t> xquad_;
  std::vector<real_t> wquad_;
};

template<>
class Simplex<Lagrange> : public SimplexBase<Lagrange>
{
public:
  using SimplexBase<Lagrange>::SimplexBase;

  void eval() const { printf("calling lagrange simplex eval\n"); }
  void eval() {}
};

template<>
class Simplex<Bezier> : public SimplexBase<Bezier>
{
public:
  using SimplexBase<Bezier>::SimplexBase;

  void eval() const { printf("calling bezier simplex eval\n"); }
  void eval() {}
};

} // ursa

#endif
