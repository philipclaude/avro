#ifndef URSA_LIB_MASTER_MASTER_H_
#define URSA_LIB_MASTER_MASTER_H_

#include "common/error.h"
#include "common/types.h"

#include "master/basis.h"

#include "numerics/matrix.h"
#include "numerics/types.h"

#include <vector>

namespace ursa
{

class Quadrature;
template<typename type> class Topology;
class Vertices;

class Master
{
public:

  Master( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  coord_t number() const { return number_; }
  coord_t order() const {return order_; }

protected:

  coord_t number_;
  coord_t order_;

};

template<typename Basis> class Simplex;

template<typename Basis>
class SimplexBase : public Master
{
public:
  SimplexBase( const coord_t number , const coord_t order ) :
    Master(number,order)
  {
    precalculate();
  }

  void loadQuadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

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
  void eval() const { printf("calling bezier simplex eval\n"); }
  void eval() {}
};

class Polytope : public Master
{

};

} // ursa

#endif
