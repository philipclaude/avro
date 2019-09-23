#ifndef URSA_LIB_MASTER_MASTER_H_
#define URSA_LIB_MASTER_MASTER_H_

#include "common/types.h"

#include "numerics/matrix.h"
#include "numerics/types.h"

namespace ursa
{

class Quadrature;
template<typename type> class Topology;

template<typename Basis>
class Master
{
public:
  template<typename T> void eval( index_t elem , const ParaCoord& x , T& f ) const
    { derived().eval(elem,x,f); }
  template<typename T> void eval( index_t elem , const ParaCoord& x , Gradient<T>& g ) const
    { derived().eval(elem,x,g); }
  template<typename T> void eval( index_t elem , const ParaCoord& x , Hessian<T>& h ) const
    { derived().eval(elem,x,h); }

  template<typename T> void eval( index_t elem , const PhysCoord& x , T& f ) const
    { derived().eval(elem,x,f); }
  template<typename T> void eval( index_t elem , const PhysCoord& x , Gradient<T>& g ) const
    { derived().eval(elem,x,g); }
  template<typename T> void eval( index_t elem , const PhysCoord& x , Hessian<T>& h ) const
    { derived().eval(elem,x,h); }

private:
  Basis& derived() { return static_cast<Basis>(*this); }
  const Basis& derived() const { return static_cast<const Basis>(*this); }

  Topology<Master<Basis>>* topology_;
};

template<typename Basis>
class Simplex : public Master<Simplex<Basis>>
{
public:
  template<typename T> void eval( index_t elem , const ParaCoord& x , T& f ) const
    { derived().eval(elem,x,f); }
  template<typename T> void eval( index_t elem , const ParaCoord& x , Gradient<T>& g ) const
    { derived().eval(elem,x,g); }
  template<typename T> void eval( index_t elem , const ParaCoord& x , Hessian<T>& h ) const
    { derived().eval(elem,x,h); }

  Simplex( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {
    precalculate();
  }

  index_t nb_poly() const { return phi_.m(); }
  index_t nb_quad() const { return phi_.n(); }

  void loadQuadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

  coord_t number_;
  index_t order_;

private:
  Basis& derived() { return *static_cast<Basis*>(this); }
  const Basis& derived() const { return *static_cast<const Basis*>(this); }

  void precalculate();

  numerics::MatrixD<real> phi_;
  std::vector< numerics::MatrixD<real> > dphi_;

  std::vector<real> xquad_;
  std::vector<real> wquad_;
};

class LagrangeSimplex : public Simplex<LagrangeSimplex>
{
public:
  using Simplex<LagrangeSimplex>::Simplex;

  template<typename T> void eval( index_t elem , const ParaCoord& x , T& f ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Gradient<T>& g ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Hessian<T>& h ) const;
};

class BezierSimplex : public Simplex<BezierSimplex>
{
public:
  template<typename T> void eval( index_t elem , const ParaCoord& x , T& f ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Gradient<T>& g ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Hessian<T>& h ) const;
};

class Polytope : public Master<Polytope>
{

};

} // ursa

#endif
