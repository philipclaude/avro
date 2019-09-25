#ifndef URSA_LIB_MASTER_MASTER_H_
#define URSA_LIB_MASTER_MASTER_H_

#include "common/error.h"
#include "common/types.h"

#include "numerics/matrix.h"
#include "numerics/types.h"

#include <vector>

namespace ursa
{

class Quadrature;
template<typename type> class Topology;
class Vertices;

template<typename ShapeBasis>
class Master
{
public:
  void eval( index_t elem , const ParaCoord& x , real& f ) const
  {
    ursa_assert( topology_!=nullptr );
    std::vector<index_t> indices;
    //topology_->retrieve( elem , indices );
    //derived().eval( topology_->vertices() , indices.data() , indices.size() , x , &f , NULL , NULL );
  }

  void setTopology( Topology< Master<ShapeBasis> >* topology_ );

  Master( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  coord_t order() const { return order_; }

protected:

  Topology< Master<ShapeBasis> >* topology_;

  coord_t number_;
  coord_t order_;

private:
  ShapeBasis& derived() { return *static_cast<ShapeBasis*>(this); }
  const ShapeBasis& derived() const { return *static_cast<const ShapeBasis*>(this); }

};

class Simplex
{
public:
  void loadQuadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

protected:

  void eval( const Vertices& p , const index_t* v , index_t nv , ParaCoord& x , real* f , Gradient<double>* g , Gradient<double>* H );

  Simplex( const coord_t& number , const coord_t& order ) :
    number_(number), order_(order)
  {
    precalculate();
  }

  index_t nb_poly() const { return phi_.m(); }
  index_t nb_quad() const { return phi_.n(); }

private:

  void precalculate();

  numerics::MatrixD<real> phi_;
  std::vector< numerics::MatrixD<real> > dphi_;

  std::vector<real> xquad_;
  std::vector<real> wquad_;

  const coord_t& number_; // references inherited from Master<Basis>
  const coord_t& order_;
};

class LagrangeSimplex : public Simplex, public Master<LagrangeSimplex>
{
public:
  using Master<LagrangeSimplex>::number_;
  using Master<LagrangeSimplex>::order_;

  LagrangeSimplex( coord_t number , coord_t order ) :
    Master<LagrangeSimplex>(number,order),
    Simplex(number_,order_)
  {}

  template<typename T> void eval( index_t elem , const ParaCoord& x , T& f ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Gradient<T>& g ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Hessian<T>& h ) const;
};

/*
class BezierSimplex : public Simplex<BezierSimplex>
{
public:
  template<typename T> void eval( index_t elem , const ParaCoord& x , T& f ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Gradient<T>& g ) const;
  template<typename T> void eval( index_t elem , const ParaCoord& x , Hessian<T>& h ) const;
};
*/

class Polytope : public Master<Polytope>
{

};

} // ursa

#endif
