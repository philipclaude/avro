#ifndef LUNA_LIB_MASTER_MASTER_H_
#define LUNA_LIB_MASTER_MASTER_H_

#include "common/types.h"

#include "master/basis.h"
#include "master/reference.h"

#include "numerics/matrix.h"

#include <memory>
#include <string>
#include <vector>

namespace luna
{

class Quadrature;
template<typename Shape> class Basis;

template<typename Shape>
class Master
{
public:

  Master( coord_t number , coord_t order ) :
    number_(number),
    order_(order),
    reference_(number_,order_),
    basis_(nullptr)
  {}

  Master( coord_t number , coord_t order , const std::string& name ) :
    number_(number),
    order_(order),
    name_(name),
    reference_(number,order_),
    basis_(nullptr)
  {}

  void set_basis( BasisFunctionCategory category );

  coord_t number() const { return number_; }
  coord_t order() const {return order_; }
  const std::string name() const { return name_; }

  void load_quadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

  void eval( const double* x , double* phi )
  {
    luna_assert( basis_!=nullptr );
    basis_->evaluate(x,phi);
  }

protected:
  const Basis<Shape>& basis() const { return *basis_.get(); }
  Basis<Shape>& basis() { return *basis_.get(); }

protected:
  coord_t number_;
  coord_t order_;
  std::string name_;

  numerics::MatrixD<real_t> phi_;
  std::vector< numerics::MatrixD<real_t> > dphi_;

  std::vector<real_t> xquad_;
  std::vector<real_t> wquad_;

  ReferenceElement<Shape> reference_;

private:
  std::shared_ptr< Basis<Shape> > basis_;
};

typedef struct
{
  std::vector<index_t> indices;
  coord_t dim;
  bool sorted = true;
} Element;

bool operator< ( const Element& f , const Element& g );
bool operator== ( const Element& f , const Element& g );

} // luna

#endif
