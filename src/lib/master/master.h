#ifndef LUNA_LIB_MASTER_MASTER_H_
#define LUNA_LIB_MASTER_MASTER_H_

#include "common/types.h"

#include <memory>
#include <string>
#include <vector>

namespace luna
{

//template<typename Shape>
class Master
{
public:

  Master( coord_t number , coord_t order ) :
    number_(number),
    order_(order)
  {}

  Master( coord_t number , coord_t order , const std::string& name ) :
    number_(number),
    order_(order),
    name_(name)
  {}

  coord_t number() const { return number_; }
  coord_t order() const {return order_; }
  const std::string name() const { return name_; }

protected:
//  const Basis<Shape>& basis() const { return *basis_.get(); }
//  Basis<Shape>& basis() { return *basis_.get(); }

protected:
  coord_t number_;
  coord_t order_;
  std::string name_;

private:
//  std::unique_ptr< Basis<Shape> > basis_;
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
