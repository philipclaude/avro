#ifndef URSA_LIB_MASTER_MASTER_H_
#define URSA_LIB_MASTER_MASTER_H_

#include "common/types.h"

#include <string>
#include <vector>

namespace ursa
{

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
  coord_t number_;
  coord_t order_;
  std::string name_;
};

typedef struct
{
  std::vector<index_t> indices;
  coord_t dim;
  bool sorted = true;
} Element;

bool operator< ( const Element& f , const Element& g );
bool operator== ( const Element& f , const Element& g );

} // ursa

#endif
