#ifndef avro_LIB_GEOMETRY_PSC_NODE_H_
#define avro_LIB_GEOMETRY_PSC_NODE_H_

#include "geometry/psc/object.h"

namespace avro
{

namespace PSC
{

class Node : public Object
{
public:
  Node( Body* body , real_t* data ) :
    Object(0,body->dim()),
    x_(data,data+dim_)
  {}

  real_t operator()(coord_t d) const
  {
    avro_assert( d<dim_ );
    return x_[d];
  }

  void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const
  { inverse(x,u); }
  void evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const;


private:
  std::vector<real_t> x_;
};

} // PSC

} // avro

#endif
