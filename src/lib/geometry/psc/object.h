#ifndef avro_LIB_GEOMETRY_PLC_H_
#define avro_LIB_GEOMETRY_PLC_H_

#include "common/types.h"

#include "geometry/body.h"
#include "geometry/entity.h"

#include "numerics/matrix.h"

namespace avro
{

namespace PSC
{

class Object : public Entity
{
public:
  Object( coord_t num , coord_t dim );

  void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const;

  numerics::MatrixD<real_t> basis_;

  void print(bool with_children) const;
  void build_hierarchy();


protected:
  coord_t dim_;
};

class Body : public avro::Body
{
public:
  Body( coord_t number , coord_t dim ) :
    avro::Body(number),
    dim_(dim)
  {}

  coord_t dim() const { return dim_; }

private:
  coord_t dim_;
};

class Facet : public Object
{
public:
  Facet( Body* body , std::vector<std::shared_ptr<Entity>>& facets ) :
    Object(facets.size()/2,body->dim())
  {
    avro_assert( number_>0 );
    for (index_t k=0;k<facets.size();k++)
      add_child(facets[k]);
  }
};

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

private:
  std::vector<real_t> x_;
};

} // PLC

} // avro

#endif
