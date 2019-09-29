#ifndef URSA_LIB_GEOMETRICS_PRIMITIVE_H_
#define URSA_LIB_GEOMETRICS_PRIMITIVE_H_

#include "common/tree.h"
#include "common/types.h"

#include <string>

namespace ursa
{

namespace numerics
{
  class Coordinate;
}

namespace geometrics
{

class Body;

class Primitive : public Tree<Primitive>
{
public:
  coord_t number() const { return number_; }
  std::string name() const { return name_; }

  index_t identifier() const { return identifier_; }
  Body* body() const { return body_; }

protected:
  Primitive( coord_t number );
  Primitive( coord_t number , const std::string& name );

  virtual ~Primitive() {}

  virtual void inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const = 0;
  virtual void evaluate( const numerics::Coordinate& u , numerics::Coordinate& p ) const = 0;

  Primitive* intersect( Primitive* e1 );
  Primitive* intersect( Primitive* e1 , Primitive* e2 , bool only_check=false );
  Primitive* intersect( Primitive* e1 , Primitive* e2 , Primitive* e3 );

  bool interior() const { return interior_; }
  void interior( bool x ) { interior_ = x; }

  bool sense() const { return sense_; }
  void sense( bool x ) { sense_ = x; }

  coord_t number_;
  std::string name_;
  index_t identifier_;
  Body* body_;

  bool interior_;
  bool sense_;
};

} // geometrics

} // ursa

#endif
