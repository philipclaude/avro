#ifndef URSA_LIB_LIBRARY_EGADS_H_
#define URSA_LIB_LIBRARY_EGADS_H_

#include "geometrics/body.h"
#include "geometrics/primitive.h"

struct egObject;
typedef egObject* ego;

namespace ursa
{

namespace numerics
{
  class Coordinate;
}

namespace geometrics
{

namespace EGADS
{

typedef struct
{

} egoData;

class Context
{
public:
  Context();
  Context( ego* context );
  ~Context();

  ego* get();

private:
  ego* context_;
  bool mine_;
};

class Object : public Primitive
{
public:
  Object( const Context& context , ego* object );

  void inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const;
  void evaluate( const numerics::Coordinate& u , numerics::Coordinate& p ) const;

  ego* object();
  const ego* object() const;

private:
  const Context& context_;
  ego* ego_;

  egoData data_;
};

class Body : public ursa::geometrics::Body
{
public:
  Body( const Context& context , ego* obj );

  void build();

private:
  const Context& context_;
  ego* ego_;

};

} // EGADS

} // geometrics

} // ursa

#endif
