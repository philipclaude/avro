#ifndef LUNA_LIB_GEOMETRY_EGADS_H_
#define LUNA_LIB_GEOMETRY_EGADS_H_

#include "geometry/body.h"
#include "geometry/entity.h"

struct egObject;
typedef egObject* ego;

namespace luna
{

namespace numerics
{
  class Coordinate;
}

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

class Object : public Entity
{
public:
  Object( const Context& context , ego* object );

  void inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const;
  void evaluate( const numerics::Coordinate& u , numerics::Coordinate& p ) const;

  ego* object();
  const ego* object() const;

  void print() const { luna_implement; }

private:
  const Context& context_;
  ego* ego_;

  egoData data_;
};

class Body : public luna::Body
{
public:
  Body( const Context& context , ego* obj );

  void build();

private:
  const Context& context_;
  ego* ego_;

};

} // EGADS

} // luna

#endif
