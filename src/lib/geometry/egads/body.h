#ifndef LUNA_LIB_GEOMETRY_EGADS_BODY_H_
#define LUNA_LIB_GEOMETRY_EGADS_BODY_H_

#include "geometry/body.h"
#include "geometry/egads/egads_fwd.h"
#include "geometry/egads/model.h"

namespace luna
{

namespace EGADS
{

class Context;

class Body : public luna::Body
{
public:
  Body( const Context& context , ego* obj ) :
    luna::Body(0),
    context_(context),
    model_(NULL)
  {}

  Body( ego* ego , Model* model ) :
    luna::Body(0),
    context_(model->context()),
    model_(model)
  {}

  void build_hierarchy();
  void tessellate();

private:
  const Context& context_;
  Model* model_;
  ego* ego_;

};

} // EGADS

} // luna

#endif
