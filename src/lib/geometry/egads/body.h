#ifndef LUNA_LIB_GEOMETRY_EGADS_BODY_H_
#define LUNA_LIB_GEOMETRY_EGADS_BODY_H_

#include "geometry/body.h"
#include "geometry/egads/data.h"
#include "geometry/egads/egads_fwd.h"
#include "geometry/egads/model.h"

#include <map>

namespace luna
{

namespace EGADS
{

class Context;

class Body : public luna::Body
{
public:
  Body( const Context& context , ego* object ) :
    luna::Body(0),
    context_(context),
    model_(NULL),
    object_(object)
  {}

  Body( ego* object , Model* model ) :
    luna::Body(0),
    context_(model->context()),
    model_(model),
    object_(object)
  {}

  const Context& context() const { return context_; }

  void build_hierarchy();
  void tessellate();

  Entity_ptr child( index_t k );
  void add_child( ego object , Entity_ptr entity );

  Entity_ptr lookup( ego object ) const;

  void print() const;

  ego* object() { return object_; }
  const ego* object() const { return object_; }

private:
  const Context& context_;
  EGADS::Model* model_;
  ego* object_;
  egoData data_;

  std::map<ego,Entity_ptr> children_;
};

} // EGADS

} // luna

#endif
