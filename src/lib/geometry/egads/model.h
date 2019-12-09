#ifndef LUNA_LIB_GEOMETRY_EGADS_MODEL_H_
#define LUNA_LIB_GEOMETRY_EGADS_MODEL_H_

#include "geometry/model.h"

#include "geometry/egads/egads_fwd.h"

#include <string>

namespace luna
{

namespace EGADS
{

class Body;
class Context;

class Model : public luna::Model
{
public:
  Model( const Context& context );
  Model( const Context& context , const std::string& filename , bool split=false );

  Body& body(index_t k) { return *body_[k].get(); }
  const Body& body(index_t k) const { return *body_[k].get(); }

  index_t nb_bodies() const { return body_.size(); }

  void determine_number();

  const Context& context() const { return context_; }

  Entity* find_entity( index_t id , int object_class ) const;

  void print() const;

private:
   const Context& context_;
   ego object_;

   std::vector<std::shared_ptr<EGADS::Body>> body_;
};

} // EGADS

} // luna

#endif
