#ifndef luma_LIB_GEOMETRY_EGADS_MODEL_H_
#define luma_LIB_GEOMETRY_EGADS_MODEL_H_

#include "geometry/model.h"

#include "geometry/egads/egads_fwd.h"

#include <string>

namespace luma
{

namespace EGADS
{

class Body;
class Context;

class Model : public luma::Model
{
public:
  Model( const Context& context );
  Model( const Context& context , const std::string& filename , bool split=false );

  void determine_number();

  const Context& context() const { return context_; }

  Entity* find_entity( index_t id , int object_class ) const;

  void print() const;

private:
   const Context& context_;
   ego object_;
};

} // EGADS

} // luma

#endif
