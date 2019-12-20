#ifndef LUNA_LIB_GEOMETRY_EGADS_H_
#define LUNA_LIB_GEOMETRY_EGADS_H_

#include "geometry/egads/body.h"
#include "geometry/egads/data.h"
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

class Context;

class Object : public Entity
{
public:
  Object( const Context& context , ego* object );
  Object( ego* object , EGADS::Body* body );

  void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const;
  void evaluate( const std::vector<real_t>& u , std::vector<real_t>& p ) const;

  void project( std::vector<real_t>& x , std::vector<real_t>& u ) const;

  void set_object( ego* object );

  void build_hierarchy();

  ego* object();
  const ego* object() const;

  const ego egchild( index_t k ) const;

  int object_class() const { return data_.object_class; }


  void print(bool with_children=true) const;

private:
  EGADS::Body* body_;
  const Context& context_;
  ego* object_;

  egoData data_;

  int body_index_;
};

} // EGADS

} // luna

#endif
