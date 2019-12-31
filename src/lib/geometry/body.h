#ifndef avro_LIB_GEOMETRY_BODY_H_
#define avro_LIB_GEOMETRY_BODY_H_

#include "common/types.h"

#include <memory>
#include <vector>

namespace avro
{

class Entity;

class Body
{
protected:
  typedef std::shared_ptr<Entity> Entity_ptr;
public:
  virtual ~Body() {}

  coord_t number() const { return number_; }
  void add( Entity_ptr prim );

  index_t nb_entities() const { return entity_.size(); }

  void build_parents();

  virtual void print() const = 0;


  void get_entities( std::vector<Entity*>& entities ) const;
  void get_tessellatable( std::vector<Entity*>& entities ) const;

protected:

  Body( coord_t number );

  std::vector<Entity_ptr> entity_;

  coord_t number_;
};

} // avro

#endif
