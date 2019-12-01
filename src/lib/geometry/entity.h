#ifndef LUNA_LIB_GEOMETRICS_ENTITY_H_
#define LUNA_LIB_GEOMETRICS_ENTITY_H_

#include "common/tree.h"
#include "common/types.h"

#include <string>

namespace luna
{

namespace numerics
{
  class Coordinate;
}

class Body;

class Entity : public Tree<Entity>
{
public:
  coord_t number() const { return number_; }
  std::string name() const { return name_; }

  index_t identifier() const { return identifier_; }
  Body* body() const { return body_; }

  Entity* intersect( Entity* e1 );
  Entity* intersect( Entity* e1 , Entity* e2 , bool only_check=false );
  Entity* intersect( Entity* e1 , Entity* e2 , Entity* e3 );

  index_t nb_parents() const { return parents_.size(); }
  Entity* parents( index_t k ) { return parents_[k]; }
  const Entity* parent( index_t k ) const { return parents_[k]; }

  const std::vector<Entity*>& parents() const { return parents_; }

  bool tessellatable() const { return tessellatable_; }

  bool above( const Entity* entity ) const;

  bool interior() const { return interior_; }
  void set_interior( bool x ) { interior_ = x; }

  virtual void print() const = 0;

protected:
  Entity( coord_t number );
  Entity( coord_t number , const std::string& name );

  virtual ~Entity() {}

  virtual void inverse( numerics::Coordinate& x , numerics::Coordinate& u ) const = 0;
  virtual void evaluate( const numerics::Coordinate& u , numerics::Coordinate& p ) const = 0;

  bool sense() const { return sense_; }
  void sense( bool x ) { sense_ = x; }

  coord_t number_;
  std::string name_;
  index_t identifier_;
  Body* body_;

  bool interior_;
  bool sense_;

  std::vector<Entity*> parents_;
  bool tessellatable_;
};

} // luna

#endif