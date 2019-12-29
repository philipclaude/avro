#ifndef luma_LIB_GEOMETRICS_ENTITY_H_
#define luma_LIB_GEOMETRICS_ENTITY_H_

#include "common/tree.h"
#include "common/types.h"

#include <string>

namespace luma
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
  void set_parent( Entity* parent ) { parent_ = parent; }

  const std::vector<Entity*>& parents() const { return parents_; }

  bool tessellatable() const { return tessellatable_; }

  bool above( const Entity* entity ) const;

  bool interior() const { return interior_; }
  void set_interior( bool x ) { interior_ = x; }

  bool sense_required() const { return sense_required_; }
  void set_sense_required( bool x ) { sense_required_ = x; }

  virtual void print(bool with_children=false) const = 0;
  virtual void build_hierarchy() = 0;

  virtual void inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const = 0;
  virtual void inverse_guess( std::vector<real_t>& x , std::vector<real_t>& u ) const = 0;
  virtual void evaluate( const std::vector<real_t>& u , std::vector<real_t>& p ) const = 0;

protected:
  Entity( coord_t number );
  Entity( coord_t number , const std::string& name );

  virtual ~Entity() {}

  coord_t number_;
  std::string name_;
  index_t identifier_;
  Body* body_;

  bool interior_;
  bool sense_required_;

  Entity* parent_;
  bool tessellatable_;
};

} // luma

#endif
