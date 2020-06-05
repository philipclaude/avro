//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRICS_ENTITY_H_
#define avro_LIB_GEOMETRICS_ENTITY_H_

#include "common/tree.h"
#include "common/types.h"

#include <string>

namespace avro
{

class Body;
class BodyTessellation;
template<typename type> class Topology;

class Entity : public Tree<Entity>
{
public:
  coord_t number() const { return number_; }
  std::string name() const { return name_; }

  void set_identifier( index_t id ) { identifier_ = id; }
  index_t identifier() const { return identifier_; }
  Body* body() const { return body_; }
  void set_body( Body* body ) { body_ = body; }

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

  int sign() const { return sign_; }

  bool egads() const { return egads_; }

  template<typename type> void construct( std::shared_ptr<Topology<type>>& node , Topology<type>& root ) const;

  void print_header() const;

protected:
  virtual ~Entity() {}
  Entity( coord_t number );
  Entity( coord_t number , const std::string& name );

  coord_t number_;
  std::string name_;
  Entity* parent_;
  int identifier_;
  Body* body_;

  bool interior_;
  bool sense_required_;
  bool tessellatable_;

  bool egads_;

  int sign_;
};

} // avro

#endif
