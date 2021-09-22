//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_MESH_points_H_
#define avro_MESH_points_H_

#include "common/array.h"
#include "common/table.h"
#include "common/error.h"
#include "common/json.h"
#include "avro_types.h"

#include "mesh/dof.h"

#include <memory>
#include <vector>

#if AVRO_MPI
#include "common/mpi.hpp"
#endif

namespace avro
{

class Entity;
class Body;
class Model;

class Points : public DOF<real_t>
{
public:
  Points();
  Points( const coord_t dim , const coord_t udim );
  Points( const coord_t dim );
  ~Points();

  void copy( Points& v , const bool ghosts=true ) const;

  coord_t dim() const { return dim_; }
  coord_t udim() const { return udim_; }
  index_t nb() const { return DOF<real_t>::nb(); }

  void set_dim( coord_t _dim )
    { dim_ = _dim; DOF<real_t>::set_rank(dim_); }
  void set_parameter_dim( coord_t _udim )
    { udim_ = _udim; u_.set_rank(udim_); }

  // vertex creation
  void create( const std::vector<real_t>& x );
  void create( const real_t* x );
  void create_ghost();

  // parameter coordinate access
  const real_t* u( const index_t k ) const { return u_[k]; }
  real_t* u( const index_t k ) { return u_[k]; }
  real_t& u( const index_t k , const index_t d )
    { return u_(k,d); }

  void print( bool info=false ) const;
  void print( index_t k , bool info=false ) const;

	void clear();

  int& body( const index_t k );
  void set_entity( const index_t k , Entity* e );
  void set_param( const index_t k , const std::vector<real_t>& u );
  void set_param( const index_t k , const real_t* u );
  void set_fixed( const index_t k , bool f ) { fixed_.set(k,f); }
  void set_global( const index_t k , index_t g ) { global_.set(k,g); }
  void set_age( const index_t k , index_t a ) { age_.set(k,a); }

  Entity* entity( const index_t k ) const { return primitive_[k]; }
  bool fixed( const index_t k ) const { return fixed_[k]; }
  index_t global( const index_t k ) const { return global_[k]; }
  index_t age( const index_t k ) const { return age_[k]; }
  index_t& age( const index_t k ) { return age_[k]; }

  // boundary vertex query function
  bool boundary( const index_t k ) const;

  bool ghost( const index_t k ) const { return k<nb_ghost_; }
  index_t nb_ghost() const { return nb_ghost_; }

  void remove( const index_t k );

  Table<int>& incidence() { return incidence_; }
  const Table<int>& incidence() const { return incidence_; }

  // duplicates detection (useful when the mesh may have come from an OBJ file)
  // or given a vertex-facet incidence matrix F so that merging is done symbolically
  void duplicates( std::vector<index_t>& idx ,real_t tol=1e-16 ) const;
  void duplicates( std::vector<index_t>& idx , const Table<int>& F ) const;

  void dump( const std::string& filename ) const;

  void attach( const Body& body , index_t ibody=1 ,real_t tol=1e-12 );
  void attach( const Model& model ,real_t tol=1e-12 );

  void to_json( json& J ) const;
  void from_json( const json& J , const Model* model=NULL );

  // computing parameters from inverse evaluations
  void compute_param( index_t k );
  void compute_params();

  // extracting parameters for a point on an entity
  void extract_params( index_t k , Entity* face , real_t* u ) const;

  real_t INFTY = 1.0e20;
  index_t GLOBAL_UNSET = 0;

  bool parameter_space() const { return parameter_space_; }
  void set_parameter_space( bool x ) { parameter_space_ = x; }

  #if AVRO_MPI
  void send( mpi::communicator& comm , index_t destination ) const;
  void receive( mpi::communicator& comm , index_t sender );
  #endif

  void move_to( index_t k0 , index_t k1 );
  void move_to_front( const std::vector<index_t>& idx );

  void reserve( index_t n );
  void shrink_to_fit();
  void batch_erase( index_t n );

protected:

  coord_t dim_;  // dimension of ambient space the vertices live in
  coord_t udim_; // dimension of the parameter space of the geometry

  DOF<real_t>    u_;                           // geometry parameter coordinates
  Array<Entity*> primitive_; // which geometric primitive this vertex is on
  Array<int>     body_;      // which geometry body this vertex is on
  Array<bool>    fixed_;     // whether this vertex is tagged as fixed
  Array<index_t> global_;    // global id of this vertex (used for parallel algorithms), 1-bias since 0 = GLOBAL_UNSET
  Table<int>     incidence_; // vertex-facet/bisector incidence matrix
  Array<index_t> age_;

  index_t nb_ghost_; // how many ghost vertices

  bool parameter_space_;
};

} // avro

#endif
