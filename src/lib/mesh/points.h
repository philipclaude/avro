#ifndef LUNA_MESH_points_H_
#define LUNA_MESH_points_H_

#include "common/array.h"
#include "common/table.h"
#include "common/error.h"
#include "common/json.h"
#include "common/types.h"

#include "mesh/dof.h"

#include <memory>
#include <vector>

namespace luna
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

  Entity* entity( const index_t k ) const { return primitive_[k]; }
  bool fixed( const index_t k ) const { return fixed_[k]; }

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

  void compute_param( index_t k );
  void compute_params();

  real_t INFTY = 1.0e20;

protected:

  coord_t dim_;  // dimension of ambient space the vertices live in
  coord_t udim_; // dimension of the parameter space of the geometry

  DOF<real_t>    u_;                           // geometry parameter coordinates
  Array<Entity*> primitive_; // which geometric primitive this vertex is on
  Array<int>     body_;      // which geometry body this vertex is on
  Array<bool>    fixed_;     // whether this vertex is tagged as fixed
  Table<int>     incidence_; // vertex-facet/bisector incidence matrix

  index_t nb_ghost_; // how many ghost vertices
};

} // luna

#endif
