#ifndef AVRO_API_AVRO_H_
#define AVRO_API_AVRO_H_

#include "avro_types.h"
#include "avro_params.h"

#include <map>
#include <memory>
#include <set>
#include <vector>
#include <iostream>

// forward declaration of ego
struct egObject;
typedef egObject* ego;

namespace avro
{

static const real_t unset_value = 1e20;

// forward declarations of the avro classes
class Entity;
class Model;
class Body;
class Simplex;
template <typename type> class Topology;
class TopologyBase;
class Points;

namespace EGADS {
class Context;
}

class EGADSGeneralGeometry {

public:
  EGADSGeneralGeometry( coord_t number );

  void add_body( ego object );
  void add_object( ego body , ego object );
  void add_child( ego parent , ego child );
  void set_interior( ego object );
  void finalize();

private:
  coord_t number_;
  std::map< ego , std::shared_ptr<Body> > ego2body_;
  std::map< ego , std::shared_ptr<Entity> > ego2entity_;

  std::shared_ptr<EGADS::Context> context_;
};

class Context
{
public:
  Context( coord_t number , coord_t dim , coord_t udim );
  Context( const Context& context );

  coord_t number() const { return number_; }
  coord_t dim() const { return dim_; }
  coord_t udim() const { return udim_; }

  int facet_geometry( const index_t* v , index_t nv ) const;
  void get_geometry_ids( std::vector<int>& ids ) const;
  void get_geometry_idx2ego( std::map<int,ego*>& idx2ego );

  // this should be called once
  void define_geometry( const std::string& geometry );
  void define_geometry( const Model& model );
  void attach_geometry(); // only at the beginning
  void define_geometry( const EGADSGeneralGeometry& geometry );

  // these should be called anytime an adaptation problem needs to be defined
  void define_mesh( const std::string& mesh );
  void load_mesh( const std::vector<real_t>& x , const std::vector<index_t>& s );
  void load_coordinates( const std::vector<real_t>& x );
  void load_simplices( const std::vector<index_t>& s );
  void load_polytopes( const std::vector<index_t>& indices , const std::vector<index_t>& nv_per_elem );
  void load_local2global( const std::vector<index_t>& g );
  void load_geometry( const std::vector<int>& g , const std::vector<real_t>& u );

  // partition function
  void partition();

  // the actual adaptation function
  int adapt( const std::vector<real_t>& metric );
  int adapt_parallel( const std::vector<real_t>& metric );

  // these should be called anytime a mesh needs to be retrieved (e.g. after an adaptation)
  void retrieve_mesh( std::vector<real_t>& x , std::vector<index_t>& s ) const;
  void retrieve_local2global( std::vector<index_t>& local2global ) const;
  void retrieve_geometry( std::vector<int>& g , std::vector<real_t>& u ) const;
  void retrieve_boundary( std::vector<std::vector<index_t>>& faces ,
                          std::vector<int>& geometry , bool interior=false ) const;
  void retrieve_boundary_parallel( std::vector<std::vector<index_t>>& faces ,
                                   std::vector<int>& geometry ) const;

  // some applications require an interface to the neighbours and inverse (e.g. SANS wake geometries)
  void retrieve_neighbours( std::vector<int>& neighbours ) const;
  void retrieve_inverse( std::vector<std::set<index_t>>& inverse ) const;
  ego  get_vertex_ego( index_t k ) const;

  // polytopal mesh retrieval functions
  void retrieve_polytopes( std::vector<real_t>& x , std::vector<index_t>& indices , std::vector<index_t>& nv_per_elem ) const;
  void retrieve_polytope_facets( std::vector<index_t>& indices , std::vector<index_t>& nv_per_elem ) const;
  std::vector<real_t> compute_laguerre( const std::vector<real_t>& sites , const std::vector<real_t>& weights , index_t nb_iter = 0 );
  std::vector<real_t> compute_optimal_transport( const std::vector<real_t>& sites , const std::vector<real_t>& mass , const std::vector<real_t>& initial_weights , index_t nb_iter = 10 );

  ParameterSet& parameters() { return parameters_; }
  const ParameterSet& parameters() const { return parameters_; }

  index_t nb_bodies() const;

  void plot() const;

private:
  Entity* id2geometry( int id ) const;
  void import_model();
  const std::shared_ptr<Model>& model() const { return model_; }

  coord_t number_;
  coord_t dim_;
  coord_t udim_;

  std::shared_ptr<Model> model_;
  std::map<int,Entity*> id2entity_;
  std::map<Entity*,int> entity2id_;
  ParameterSet parameters_;

  std::shared_ptr<Points> points_;
  std::shared_ptr<TopologyBase> topology_;
  std::vector<real_t> metric_;

};

} // avro

#endif
