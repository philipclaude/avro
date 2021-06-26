#ifndef AVRO_API_AVRO_H_
#define AVRO_API_AVRO_H_

#include "avro_types.h"
#include "avro_params.h"

#include <memory>
#include <vector>
#include <iostream>

namespace avro
{

static const real_t unset_value = 1e20;

// forward declarations of the avro classes
class Entity;
class Model;
class Simplex;
template <typename type> class Topology;
class Points;

// forward declaration of ego
struct egObject;
typedef egObject* ego;

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

  // polytopal mesh retrieval functions
  void retrieve_polytopes( std::vector<index_t>& indices , std::vector<index_t>& nv_per_elem ) const;
  void retrieve_polytope_facets( std::vector<index_t>& indices , std::vector<index_t>& nv_per_elem ) const;
  void compute_laguerre( const std::vector<real_t>& sites , const std::vector<real_t>& weights );
  void compute_optimal_transport( const std::vector<real_t>& sites , const std::vector<real_t>& mass , const std::vector<real_t>& initial_weights );

  ParameterSet& parameters() { return parameters_; }
  const ParameterSet& parameters() const { return parameters_; }

  index_t nb_bodies() const;

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
  std::shared_ptr<Topology<Simplex>> topology_;
  std::vector<real_t> metric_;

};

} // avro

#endif
