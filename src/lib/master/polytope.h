#ifndef avro_LIB_MASTER_POLYTOPE_H_
#define avro_LIB_MASTER_POLYTOPE_H_

#include "common/table.h"

#include "master/master.h"
#include "master/simplex.h"

namespace avro
{

template<typename type> class Topology;
template<typename type> class Triangulation;

class Polytope : public Master<Polytope>
{

public:
  Polytope( coord_t number , coord_t order , const Table<int>& incidence );
  Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence );

  static std::string type_name() { return "polytope"; }

  const Table<int>& incidence() const { return incidence_; }

  void get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const;

  real_t volume( const Points& points , const index_t* v , index_t nv ) const;
  real_t volume( const std::vector<const real_t*>& x , const coord_t dim ) const;

  void facet( const index_t* v , index_t j , std::vector<index_t>& f ) const
    { avro_implement; }

    // for switching between the hrep and vrep
  void hrep( const index_t* v , index_t nv , std::vector<int>& facets ) const;
  void vrep( const index_t* v , index_t nv , int facet , std::vector<index_t>& points ) const;

  // edge retrieval
  void edges( const index_t* v , index_t nv , std::vector<index_t>& e  ) const;
  bool is_edge( index_t v0 , index_t v1 ) const;
  bool is_edge( const std::vector<index_t>& b0 , const std::vector<index_t>& b1 ) const;
  bool is_edge( const std::vector<int>& b0 , const std::vector<int>& b1 ) const;
  index_t edge( index_t k , index_t j ) const
    { avro_assert_not_reached; return 0; }

  void triangulate( const index_t* v , index_t nv , Triangulation<Polytope>& triangulation ) const;

  bool& fullmesh() { return fullmesh_; }

private:
  void triangulate( coord_t number , Table<index_t>& table , Points& points , const index_t* v , index_t nv ) const;

  Simplex simplex_;
  const Table<int>& incidence_;

  bool fullmesh_; // whether the master polytope is responsible for a full mesh
};

} // avro

#endif
