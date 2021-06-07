//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_MESH_DELAUNAY_VORONOI_H_
#define AVRO_MESH_DELAUNAY_VORONOI_H_

#include "common/nearest_neighbours.h"

#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/field.hpp"
#include "mesh/mesh.h"
#include "mesh/topology.h"

#include "numerics/predicates.h"

namespace avro
{

class Delaunay;
class RestrictedVoronoiSimplex;
class Entity;

namespace delaunay
{

class Bisector
{
public:
  Bisector( index_t q0 , index_t q1 )
  {
    if (q0<q1)
    {
      p0 = q0;
      p1 = q1;
    }
    else
    {
      p0 = q1;
      p1 = q0;
    }
  }
  int p0;
  int p1;
};

// needed to create a set/map of elements
bool operator==( const Bisector& bx , const Bisector& by );
bool operator<( const Bisector& f , const Bisector& g );

typedef struct
{
  std::vector<int> indices;
} SymbolicVertex;

// needed to create a set/map of elements
inline bool
operator==( const SymbolicVertex& fx , const SymbolicVertex& fy )
{
  // assumes fx and fy have the same topological dimension
  // and that the indices are sorted
  avro_assert( fx.indices.size()==fy.indices.size() );
  for (index_t j=0;j<fx.indices.size();j++)
    if (fx.indices[j]!=fy.indices[j])
      return false;
  return true;
}

// needed to create a map of elements
inline bool
operator<( const SymbolicVertex& f , const SymbolicVertex& g )
{
  // lexicographically compare the indices
  return std::lexicographical_compare(f.indices.begin(), f.indices.end(),
                                      g.indices.begin(), g.indices.end());
}

class Vertex
{
  friend class VoronoiVertex_tester;

public:
  // constructors/destructor
  Vertex() {}
  explicit Vertex( const coord_t _dim );
  Vertex( const coord_t _dim , const coord_t _number );
  Vertex( const Vertex& v0 );
  Vertex(Vertex&& v);
  Vertex( const std::vector<real_t> _x , const coord_t _number=0 );
  Vertex( const real_t* x, const coord_t dim );

  Vertex& operator=( const Vertex& rhs) { return *this; }

  // coordinate/topology set/get functions
  void init();
  coord_t dim() const { return dim_; }
  coord_t number() const { return number_; }
  void setNumber( coord_t _number ) { number_ = _number; }
  const std::vector<real_t>& x() const { return x_; }
  const real_t* X() const { return x_.data(); }
  const real_t& operator[] ( const index_t d ) const { return x_[d]; }
  const real_t& coord( const index_t d ) const { return x_[d]; }
  void setCoordinates( const real_t* x , const coord_t _dim )
  {
    dim_ = _dim;
    init();
    for (coord_t d=0;d<dim_;d++)
      x_[d] = x[d];
  }

  // bisector retrieval/addition functions
  void addBisector( int b ) { bisector_.push_back(b); }
  const std::vector<int>& bisectors() const { return bisector_; }
  index_t bisector( const index_t k ) const { return bisector_[k]; }
  index_t nb_bisectors() const { return bisector_.size(); }

  // topology retrieval/addition functions
  void addTopology( Topology<Simplex>* m ) { topology_.push_back(m); }
  const std::vector< Topology<Simplex>* >& topologies() const { return topology_; }
  Topology<Simplex>* topology( const index_t k ) const { return topology_[k]; }

  // simplex retrieval/addition functions
  void addSimplexVertex( const real_t* v ) { simplex_.push_back(v); }
  const std::vector<const real_t*>& simplices() const { return simplex_; }
  index_t nb_simplices() const { return simplex_.size(); }
  const real_t* simplex( const index_t k ) const { return simplex_[k]; }

  // site addition function
  void addSite( const real_t* zj );
  std::vector<const real_t*> sites() const { return site_; }
  const real_t* site( const index_t k ) const { return site_[k]; }
  index_t nb_sites() const { return site_.size(); }

  // intersection functions
  void intersectSymbolic( const Vertex* v0 , const Vertex* v1 ,
      const Delaunay& delaunay );
  void intersectGeometric( const real_t* q1 , const real_t* q2 , const real_t* p1 , const real_t* p2 );
  void intersectBisectors( const Vertex* v0 , const Vertex* v1 );
  void intersectMeshes( const Vertex* v0 , const Vertex* v1 );
  void intersectSimplices( const Vertex* v0 , const Vertex* v1 );
  void setSites( const Delaunay& delaunay );
  void setSites( const Delaunay& delaunay , const std::map<int,Bisector>& B );
  void setBaseSite( const real_t* z0 ) { z0_ = z0; }
  void setDelaunaySite( const index_t k , const real_t* z )
    { (k==0) ? z0_ = z : site_[k-1] = z; }

  // side query relative to a bisector
  GEO::Sign sideFast( const real_t* zi , const real_t *zj );
  GEO::Sign side( const real_t* zi , const real_t* zj , const bool exact = true );
  GEO::Sign side_meshless( const real_t* zi , const real_t* zj , const bool exact = true );

  // print function
  void print( const std::string& pre , const bool symbolic=false ) const;

private:
  coord_t dim_;
  coord_t number_;
  std::vector<real_t>  x_;

  const real_t* z0_;
  std::vector<int> bisector_;
  std::vector< Topology<Simplex>* > topology_;
  std::vector<const real_t*> simplex_;
  std::vector<const real_t*> site_;
};

class RVDFacets
{
public:
  RVDFacets( const Topology<Simplex>& topology );

  void create();
  int facet( const std::vector<index_t>& f ) const;
  void print() const;

private:
  std::string generate( const std::vector<index_t>& f ) const;
  int lookup( const std::string& s , int& id ) const;

  std::map<std::string,int> store_;
  const Topology<Simplex>& topology_;

};


class RestrictedVoronoiSimplex : public Topology<Polytope>
{
  friend class RestrictedSimplex_tester;

public:
  RestrictedVoronoiSimplex(
    const index_t k , const Topology<Simplex>& mesh , RVDFacets& facets ,
    Delaunay& _delaunay , NearestNeighbours& _neighbours , bool _exact );

  void set_entity( Entity* entity );

  void reset();

  // clip the simplex by the entire Voronoi diagram defined by delaunay_
  void clip();

  // clip the simplex by the cell generated by seed i of the voronoi diagram
  void clip( const index_t i );

  // clip the current polytope by the bisector zi_ to neighbour_[zi][j]
  void clipPolytope( index_t j );

  // clip the edge e0-e1 by the bisector b into the clipped polytope q
  void clipEdge( const index_t e0 , const index_t e1 , const int b ,
    std::vector<index_t>& q );

  // delaunay site functions
  index_t nb_sites() const;
  void nextSite();
  index_t site() const { return site_; }
  void addSite( const index_t zj );

  bool securityRadiusReached( const real_t* zi , const real_t* zj ) const;

  // function to send the points into their location after merging and stuff
  void finalize();

  index_t seed( const index_t k ) const { return seed_[k]; }

  Topology<Polytope>& topology() { return topology_; }

private:
  index_t site_;
  std::vector<index_t> sites_;
  std::vector<bool>    clipped_;

  // yes, this is a vector of objects, but the invocation of the copy
  // constructor is avoided by pushing a default VoronoiVertex
  // and then populating its fields from the created one
  std::vector<Vertex> vertex_;

  std::vector<index_t> simplex_;   // original simplex points
  std::vector<index_t> polytope_;  // current polytope points

  // the delaunay neighbours for each cell in this RVS
  Table<index_t> region_;

  const RVDFacets& facets_;
  const Delaunay& delaunay_;
  const NearestNeighbours& neighbours_;
  const bool exact_;

  std::vector<index_t> seed_;

  Topology<Polytope> topology_;

  Points points_;

  Entity* entity_;
};

class VoronoiSites;

class RestrictedVoronoiDiagram : public Topology<Polytope>
{
public:

  typedef RestrictedVoronoiDiagram thisclass;

  RestrictedVoronoiDiagram( const Topology<Simplex>& _mesh , Delaunay& _delaunay , Entity* entity=nullptr );
  RestrictedVoronoiDiagram( Delaunay& _delaunay );

  void compute( const bool exact = true );

  void accumulate();

  index_t nb_simplices() const { return simplices_.size(); }
  RestrictedVoronoiSimplex* simplex( const index_t k ) const
    { return simplices_[k].get(); }

  real_t compute_centroids( Points& centroids );

  real_t energy();
  real_t energy() const { return energy_; }

  real_t energy_nd();

  bool& parallel() { return parallel_; }
  bool& gpu() { return gpu_; }

  void optimise( const index_t nb_iter , bool exact=false );

  std::string& outdir() { return outdir_; }

  void clip( const index_t k )
    { simplices_[k]->clip(); }

  void extract( Topology<Simplex>& triangulation ) const;

  VoronoiSites& sites() { return *sites_.get(); }
  std::shared_ptr<VoronoiSites>& sites_ptr() { return sites_; }

private:

  //Points points_;
  Points vertices_;

  const Topology<Simplex>& mesh_;
  Delaunay& delaunay_;
  NearestNeighbours neighbours_;

  std::vector<std::shared_ptr<RestrictedVoronoiSimplex>> simplices_;

  bool parallel_;
  bool gpu_;
  real_t energy_;
  std::string outdir_;

  std::shared_ptr<VoronoiSites> sites_;

  Entity* entity_;

};

class VoronoiSites : public Field<Polytope,real_t>
{
public:
  VoronoiSites( Topology<Polytope>& rvd ) :
    Field<Polytope,real_t>(rvd,0,DISCONTINUOUS)
  {}

  void add_cell( const real_t z )
  {
    //Field<real_t>::add(z+1);
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const
   {std::vector<std::string> result; result.push_back("sites"); return result;}
};

} // delaunay

} // avro

#endif
