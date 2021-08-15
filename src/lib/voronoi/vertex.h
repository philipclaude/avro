#ifndef AVRO_LIB_VORONOI_NEW_VERTEX_H_
#define AVRO_LIB_VORONOI_NEW_VERTEX_H_

#include "avro_types.h"

#include "numerics/predicates.h"

#include "voronoi/facets.h"

namespace avro
{

namespace voronoi
{

class Vertex
{

public:
  // constructors/destructor
  Vertex() {}
  explicit Vertex( const coord_t _dim );
  Vertex( const coord_t _dim , const coord_t _number );
  Vertex( const Vertex& v0 );
  Vertex( Vertex&& v );
  Vertex( const std::vector<real_t>& _x , const coord_t _number=0 );
  Vertex( const real_t* x, const coord_t dim );

  Vertex& operator=( const Vertex& rhs) { return *this; }

  coord_t dim() const { return dim_; }
  coord_t number() const { return number_; }
  void set_number( coord_t _number ) { number_ = _number; }

  // coordinates
  const std::vector<real_t>& x() const { return x_; }
  const real_t* X() const { return x_.data(); }
  const real_t& operator[] ( const index_t d ) const { return x_[d]; }
  const real_t& coord( const index_t d ) const { return x_[d]; }
  void set_coordinates( const real_t* x , const coord_t _dim ) {
    dim_ = _dim;
    x_.resize(dim_);
    for (coord_t d = 0; d < dim_; d++)
      x_[d] = x[d];
  }

  // bisector retrieval/addition functions
  index_t nb_bisectors() const { return bisector_.size(); }
  void add_bisector( int b ) { bisector_.push_back(b); }
  const std::vector<int>& bisectors() const { return bisector_; }
  index_t bisector( const index_t k ) const { return bisector_[k]; }

  // simplex retrieval/addition functions
  index_t nb_simplices() const { return simplex_.size(); }
  void add_simplex_vertex( const real_t* v ) { simplex_.push_back(v); }
  const std::vector<const real_t*>& simplices() const { return simplex_; }
  const real_t* simplex( const index_t k ) const { return simplex_[k]; }

  // site addition function
  index_t nb_sites() const { return site_.size(); }
  void add_site( const real_t* zj );
  const std::vector<const real_t*>& sites() const { return site_; }
  const real_t* site( const index_t k ) const { return site_[k]; }

  // intersection functions
  void intersect_geometric( const real_t* q1 , const real_t* q2 , const real_t* p1 , const real_t* p2 );
  void intersect_symbolic( const Vertex* v0 , const Vertex* v1 );
  void intersect_bisectors( const Vertex* v0 , const Vertex* v1 );
  void intersect_simplices( const Vertex* v0 , const Vertex* v1 );

  void set_sites( const Points& delaunay , const std::map<int,Bisector>& B );
  void set_base_site( const real_t* z0 ) { z0_ = z0; }
  void set_delaunay_site( const index_t k , const real_t* z )
    { (k==0) ? z0_ = z : site_[k-1] = z; }

  // side query relative to a bisector
  GEO::Sign side_inexact( const real_t* zi , const real_t *zj );
  GEO::Sign side( const real_t* zi , const real_t* zj , const bool exact = true );

  // print function
  void print( const std::string& pre , const bool symbolic=false ) const;

private:
  coord_t dim_;
  coord_t number_;
  std::vector<real_t>  x_;

  const real_t* z0_;
  std::vector<int>           bisector_;
  std::vector<const real_t*> simplex_;
  std::vector<const real_t*> site_;
};

} // voronoi

} // avro

#endif
