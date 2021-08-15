#ifndef AVRO_LIB_VORONOI_CELL_ND_H_
#define AVRO_LIB_VORONOI_CELL_ND_H_

#include "element/polytope.h"

#include "mesh/points.h"
#include "mesh/topology.h"

#include "voronoi/vertex.h"

#include <geogram/nn_search.h>

namespace avro
{

namespace voronoi
{

class IntegrationSimplices : public Topology<Simplex>
{
public:
  IntegrationSimplices( coord_t number , coord_t dim ) :
    Topology<Simplex>(points_,number),
    points_(dim)
  {}

  void add_simplex( const std::vector<index_t>& simplex , index_t elem , index_t site )
  {
    Topology<Simplex>::add( simplex.data() , simplex.size() );
    simplex2elem_.push_back(elem);
    simplex2site_.push_back( site );
  }

  void add_point( const real_t* x , index_t elem , index_t site )
  {
    points().create(x);
    point2elem_.push_back(elem);
    point2site_.push_back(site);
  }

  index_t simplex2elem( index_t k ) const { return simplex2elem_[k]; }
  index_t simplex2site( index_t k ) const { return simplex2site_[k]; }
  index_t point2elem( index_t k ) const { return point2elem_[k]; }
  index_t point2site( index_t k ) const { return point2site_[k]; }

  const std::vector<index_t>& simplex2site() const { return simplex2site_; }

  void clear()
  {
    Topology<Simplex>::clear();
    points_.clear();
    point2elem_.clear();
    point2site_.clear();
    simplex2elem_.clear();
    simplex2site_.clear();
  }

private:
  Points points_;
  std::vector<index_t> point2elem_;   // triangulation point inside which domain element?
  std::vector<index_t> point2site_;   // triangulation point inside which voronoi cell?
  std::vector<index_t> simplex2elem_;  // triangulation cell inside which domain element?
  std::vector<index_t> simplex2site_; // triangulation cell inside which voronoi cell?
};

template<typename type>
class LaguerreCellBase : public Topology<Polytope>
{
protected:
  LaguerreCellBase(const Points& delaunay , GEO::NearestNeighborSearch& nns ,
               const Topology<type>& domin , bool exact , index_t nb_nns=50 );

  // for passing initial guess back to neighbour reconstruction
  index_t nb_neighbours() const { return neighbours_.size(); }
  const std::vector<index_t>& neighbours() const { return neighbours_; }

  // decomposition-related functions
  void generate_simplices();

public:
  const IntegrationSimplices& simplices() const { return simplices_; }
  void set_edges( const std::vector<index_t>& edges ) { domain_edges_ = edges; }
  void set_facets( RVDFacets* facets ) { domain_facets_ = facets; }

  void get_bisector( int b , index_t& p0 , index_t& p1 ) const;

  index_t cell2site( index_t j ) const { return cell2site_[j]; }

  real_t time_neighbours() const { return time_neighbours_; }
  real_t time_clip() const { return time_clip_; }
  real_t time_decompose() const { return time_decompose_; }
  void print() const;

protected:

  void reset();
  void enlarge_neighbours();

  void set_site( index_t site ) { site_ = site; }
  void add_site( index_t site ) { sites_.push_back(site); }
  void clip_by_bisector( index_t j , index_t bj );
  int  clip_edge( index_t e0 , index_t e1 , const int b , std::vector<index_t>& q ,int& q0, int& q1  );
  bool security_radius_reached( index_t bj ) const;
  int  add_bisector( index_t p0 , index_t p1 );

protected:
  Points points_;                   // points stored for voronoi vertices
  const Points& delaunay_;        // the delaunay vertices/voronoi sites/dirac masses
  GEO::NearestNeighborSearch& neighbour_search_; // search structure through delaunay points
  std::vector<index_t> neighbours_;         // current list of nearest neighbours
  const Topology<type>& domain_;    // domain over which we compute the diagram
  const bool exact_;                // whether to be in exact or inexact mode
  index_t nb_neighbours_;           // number of nearest neighbors to start off with in the search structure

  index_t site_;
  std::vector<index_t> sites_; // sites to process

  std::vector<index_t> cell2elem_;
  std::vector<index_t> cell2site_;
  std::vector<index_t> point2elem_;
  std::vector<index_t> point2site_;

  // domain data (facets & edges)
  RVDFacets* domain_facets_;
  std::vector<index_t> domain_edges_;

  std::map<voronoi::Bisector,int> bisector_;
  std::map<int,voronoi::Bisector> ids_;

  // decomposition-related structures
  IntegrationSimplices simplices_;

  // polytope manipulation structures
  std::vector<Vertex>  vertex_;    // list of vertices
  std::vector<index_t> polytope_;  // current polytope points
  std::vector<index_t> qpolytope_;
  std::vector<index_t> pedges_;
  std::vector<index_t> qedges_;
  std::vector<index_t> qplane_;

  // timing stuff
  real_t time_neighbours_;
  real_t time_clip_;
  real_t time_decompose_;
};

template<typename type> class LaguerreCell;

template<>
class LaguerreCell<Polytope> : public LaguerreCellBase<Polytope>
{
public:

  LaguerreCell(index_t site , const Points& delaunay , GEO::NearestNeighborSearch& nns ,
               const Topology<Polytope>& domain , bool exact , index_t nb_nns=50 );

  void initialize();
  void clip();
  void compute();
};

template<>
class LaguerreCell<Simplex> : public LaguerreCellBase<Simplex>
{
public:
  LaguerreCell(index_t elem , const Points& delaunay , GEO::NearestNeighborSearch& nns ,
               const Topology<Simplex>& domin , RVDFacets* facets , bool exact , index_t nb_nns=50 );

   void clip();
   void clip( index_t i );
   void compute();

private:
   void next_site();
   index_t nb_sites() const;

private:
  index_t elem_;
  std::vector<bool> clipped_;
};


} // voronoi

} // avro


#endif
