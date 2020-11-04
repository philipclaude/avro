#ifndef AVRO_LIB_VORONOI_OPTIMAL_TRANSPORT_H_
#define AVRO_LIB_VORONOI_OPTIMAL_TRANSPORT_H_

#include "element/polytope.h"

#include "mesh/topology.h"
#include "numerics/integration.h"

#include "voronoi/voronoi.h"

#include <memory>
#include <queue>

#include "nnsearch/nn_search.h"

namespace avro
{

namespace delaunay
{

class RVDFacets;

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

class DensityMeasure
{
public:
  DensityMeasure() {}
  virtual ~DensityMeasure() {}

  virtual real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const = 0;
};

class DensityMeasure_Uniform : public DensityMeasure
{
public:
  DensityMeasure_Uniform( real_t rho = 1.0 ) :
    rho_(rho)
  {}
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const { return rho_; }

private:
  real_t rho_;
};

class DensityMeasure_Example : public DensityMeasure
{
public:
  real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const { return x[0]*x[0]; }
};

template<typename type>
class DensityMeasure_Interpolated : public DensityMeasure
{
public:
  DensityMeasure_Interpolated( const Field<type,real_t>& field ) :
    field_(field)
  {}

  real_t evaluate( index_t elem , const real_t* xref , const real_t* p ) const
  {
    // look up the value using the field, with an initial guess of elem for the search
    avro_implement;
    return 0.0;
  }

private:
  const Field<type,real_t>& field_;
};

template<typename type>
class LaguerreCellBase : public Topology<Polytope>
{
protected:
  LaguerreCellBase(const Delaunay& delaunay , GEO::NearestNeighborSearch& nns ,
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
  const Delaunay& delaunay_;        // the delaunay vertices/voronoi sites/dirac masses
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

  std::map<Bisector,int> bisector_;
  std::map<int,Bisector> ids_;

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

  LaguerreCell(index_t site , const Delaunay& delaunay , GEO::NearestNeighborSearch& nns ,
               const Topology<Polytope>& domain , bool exact , index_t nb_nns=50 );

  void initialize();
  void clip();
  void compute();
};

template<>
class LaguerreCell<Simplex> : public LaguerreCellBase<Simplex>
{
public:
  LaguerreCell(index_t elem , const Delaunay& delaunay , GEO::NearestNeighborSearch& nns ,
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



template<typename type>
class LaguerreDiagram : public Topology<Polytope>
{
public:
  typedef LaguerreDiagram thisclass;

  LaguerreDiagram( Delaunay& delaunay , const Topology<type>& domain );

  void initialize();
  void create( bool exact , index_t nb_nns );
  void compute( bool exact , IntegrationSimplices* triangulation=nullptr );

  void clip( const index_t k )
  {
    cells_[k]->set_edges(domain_edges_);
    cells_[k]->compute();
  }

  void set_elements( const std::vector<index_t>& elems ) { elem_ = elems; }

  const std::vector<index_t>& sites() const { return sites_; }
  const LaguerreCell<type>& cell( index_t k ) const { return *cells_[k].get(); }

  real_t time_decompose() const { return time_decompose_; }
  real_t time_neighbours() const { return time_neighbours_; }
  real_t time_voronoi() const { return time_voronoi_; }

private:
  Points points_;                 // voronoi vertices
  Delaunay& delaunay_;            // delaunay vertices/voronoi sites/dirac masses
  const Topology<type>& domain_;  // background mesh

  std::vector<std::shared_ptr<LaguerreCell<type>>> cells_;
  std::vector<index_t> sites_;
  std::vector<index_t> elem_; // an element in the domain containing each delaunay vertex
  GEO::NearestNeighborSearch* neighbour_search_;

  // domain data
  std::shared_ptr<RVDFacets> domain_facets_;
  std::vector<index_t> domain_edges_;

  real_t time_decompose_;
  real_t time_neighbours_;
  real_t time_voronoi_;
};

class TriangulationCells : public Field<Simplex,real_t>
{
public:
  TriangulationCells( delaunay::IntegrationSimplices& t ) :
    Field<Simplex,real_t>(t,0,DISCONTINUOUS)
  {
    this->build();
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    for (index_t k=0;k<t.nb();k++)
    {
      this->value(k) = t.simplex2site(k);
    }
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const
   {std::vector<std::string> result; result.push_back("cells"); return result;}
};

class TriangulationElements : public Field<Simplex,real_t>
{
public:
  TriangulationElements( delaunay::IntegrationSimplices& t ) :
    Field<Simplex,real_t>(t,0,DISCONTINUOUS)
  {
    this->build();
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    for (index_t k=0;k<t.nb();k++)
    {
      this->value(k) = t.simplex2elem(k);
    }
  }

  index_t nb_rank() const { return 1; }

  std::vector<std::string> ranknames() const
   {std::vector<std::string> result; result.push_back("elems"); return result;}
};

class Integrand_Transport_Energy : public Integrand<Integrand_Transport_Energy>
{
public:

  typedef real_t T;

  Integrand_Transport_Energy( const Points& delaunay , const IntegrationSimplices& simplices , const DensityMeasure& density , coord_t dim) :
    delaunay_(delaunay),
    simplices_(simplices),
    density_(density),
    dim_(dim)
  {}

  bool needs_solution() const { return false; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( index_t k , const real_t* xref , const real_t* x ) const;

private:
  const Points& delaunay_;
  const IntegrationSimplices& simplices_;
  const DensityMeasure& density_;
  coord_t dim_;
};

template<typename type>
class SemiDiscreteOptimalTransport
{
public:
  SemiDiscreteOptimalTransport( const Topology<type>& domain , DensityMeasure* density=nullptr );

  void set_density( DensityMeasure* density ) { density_ = density; }
  void set_exact( bool x ) { exact_ = x; }

  void sample( index_t nb_samples );
  void optimize_points( index_t nb_iter=10 );
  void optimize_weights( index_t nb_iter=10 );
  void stochastic_gradient_descent();
  void optimize_points_lloyd( index_t nb_iter=10 );

  void compute_laguerre();
  real_t evaluate( index_t iter , index_t mode , real_t* dc_dx=nullptr , real_t* dc_dw=nullptr );

  void set_delaunay( const real_t* x , coord_t dim );
  void set_weights( const real_t* w );

  const Topology<type>& domain() const { return domain_; }
  Delaunay& delaunay() { return delaunay_; }
  IntegrationSimplices& simplices() { return simplices_; }
  LaguerreDiagram<type>& diagram() { return diagram_; }

  std::vector<real_t>& mass() { return mass_; }
  const std::vector<real_t>& mass() const { return mass_; }
  void set_nu( const std::vector<real_t>& nu ) { nu_ = nu; }
  std::vector<real_t>& nu() { return nu_; }

private:
  const Topology<type>& domain_;
  DensityMeasure* density_;

  Delaunay delaunay_;
  LaguerreDiagram<type> diagram_;
  IntegrationSimplices simplices_;
  bool exact_;

  std::vector<real_t> nu_; // the target mass
  std::vector<real_t> weight_; // the weight on each voronoi cell
  std::vector<real_t> mass_;
  std::vector<real_t> centroid_;
};

} // delaunay

} // avro

#endif
