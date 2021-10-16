#ifndef AVRO_LIB_VORONOI_OPTIMAL_TRANSPORT_H_
#define AVRO_LIB_VORONOI_OPTIMAL_TRANSPORT_H_

#include "element/polytope.h"

#include "mesh/field.hpp"
#include "mesh/points.h"
#include "mesh/topology.h"
#include "numerics/integration.h"

//#include "voronoi/vertex.h"
#include "voronoi/cell_nd.h"

#include <memory>
#include <queue>

#include <geogram/nn_search.h>

namespace avro
{

namespace voronoi
{

class RVDFacets;

class DensityMeasure
{
public:
  DensityMeasure() {}
  virtual ~DensityMeasure() {}

  virtual real_t evaluate( index_t elem , const real_t* xref , const real_t* x ) const = 0;

  const std::string& name() const { return name_; }

protected:
  std::string name_;
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
class LaguerreDiagram : public Topology<Polytope>
{
public:
  typedef LaguerreDiagram thisclass;

  LaguerreDiagram( Points& delaunay , const Topology<type>& domain );
  ~LaguerreDiagram();

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

  index_t minimum_neighbours() const { return minimum_neighbours_; }
  index_t maximum_neighbours() const { return maximum_neighbours_; }
  real_t average_neighbours() const { return average_neighbours_; }
  const std::vector<index_t>& neighbour_counts() const {return neighbour_counts_; }

private:
  void compute_neighbour_properties();
  std::vector<index_t> neighbour_counts_;

  Points points_;                 // voronoi vertices
  Points& delaunay_;            // delaunay vertices/voronoi sites/dirac masses
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

  index_t minimum_neighbours_;
  index_t maximum_neighbours_;
  real_t average_neighbours_;
};

class TriangulationCells : public Field<Simplex,real_t>
{
public:
  TriangulationCells( voronoi::IntegrationSimplices& t ) :
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
  TriangulationElements( voronoi::IntegrationSimplices& t ) :
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

  T operator()( index_t k , const QuadraturePoint& point , const real_t* x ) const;

private:
  const Points& delaunay_;
  const IntegrationSimplices& simplices_;
  const DensityMeasure& density_;
  coord_t dim_;
};


struct SDOT_Properties
{
  std::vector<real_t> mass_min;
  std::vector<real_t> mass_max;
  std::vector<real_t> time_voronoi;
  std::vector<real_t> time_neighbours;
  std::vector<real_t> time_integration;
  std::vector<real_t> time_triangulation;
  std::vector<real_t> neighbours_average;
  std::vector<index_t> neighbours_minimum;
  std::vector<index_t> neighbours_maximum;
  std::vector<real_t> energy;
  std::vector<real_t> gradient;

  void save( const std::string& filename ) const;
  void clear()
  {
    mass_min.clear();
    mass_max.clear();
    time_voronoi.clear();
    time_neighbours.clear();
    time_integration.clear();
    time_triangulation.clear();
    neighbours_average.clear();
    neighbours_maximum.clear();
    energy.clear();
    gradient.clear();
  }
};

struct SDOT_Snapshot
{
  coord_t dim;
  std::vector<real_t> weights;
  std::vector<real_t> points;
  std::vector<index_t> neighbours;
  std::vector<real_t> mass;
  std::string name;

  void save( const std::string& filename ) const;
  void clear()
  {
    weights.clear();
    points.clear();
  }
};

class OptimalTransportBase
{
public:
  virtual real_t transport_objective( index_t n , const real_t* x , real_t* grad ) = 0;
  virtual ~OptimalTransportBase() {}

  index_t& iteration() { return iteration_; }

  virtual void sample( index_t nb_sample ) = 0;
  virtual void optimize_points(index_t) = 0;
  virtual void optimize_points_lloyd(index_t) = 0;
  virtual void optimize_weights(index_t) = 0;
  virtual void set_delaunay( const real_t* x , coord_t dim ) = 0;
  virtual void set_weights( const real_t* w ) = 0;
  virtual void set_nu( const std::vector<real_t>& nu ) = 0;

  virtual const Topology<Polytope>& get_diagram() const = 0;
  virtual std::vector<real_t> get_sites() const = 0;
  virtual std::vector<real_t> get_weights() const = 0;

protected:
  index_t iteration_;
};

template<typename type>
class SemiDiscreteOptimalTransport : public OptimalTransportBase
{
public:
  SemiDiscreteOptimalTransport( const Topology<type>& domain , DensityMeasure* density=nullptr );

  void set_density( DensityMeasure* density ) { density_ = density; }
  void set_exact( bool x ) { exact_ = x; }

  // applications
  void sample( index_t nb_samples );
  void optimize_points( index_t nb_iter=10 );
  void optimize_weights( index_t nb_iter=10 );
  void stochastic_gradient_descent( index_t nb_iter=10 );
  void optimize_points_lloyd( index_t nb_iter=10 );
  void generate_bluenoise();

  void compute_laguerre();
  real_t evaluate( real_t* dc_dx=nullptr , real_t* dc_dw=nullptr );
  real_t transport_objective( index_t n , const real_t* x , real_t* grad );

  void set_delaunay( const real_t* x , coord_t dim );
  void set_weights( const real_t* w );

  const Topology<type>& domain() const { return domain_; }
  Points& delaunay() { return delaunay_; }
  LaguerreDiagram<type>& diagram() { return diagram_; }
  IntegrationSimplices& simplices() { return simplices_; }

  std::vector<real_t>& mass() { return mass_; }
  const std::vector<real_t>& mass() const { return mass_; }
  void set_nu( const std::vector<real_t>& nu ) { nu_ = nu; }
  std::vector<real_t>& nu() { return nu_; }
  const std::vector<real_t>& weights() const { return weight_; }
  std::vector<real_t> get_sites() const;
  std::vector<real_t> get_weights() const { return weight_; }

  void start();
  int mode() const { return mode_; }
  void set_mode( int x ) { mode_ = x; }

  index_t& quad_order() { return quad_order_; }
  real_t& weight_max() { return weight_max_; }
  real_t time_integrate() const { return time_integrate_; }

  const SDOT_Properties& properties() const { return properties_; }

  void save_every( index_t x , const std::string& prefix )
  {
    save_every_ = x;
    prefix_     = prefix;
  }
  void save_snapshot();

  const Topology<Polytope>& get_diagram() const { return diagram_; }

private:
  const Topology<type>& domain_;
  DensityMeasure* density_;

  Points delaunay_;
  LaguerreDiagram<type> diagram_;
  IntegrationSimplices simplices_;
  bool exact_;
  int mode_;

  std::vector<real_t> nu_; // the target mass
  std::vector<real_t> weight_; // the weight on each voronoi cell
  std::vector<real_t> mass_;
  std::vector<real_t> centroid_;

  bool print_;
  index_t quad_order_;
  real_t weight_max_;
  index_t save_every_;
  std::string prefix_;

  SDOT_Properties properties_;
  SDOT_Snapshot snapshot_;
  real_t time_integrate_;
};

} // voronoi

} // avro

#endif
