#include "common/parallel_for.h"

#include "element/quadrature.h"

#include "numerics/geometry.h"
#include "numerics/integration_rank.h"
#include "numerics/nlopt_result.h"

#include "voronoi/delaunay.h"
#include "voronoi/optimal_transport.h"
#include "voronoi/voronoi_cell.h"

#include <nlopt.hpp>
#include <HLBFGS/HLBFGS.h>
#include "json/json.hpp"

#include <fstream>
#include <iomanip>

namespace avro
{

namespace delaunay
{

template<typename type>
LaguerreCellBase<type>::LaguerreCellBase( const Delaunay& delaunay, GEO::NearestNeighborSearch& nns,
                          const Topology<type>& domain, bool exact, index_t nb_nns ) :
  Topology<Polytope>(points_,domain.number()),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  neighbour_search_(nns),
  neighbours_(nb_nns,0),
  domain_(domain),
  exact_(exact),
  site_(delaunay.nb()),
  domain_facets_(nullptr),
  simplices_(domain.number(),domain.number())
{}

template<typename type>
void
LaguerreCellBase<type>::reset() {

  // reset the mesh data
  clear();
  points_.clear();
  simplices_.clear();

  // reset the timing data
  time_neighbours_ = 0;
  time_clip_       = 0;
  time_decompose_  = 0;

  // the bisectors need to be cleared everytime we reclip
  bisector_.clear();
  ids_.clear();
  sites_.clear();

  cell2elem_.clear();
  cell2site_.clear();
  point2elem_.clear();
  point2site_.clear();
}

template<typename type>
void
LaguerreCellBase<type>::enlarge_neighbours() {

  index_t nb_nns = neighbours_.size();
  nb_nns += 10;
  if (nb_nns > delaunay_.nb()) nb_nns = delaunay_.nb();
  std::vector<real_t> dist2(nb_nns); // not used
  neighbours_.resize( nb_nns );
  clock_t t0 = clock();
  neighbour_search_.get_nearest_neighbors( nb_nns , delaunay_[site_] , neighbours_.data() , dist2.data()  );
  time_neighbours_ += real_t(clock()-t0)/real_t(CLOCKS_PER_SEC);
}

LaguerreCell<Polytope>::LaguerreCell(index_t site , const Delaunay& delaunay ,
                                     GEO::NearestNeighborSearch& nns ,
                                     const Topology<Polytope>& domain , bool exact , index_t nb_nns ) :
  LaguerreCellBase(delaunay,nns,domain,exact,nb_nns) {
  set_site(site);
}

void
LaguerreCell<Polytope>::initialize() {

  index_t elem_ = 0;
  avro_assert( domain_.nb() == 1 );

  // create the initial points
  vertex_.clear();
  vertex_.resize( domain_.nv(elem_) );
  for (index_t j = 0; j < domain_.nv(elem_); j++) {
    vertex_[j] = Vertex(delaunay_.dim(),domain_.number());
    vertex_[j].setCoordinates( domain_.points()[ domain_(elem_)[j] ] , domain_.points().dim() );
    vertex_[j].setNumber( domain_.number() );
    vertex_[j].addSimplexVertex( domain_.points()[ domain_(elem_,j) ] );
    std::vector<int> b = domain_.points().incidence().get( domain_(elem_,j) );
    for (index_t i=0;i<b.size();i++)
      vertex_[j].addBisector(b[i]);
  }

  polytope_ = linspace(domain_.nv(elem_));

  // initialize the edges
  if (domain_edges_.size() > 0) {
    pedges_ = domain_edges_;
  }
  else {
    pedges_.clear();
    for (index_t ii=0;ii<polytope_.size();ii++)
    for (index_t jj=ii+1;jj<polytope_.size();jj++) {
      index_t e0 = polytope_[ii];
      index_t e1 = polytope_[jj];
      if (element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) ) {
        // clip the edge and save the result into q
        pedges_.push_back(e0);
        pedges_.push_back(e1);
      }
    }
  }
}

void
LaguerreCell<Polytope>::compute() {

  clock_t t0;

  // reset everything
  reset();
  initialize();

  // initialize the neighbours
  std::vector<double> dist2(neighbours_.size(),0.0);
  t0 = clock();
  neighbour_search_.get_nearest_neighbors( neighbours_.size() , site_ , neighbours_.data() , dist2.data() );
  avro_assert( neighbours_[0] == site_ );
  time_neighbours_ += real_t(clock() -t0) / real_t(CLOCKS_PER_SEC);

  // do the clipping
  t0 = clock();
  clip();
  time_clip_ += real_t(clock() - t0)/real_t(CLOCKS_PER_SEC);

  // option to decompose the simplices
  t0 = clock();
  if (number_ < 5) generate_simplices();
  time_decompose_ = real_t(clock() - t0)/real_t(CLOCKS_PER_SEC);

  vertex_.clear();
}

void
LaguerreCell<Polytope>::clip()
{
  // start clipping with the nearest neighbour and continue until the security radius is reached
  index_t j = 1;
  index_t nb_clip = 0;
  while (true)
  {
    if (j == neighbours_.size())
    {
      // max nearest neighbours
      enlarge_neighbours();
    }

    //index_t bj = neighbours_(site_,j);
    index_t bj = neighbours_[j];
    clip_by_bisector( j , bj );
    nb_clip++;

    if (security_radius_reached(bj)) break;
    j++;
    if (j == delaunay_.nb()) break;
    //if (polytope_.size() <= number_) break;
  }

  // add the points
  for (index_t j=0;j<polytope_.size();j++)
  {
    // save the geometry and bisector information
    points_.create( vertex_[polytope_[j]].X() );
    const std::vector<int>& bisectors = vertex_[ polytope_[j]].bisectors();
    points_.incidence().add( bisectors.data() , bisectors.size() );
    point2elem_.push_back(0);
    point2site_.push_back( site_ );
  }
  polytope_ = linspace( points_.nb() );
  add( polytope_.data() , polytope_.size() );
  cell2elem_.push_back(0);
  cell2site_.push_back(site_);
}

// SIMPLEX LAGUERRE CELLS
LaguerreCell<Simplex>::LaguerreCell(index_t elem , const Delaunay& delaunay ,
                                     GEO::NearestNeighborSearch& nns , const Topology<Simplex>& domain ,
                                     RVDFacets* facets,
                                     bool exact , index_t nb_nns ) :
  LaguerreCellBase(delaunay,nns,domain,exact,nb_nns),
  elem_(elem)
{
  set_facets(facets);

  vertex_.resize( domain_.nv(elem_) );
  for (index_t j=0;j<domain_.nv(elem_);j++)
  {
    vertex_[j] = Vertex(delaunay_.dim(),number_);
    vertex_[j].setCoordinates( domain_.points()[ domain_(elem_)[j] ] ,
                domain_.points().dim() );
    vertex_[j].setNumber( domain_.number() );
    vertex_[j].addSimplexVertex( domain_.points()[ domain_(elem_,j) ] );
  }

  // create the initial simplex
  std::vector<index_t> simplex_k = domain_.get(elem_);

  // add the bounding facets as bisectors
  index_t nf = domain_.nv(elem);
  for (index_t j=0;j<nf;j++)
  {
    std::vector<index_t> facet = simplex_k;
    facet.erase( facet.begin() +j );
    std::sort( facet.begin() , facet.end() );

    // get the facet label
    int b = domain_facets_->facet(facet);

    for (index_t i=0;i<simplex_k.size();i++)
    {
      if (i==j) continue; // skip the vertex opposite the facet
      vertex_[i].addBisector(b);
    }
  }
}


index_t
LaguerreCell<Simplex>::nb_sites() const
{
  index_t count = 0;
  for (index_t k=0;k<sites_.size();k++)
  {
    if (!clipped_[sites_[k]])
      count++;
  }
  return count;
}

void
LaguerreCell<Simplex>::next_site()
{
  index_t k;
  for (k=0;k<sites_.size();k++)
  {
    // skip sites that are alraedy clipped
    if (clipped_[sites_[k]]) continue;
    break; // found one!
  }
  clipped_[sites_[k]] = true;
  set_site(sites_[k]);
}

void
LaguerreCell<Simplex>::clip( index_t i )
{
  // create the initial points
  polytope_ = linspace(domain_.nv(elem_));

  // reset the edges to those of the initial simplex (in local indexing)
  pedges_.clear();
  for (index_t i = 0; i < polytope_.size(); i++)
  for (index_t j = i+1; j < polytope_.size(); j++)
  {
    pedges_.push_back(i);
    pedges_.push_back(j);
  }

  // start the clipping with the nearest neighbour, and continue until the security radius is reached
  index_t j = 1;
  while (true)
  {
    index_t bj = neighbours_[j];
    clip_by_bisector( j , bj );

    if (security_radius_reached(bj)) break;
    j++;
    if (j == delaunay_.nb()) break;
    if (j == neighbours_.size())
    {
      // max nearest neighbours
      enlarge_neighbours();
    }
  }

  // only polytopes of an admissible size are created
  if (polytope_.size()>number_)
  {
    // add the points
    index_t nb_points = points_.nb();
    for (index_t j=0;j<polytope_.size();j++)
    {
      // save the coordinate and bisector information
      points_.create( vertex_[polytope_[j]].X() );
      const std::vector<int>& bisectors = vertex_[ polytope_[j]].bisectors();
      points_.incidence().add( bisectors.data() , bisectors.size() );
      point2elem_.push_back(elem_);
      point2site_.push_back(site_);
    }
    polytope_ = linspace( polytope_.size() );
    for (index_t j=0;j<polytope_.size();j++)
      polytope_[j] += nb_points;
    add( polytope_.data() , polytope_.size() );
    cell2elem_.push_back(elem_);
    cell2site_.push_back(site_);
  }
}

void
LaguerreCell<Simplex>::clip()
{
  clock_t t0;

  // reset everything
  reset();
  clipped_.resize(delaunay_.nb(),false);

  // find the nearest delaunay site
  index_t s = delaunay_.closest( domain_.points()[domain_(elem_,0)] );
  sites_.push_back(s);
  set_site(s);

  while (true)
  {
    // re-initialize the neighbours
    std::vector<double> dist2(neighbours_.size(),0.0);
    t0 = clock();
    neighbour_search_.get_nearest_neighbors( neighbours_.size() , site_ , neighbours_.data() , dist2.data() );
    avro_assert( neighbours_[0] == site_ );
    time_neighbours_ += real_t(clock() -t0) / real_t(CLOCKS_PER_SEC);

    // clip with the next delaunay site
    clip(site_);

    // check if we're done, or move to the next delaunay site (voronoi cell)
    if (nb_sites() == 0) break;
    next_site();
  }
}

void
LaguerreCell<Simplex>::compute()
{
  // for measuring the initial time before doing something
  clock_t t0;

  // do the clipping
  t0 = clock();
  clip();
  time_clip_ += real_t(clock() - t0)/real_t(CLOCKS_PER_SEC);

  // option to decompose the simplices
  t0 = clock();
  if (number_ < 5) generate_simplices();
  time_decompose_ = real_t(clock() - t0)/real_t(CLOCKS_PER_SEC);
}

template<typename type>
int
LaguerreCellBase<type>::add_bisector( index_t p0 , index_t p1 )
{
  Bisector b(p0,p1);
  if (bisector_.find(b) == bisector_.end())
  {
    int id = bisector_.size();
    bisector_.insert( {b,id});
    ids_.insert({id,b});
  }
  return bisector_[b];
}

template<typename type>
void
LaguerreCellBase<type>::get_bisector( int b , index_t& p0 , index_t& p1 ) const
{
  std::map<int,Bisector>::const_iterator it;
  it = ids_.find(b);
  avro_assert_msg( it != ids_.end() , "bisector %d not found" , b );
  p0 = it->second.p0;
  p1 = it->second.p1;
}

template<typename type>
void
LaguerreCellBase<type>::clip_by_bisector( index_t j , index_t bj )
{
  // retrieve the bisector
  int b = add_bisector( site_ , bj );
  add_site(bj);

  // initialize the clipped polytope
  qpolytope_.clear();
  qedges_.clear();
  qplane_.clear();

  // retrieve the edges
  // this needs to be improved
  for (index_t i=0;i<pedges_.size()/2;i++)
  {
    // retrieve edge indices
    index_t e0 = pedges_[2*i];
    index_t e1 = pedges_[2*i+1];

    avro_assert( element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) );

    int q0 = -1,q1 = -1;
    int qs = clip_edge(e0,e1,b,qpolytope_,q0,q1);
    if (q0 < 0 || q1 < 0) continue;

    qedges_.push_back( q0 );
    qedges_.push_back( q1 );

    if (qs >= 0) qplane_.push_back(qs);
  }

  // go back through the new vertices and determine the edges
  for (index_t ii=0;ii<qplane_.size();ii++)
  for (index_t jj=ii+1;jj<qplane_.size();jj++)
  {
    index_t e0 = qplane_[ii];
    index_t e1 = qplane_[jj];
    if (element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) )
    {
      qedges_.push_back(e0);
      qedges_.push_back(e1);
    }
  }

  // the current polytope becomes the clipped one
  polytope_.assign( qpolytope_.begin() , qpolytope_.end() );
  pedges_.assign( qedges_.begin() , qedges_.end() );
  uniquify(polytope_);
  for (index_t ii=0;ii<polytope_.size();ii++)
    vertex_[ polytope_[ii] ].setBaseSite( delaunay_[site_] );
}

template<typename type>
int
LaguerreCellBase<type>::clip_edge( const index_t e0 , const index_t e1 , const int b , std::vector<index_t>& q , int& q0, int& q1 )
{
  int on_bisector = -1;
  Vertex& v0 = vertex_[e0];
  Vertex& v1 = vertex_[e1];

  index_t pi,pj;
  get_bisector( b , pi , pj );

  const real_t* zi = delaunay_[site_];
  const real_t* zj;
  if (pi!=site_)
  {
    avro_assert_msg(pj==site_,"bisector = %d, pj = %lu, site = %lu, int limits = (%d,%d)",b,pj,site_,std::numeric_limits<int>::min(),std::numeric_limits<int>::max());
    zj = delaunay_[pi];
  }
  else
  {
    zj = delaunay_[pj];
  }

  GEO::Sign side1 = v0.side( zi , zj , exact_ );
  GEO::Sign side2 = v1.side( zi , zj , exact_ );

  if (side1 == GEO::ZERO || side2 == GEO::ZERO)
  {
    return on_bisector;
  }

  avro_assert( side1!=GEO::ZERO && side2!=GEO::ZERO );

  if (side1!=side2)
  {
    Vertex v2(delaunay_.dim(),domain_.number());

    // perform the symbolic intersection of the bisectors, simplices, meshes
    v2.setBaseSite(delaunay_[site_]);
    v2.intersectSymbolic( &v0 , &v1 , delaunay_ );
    v2.setSites( delaunay_ , ids_ );

    // this vertex lies exactly on the slicing bisector
    v2.addBisector( b );
    v2.addSite( zj );

    // compute the geometric intersection
    v2.intersectGeometric( v0.X() , v1.X() , zi , zj );

    avro_assert_msg( v0.nb_bisectors()==domain_.number() ,
                "nb_bisectors = %lu, number = %u",v0.nb_bisectors(),number_ );

    // add the vertex after saving its location
    index_t vs = vertex_.size();
    vertex_.emplace_back(v2);

    if (side1==GEO::POSITIVE)
    {
      q.push_back( e0 );
      q0 = e0;
    }
    else
    {
      q.push_back( e1 );
      q0 = e1;
    }

    q.push_back(vs);
    q1 = vs;
    on_bisector = vs;

    avro_assert( v2.number() == domain_.number() );
    avro_assert( vertex_[vs].number() == domain_.number() );

  }
  else
  {
    // both points are on the same side
    if (side1==GEO::POSITIVE)
    {
      // both points are in the voronoi cell
      q.push_back( e0 );
      q.push_back( e1 );

      q0 = e0;
      q1 = e1;
    }
    else
    {
      // both points are outside the voronoi cell
      q0 = -1;
      q1 = -1;
    }
  }

  return on_bisector;
}

template<typename type>
bool
LaguerreCellBase<type>::security_radius_reached( index_t bj ) const
{
  const coord_t dim = delaunay_.dim();
  real_t R = -1.;
  real_t d;
  for (index_t k=0;k<polytope_.size();k++)
  {
    d = numerics::distance2( vertex_[polytope_[k]].X() , delaunay_[site_] , dim );
    if (d>R) R = d;
  }
  // a little more than 2^2 (4.1) because we need to make sure that we clip
  // against all possible bisectors (this just means we might try to clip
  // with a non-contributing cell which will terminate anyway)
  if (numerics::distance2( delaunay_[site_] , delaunay_[bj] , dim )>4.1*R)
  {
    // the security radius has been reached, signal we need to stop clipping
    return true;
  }
  return false;
}

template<typename type>
void
LaguerreCellBase<type>::generate_simplices()
{
  simplices_.clear();

  coord_t dim = delaunay_.dim();

  for (index_t j = 0; j < points_.nb(); j++)
    simplices_.add_point( points_[j] , point2elem_[j] , point2site_[j] );

  std::vector<real_t> xc(dim,0.);
  std::vector<int> facets;
  std::vector<index_t> simplex(number_+1);

  std::vector<index_t> vf;
  std::vector<int> edges;
  std::vector<index_t> ve;

  if (number_ == 2)
  {
    // add a point for the centroid of each polygon and triangulate
    for (index_t k=0;k<nb();k++)
    {
      index_t idc = simplices_.points().nb();
      numerics::centroid( (*this)(k) , nv(k) , points_ , xc );
      simplices_.add_point( xc.data() , cell2elem_[k] , cell2site_[k] );

      // get the hrep of this polygon
      element().hrep( (*this)(k) ,  nv(k) , facets );

      avro_assert_msg( facets.size() == nv(k) , "|facets| = %lu, nv(%lu) = %lu" , facets.size(),k,nv(k) ); // for polygons
      for (index_t j=0;j<facets.size();j++)
      {
        // get the points with this bisector
        vf.clear();
        element().vrep( (*this)(k) , nv(k) , facets[j] , vf );
        avro_assert( vf.size() == 2 ); // for polygons

        simplex[0] = vf[0];
        simplex[1] = vf[1];
        simplex[2] = idc;

        simplices_.add_simplex( simplex , cell2elem_[k] , cell2site_[k] );
      }
    }
  }
  else if (number_ == 3)
  {
    // add a point for the centroid of each polygon and triangulate
    std::vector<real_t> xf(dim,0.);
    for (index_t k=0;k<nb();k++)
    {
      index_t idc = simplices_.points().nb();
      numerics::centroid( (*this)(k) , nv(k) , points_ , xc );
      simplices_.add_point( xc.data() , cell2elem_[k] , cell2site_[k] );

      // get the hrep of this polyhedron
      element().hrep( (*this)(k) ,  nv(k) , facets );

      for (index_t j=0;j<facets.size();j++)
      {
        // get the points with this bisector
        vf.clear();
        element().vrep( (*this)(k) , nv(k) , facets[j] , vf );

        // compute the centroid of the facet
        index_t idf = simplices_.points().nb();
        numerics::centroid( vf.data() , vf.size() , points_ , xf );
        simplices_.add_point( xf.data() , cell2elem_[k] , cell2site_[k] );

        // retrieve the edges
        edges.clear();
        element().hrep( vf.data() , vf.size() , edges );
        for (index_t e=0;e<edges.size();e++)
        {
          if (edges[e] == facets[j]) continue;

          ve.clear();
          element().vrep( vf.data() , vf.size() , edges[e] , ve );
          if (ve.size() != 2) continue;
          avro_assert_msg( ve.size() == 2 , "|ve| = %lu" , ve.size() );

          simplex[0] = ve[0];
          simplex[1] = ve[1];
          simplex[2] = idf;
          simplex[3] = idc;
          simplices_.add_simplex( simplex , cell2elem_[k] , cell2site_[k] );
        }
      }
    }
  }
  else if (number_ == 4)
  {
    // add a point for centroid of each polyhedron
    std::vector<real_t> xp(dim,0.0);
    std::vector<real_t> xf(dim,0.0);
    std::vector<int> h_poly;
    std::vector<index_t> v_poly;
    std::vector<int> h_face;
    std::vector<index_t> v_face;
    std::vector<int> h_edge;
    std::vector<index_t> v_edge;

    for (index_t k=0;k<nb();k++)
    {
      // compute the centroid of this polytope
      index_t idc = simplices_.points().nb();
      numerics::centroid( (*this)(k) , nv(k) , points_ , xc );
      simplices_.add_point( xc.data() , cell2elem_[k] , cell2site_[k] );

      // get the hrep of this polytope
      h_poly.clear();
      element().hrep( (*this)(k) ,  nv(k) , h_poly );

      // loop through every polyhedron
      for (index_t j=0;j<h_poly.size();j++)
      {
        // get the points with this bisector
        v_poly.clear();
        element().vrep( (*this)(k) , nv(k) , h_poly[j] , v_poly );

        // compute the centroid of the polyhedron
        index_t idp = simplices_.points().nb();
        numerics::centroid( v_poly.data() , v_poly.size() , points_ , xp );
        simplices_.add_point( xp.data() , cell2elem_[k] , cell2site_[k] );

        // retrieve the face bisectors
        h_face.clear();
        element().hrep( v_poly.data() , v_poly.size() , h_face );
        for (index_t f = 0; f < h_face.size(); f++)
        {
          // skip if this is the current bisector of the polytope
          if (h_face[f] == h_poly[j]) continue;

          v_face.clear();
          element().vrep( v_poly.data() , v_poly.size() , h_face[f] , v_face );

          // compute the centroid of this face
          index_t idf = simplices_.points().nb();
          numerics::centroid( v_face.data() , v_face.size() , points_ , xf );
          simplices_.add_point( xf.data() , cell2elem_[k] , cell2site_[k] );

          // retrieve the edges
          h_edge.clear();
          element().hrep( v_face.data() , v_face.size() , h_edge );
          for (index_t e = 0; e < h_edge.size(); e++)
          {
            // skip if this is any of the bisectors we have already visited
            if (h_edge[e] == h_poly[j] || h_edge[e] == h_face[f]) continue;

            v_edge.clear();
            element().vrep( v_face.data() , v_face.size() , h_edge[e] , v_edge );
            if (v_edge.size() != 2) continue;

            simplex[0] = v_edge[0];
            simplex[1] = v_edge[1];
            simplex[2] = idf;
            simplex[3] = idp;
            simplex[4] = idc;
            simplices_.add_simplex( simplex , cell2elem_[k] , cell2site_[k] );
          }
        }
      }
    }
  }
  else
    avro_implement;

  // orient the simplices so they all have a positive volume
  simplices_.orient();
}

template<typename type>
LaguerreDiagram<type>::LaguerreDiagram( Delaunay& delaunay , const Topology<type>& domain ) :
  Topology<Polytope>(points_,domain.number()),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  domain_(domain),
  neighbour_search_(nullptr)
{}

template<>
void
LaguerreDiagram<Simplex>::initialize()
{
  // we need the facets but not the edges
  domain_facets_ = std::make_shared<RVDFacets>(domain_);
}

template<>
void
LaguerreDiagram<Simplex>::create( bool exact , index_t nb_nns )
{
  // create a cell for every simplex in the mesh
  cells_.resize( domain_.nb() , nullptr );
  for (index_t k=0;k<domain_.nb();k++)
  {
    cells_[k] = std::make_shared<LaguerreCell<Simplex>>(k,delaunay_,*neighbour_search_,domain_,domain_facets_.get(),exact,nb_nns);
    cells_[k]->set_facets( domain_facets_.get() );
    cells_[k]->set_edges( domain_edges_ );
  }
}

template<>
void
LaguerreDiagram<Polytope>::create( bool exact , index_t nb_nns )
{
  // create a cell for every delaunay vertex
  cells_.resize( delaunay_.nb() , nullptr );
  for (index_t k=0;k<delaunay_.nb();k++)
  {
    if (cells_[k] == nullptr || neighbour_search_ == nullptr)
      cells_[k] = std::make_shared<LaguerreCell<Polytope>>(k,delaunay_,*neighbour_search_,domain_,exact,nb_nns);
    cells_[k]->set_facets( domain_facets_.get() );
    cells_[k]->set_edges( domain_edges_ );
  }
}

template<>
void
LaguerreDiagram<Polytope>::initialize()
{
  // we need the edges but not the facets
  domain_edges_.clear();
  domain_.get_edges(domain_edges_);
  elem_.resize( delaunay_.nb() , 0 );
}

template<typename type>
void
LaguerreDiagram<type>::compute( bool exact , IntegrationSimplices* triangulation )
{
  // initialize the domain data
  initialize();

  // initialize the nearest neighbours
  real_t t0 = clock();
  index_t nb_nns = 50;
  if (delaunay_.nb() < nb_nns) nb_nns = delaunay_.nb();
  const coord_t dim = delaunay_.dim();
  if (neighbour_search_ == nullptr)
    neighbour_search_ = GEO::NearestNeighborSearch::create(dim,"ANN");
  std::vector<real_t> x(delaunay_.nb()*dim);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_[k][d];
  t0 = clock();
  neighbour_search_->set_points( delaunay_.nb() , x.data() );
  time_neighbours_ = real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // create the laguerre cells
  create(exact,nb_nns);

  // clip the cells
  #if 1
  ProcessCPU::parallel_for(
    parallel_for_member_callback( this , &thisclass::clip ),
    0,cells_.size() );
  #else
  for (index_t k=0;k<cells_.size();k++) {
    cells_[k]->compute();
  }
  #endif

  time_voronoi_   = 0;
  time_decompose_ = 0;

  // accumulate the result
  sites_.clear();
  for (index_t k=0;k<cells_.size();k++) {

    // retrieve the voronoi cell
    const LaguerreCell<type>& cell = *cells_[k].get();

    // accumulate the timing information
    time_neighbours_ += cell.time_neighbours() / ProcessCPU::maximum_concurrent_threads();
    time_voronoi_    += cell.time_clip() / ProcessCPU::maximum_concurrent_threads();
    time_decompose_  += cell.time_decompose() / ProcessCPU::maximum_concurrent_threads();

    // cells could be submerged for weighted sites
    if (cell.nb() == 0) continue;

    #if 1 // turn off for lower memory usage

    // add the cells
    for (index_t j=0;j<cell.nb();j++)
    {
      sites_.push_back( k );
      std::vector<index_t> polytope = cell.get(j);

      for (index_t i=0;i<polytope.size();i++)
        polytope[i] += points_.nb();
      add( polytope.data() , polytope.size() );
    }

    // add the points
    for (index_t j=0;j<cell.points().nb();j++)
    {
      points_.create( cell.points()[j] );
    }

    // add the vertex-to-facet information
    const Table<int>& vf = cell.element().incidence();
    for (index_t j=0;j<vf.nb();j++)
    {
      std::vector<int> f = vf.get(j);
      points_.incidence().add( f.data() , f.size() );
    }
    #endif

    // option to add the triangulation data
    if (triangulation == nullptr) continue;

    const IntegrationSimplices& tk = cell.simplices();
    std::vector<index_t> simplex(tk.number()+1);
    for (index_t j = 0; j < tk.nb(); j++)
    {
      for (index_t i = 0; i < tk.nv(j); i++)
        simplex[i] = tk(j,i) + triangulation->points().nb();
      triangulation->add_simplex( simplex , tk.simplex2elem(j) , tk.simplex2site(j) );
    }

    // this will create points that have the same dimension as the requested triangulation
    // in other words, even though the laguerre cell might be in dim+1, we only copy dim coordinates
    for (index_t j = 0; j < tk.points().nb(); j++)
      triangulation->add_point( tk.points()[j] , tk.point2elem(j) , tk.point2site(j) );
  }
  //printf("--> timing: neighbours (%3.4g sec.), clipping (%3.4g sec.), decompose (%3.4g sec.)\n",time_neighbours_,time_voronoi_,time_decompose);
  compute_neighbour_properties();
}

template<typename type>
void
LaguerreDiagram<type>::compute_neighbour_properties()
{
  neighbour_counts_.clear();
  neighbour_counts_.resize( delaunay_.nb() , 0 );

  std::vector<int> hrep;
  for (index_t k = 0; k < cells_.size(); k++)
  {
    for (index_t j = 0; j < cells_[k]->nb(); j++)
    {
      hrep.clear();
      cells_[k]->element().hrep( (*cells_[k])(j) , cells_[k]->nv(j) , hrep );

      index_t count = 0;
      index_t site = cells_[k]->cell2site(j);
      for (index_t i = 0; i < hrep.size(); i++)
      {
        //if (hrep[i] < 0) continue; // mesh facet

        // decode the bisector associated with this?
        count++;
      }
      neighbour_counts_[site] += count;
    }
  }

  index_t total_neighbours = 0;
  for (index_t k = 0; k < cells_.size(); k++)
    total_neighbours += neighbour_counts_[k];

  average_neighbours_ = real_t(total_neighbours) / real_t(cells_.size());
  minimum_neighbours_ = * std::min_element( neighbour_counts_.begin() , neighbour_counts_.end() );
  maximum_neighbours_ = * std::max_element( neighbour_counts_.begin() , neighbour_counts_.end() );
}

template<typename type>
LaguerreDiagram<type>::~LaguerreDiagram() {
  if (neighbour_search_ != nullptr) delete neighbour_search_;
}

real_t
Integrand_Transport_Energy::operator()( index_t k , const QuadraturePoint& point , const real_t* x ) const
{
  avro_assert_msg( k < simplices_.nb() , "elem = %lu, nb_simplices = %lu", k , simplices_.nb() );

  // retrieve the delaunay site and coordinate
  index_t site = simplices_.simplex2site(k);
  avro_assert_msg( site < delaunay_.nb() , "elem = %lu, requested site %lu but nb_delaunay = %lu",k,site,delaunay_.nb());
  const real_t* z = delaunay_[site];

  // evaluate the density
  real_t rho = density_.evaluate( k , point.coordinate() , x );

  T f = 0;
  for (coord_t d=0;d<dim_;d++)
    f += rho*(z[d] - x[d])*(z[d] - x[d]);
  return f;
}

class Integrand_Mass : public Integrand<Integrand_Mass>
{
public:

  typedef real_t T;

  Integrand_Mass(const DensityMeasure& density) :
    density_(density)
  {}

  bool needs_solution() const { return false; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( index_t k , const QuadraturePoint& point , const real_t* x ) const
  {
    return density_.evaluate(k,point.coordinate(),x);
  }

private:
  const DensityMeasure& density_;
};

class Integrand_Moment : public Integrand<Integrand_Moment>
{
public:

  typedef real_t T;

  Integrand_Moment(coord_t dim, const DensityMeasure& density) :
    dim_(dim),
    density_(density)
  {}

  bool needs_solution() const { return false; }
  bool needs_gradient() const { return false; }
  bool needs_hessian() const { return false; }

  T operator()( index_t k , const QuadraturePoint& point , const real_t* x , std::vector<T>& I ) const
  {
    real_t rho = density_.evaluate(k,point.coordinate(),x);
    for (index_t r = 0; r < dim_; r++)
      I[r] = rho*x[r];
    return 1.0;
  }

private:
  coord_t dim_;
  const DensityMeasure& density_;
};

template<typename type>
SemiDiscreteOptimalTransport<type>::SemiDiscreteOptimalTransport( const Topology<type>& domain , DensityMeasure* density ) :
  domain_(domain),
  density_(density),
  delaunay_(domain.points().dim()),
  diagram_(delaunay_,domain_),
  simplices_(domain.number(),domain.number()),
  exact_(false),
  print_(true),
  quad_order_(4),
  weight_max_(1e6),
  save_every_(1e6),
  prefix_("avro-sdot")
{}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::sample( index_t nb_samples )
{
  coord_t dim = domain_.points().dim();

  // create samples in a dim+1 space, but leave the last dimension as 0
  delaunay_.clear();
  for (index_t k=0;k<nb_samples;k++)
  {
    std::vector<real_t> p(dim,0.);
    for (coord_t d = 0; d < domain_.number(); d++)
      p[d] = random_within(0.0,1.0);
    delaunay_.create(p.data());
  }
  nu_.resize( delaunay_.nb() , 0.0 );
  weight_.resize( delaunay_.nb() , 0.0 );
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::compute_laguerre()
{
  // compute the laguerre diagram along with the integration simplices
  diagram_.clear();
  diagram_.points().clear();
  simplices_.clear();
  diagram_.compute(exact_,&simplices_);
}

real_t
calculate_norm( const real_t* x , index_t nb )
{
  real_t n = 0.0;
  for (index_t k = 0; k < nb; k++)
    n += x[k]*x[k];
  return std::sqrt(n);
}

template<typename type>
real_t
SemiDiscreteOptimalTransport<type>::evaluate( real_t* dc_dx , real_t* dc_dw )
{
  avro_assert( density_ != nullptr );
  clock_t t0;
  real_t time_integrate = 0.0;

  coord_t dim = domain_.number();

  // compute the laguerre diagram
  compute_laguerre();

  // define the integration rules
  simplices_.element().set_basis( BasisFunctionCategory_Lagrange );

  // get the masses
  t0 = clock();
  typedef delaunay::Integrand_Mass Integrand_Mass_t;
  Integrand_Mass_t integrandm(*density_);
  Functional<Integrand_Mass_t> fm(integrandm,simplices_);
  std::vector<real_t> masses(simplices_.nb());
  fm.integrate(simplices_,masses.data());
  time_integrate += real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // get the moments
  t0 = clock();
  typedef delaunay::Integrand_Moment Integrand_Moment_t;
  Integrand_Moment_t integrandxm(dim,*density_);
  Functional_Ranked<Integrand_Moment_t> fxm(integrandxm,dim,simplices_);
  std::vector<real_t> moments( dim * simplices_.nb() );
  fxm.integrate(simplices_,moments.data());
  time_integrate += real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // first compute the mass of each voronoi cell
  mass_.clear();
  mass_.resize( delaunay_.nb() , 0.0 );
  for (index_t k = 0; k < simplices_.nb(); k++)
  {
    index_t s = simplices_.simplex2site(k);
    mass_[s] += masses[k];
  }

  if (dc_dx != nullptr)
  {
    // compute the and center of mass of each voronoi cell
    centroid_.clear();
    centroid_.resize( delaunay_.nb()*dim , 0.0 );
    for (index_t k = 0; k < simplices_.nb(); k++)
    {
      index_t s  = simplices_.simplex2site(k);
      for (coord_t d = 0; d < dim; d++)
        centroid_[ s*dim + d ] += moments[ k*dim + d];
    }

    // compute the centroid and gradient
    for (index_t k = 0; k < delaunay_.nb(); k++)
    for (coord_t d = 0; d < dim; d++)
    {
      // if weights are too high, some cells could be submerged which gives a zero volume
      if (mass_[k] < 1e-12) centroid_[k*dim+d] = delaunay_[k][d];
      else
        centroid_[k*dim+d] /= mass_[k];
      dc_dx[ k*dim + d] = 2.0*mass_[k]*( delaunay_[k][d] - centroid_[k*dim+d] );
    }
  }

  if (dc_dw != nullptr)
  {
    // compute the gradient of the transport map
    for (index_t k = 0; k < delaunay_.nb(); k++)
      dc_dw[k] = nu_[k] - mass_[k];
  }

  // get the energy
  t0 = clock();
  typedef delaunay::Integrand_Transport_Energy Integrand_t;
  Integrand_t integrand(delaunay_,simplices_,*density_,dim);
  Functional<Integrand_t> f(integrand,simplices_);
  f.integrate(simplices_);
  real_t energy = f.value();
  time_integrate += real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  // add te contribution from the weights and dirac masses
  for (index_t k=0;k<delaunay_.nb();k++)
    energy += weight_[k]*(nu_[k] - mass_[k]);

  // calculate the norm of the gradient and report convergence criteria
  real_t gnorm;
  if (mode_ == 0)
    gnorm = calculate_norm(dc_dx,delaunay_.nb()*dim);
  else
    gnorm = calculate_norm(dc_dw,delaunay_.nb());


  time_integrate /= ProcessCPU::maximum_concurrent_threads();
  real_t mass_min = * std::min_element( mass_.begin() , mass_.end() );
  real_t mass_max = * std::max_element( mass_.begin() , mass_.end() );

  if (print_)
  {
    printf("%6lu|%10.3e|%10.3e|%11.1e|%11.1e|%11.1e|%11.1e|%11.3e|%11.3e|%6lu|%6lu|%11.1e\n",
      iteration_,
      energy,
      gnorm,
      diagram_.time_voronoi(),
      diagram_.time_neighbours(),
      diagram_.time_decompose(),
      time_integrate,
      mass_min,
      mass_max,
      diagram_.minimum_neighbours(),
      diagram_.maximum_neighbours(),
      diagram_.average_neighbours()
    );
  }

  properties_.mass_min.push_back( mass_min );
  properties_.mass_max.push_back( mass_max );
  properties_.time_voronoi.push_back( diagram_.time_voronoi() );
  properties_.time_neighbours.push_back( diagram_.time_neighbours() );
  properties_.time_integration.push_back( time_integrate );
  properties_.time_triangulation.push_back( diagram_.time_decompose() );
  properties_.neighbours_average.push_back( diagram_.average_neighbours() );
  properties_.neighbours_minimum.push_back( diagram_.minimum_neighbours() );
  properties_.neighbours_maximum.push_back( diagram_.maximum_neighbours() );
  properties_.energy.push_back( energy );
  properties_.gradient.push_back( gnorm );
  time_integrate_ = time_integrate;

  if (iteration_ % save_every_ == 0) save_snapshot();

  return energy;
}

void
SDOT_Properties::save( const std::string& filename ) const
{
  json J;
  J["mass_min"] = mass_min;
  J["mass_max"] = mass_max;
  J["time_voronoi"] = time_voronoi;
  J["time_neighbours"] = time_neighbours;
  J["integration"] = time_integration;
  J["triangulation"] = time_integration;
  J["neighbours_average"] = neighbours_average;
  J["neighbours_minimum"] = neighbours_minimum;
  J["neighbours_maximum"] = neighbours_maximum;
  J["energy"] = energy;
  J["gradient"] = gradient;
  std::ofstream output(filename);
  output << std::setw(4) << J << std::endl;
  output.close();
}

void
SDOT_Snapshot::save( const std::string& filename ) const
{
  json J;
  J["dim"] = dim;
  J["points"] = points;
  J["weights"] = weights;
  J["neighbours"] = neighbours;
  J["density"] = name;
  J["mass"] = mass;
  std::ofstream output(filename);
  output << std::setw(4) << J << std::endl;
  output.close();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::save_snapshot()
{
  snapshot_.clear();
  std::vector<real_t> coordinates(delaunay_.nb() * domain_.number());
  for (index_t k = 0, j = 0 ; k < delaunay_.nb(); k++)
  for (index_t d = 0; d < domain_.number(); d++)
    coordinates[j++] = delaunay_[k][d];
  snapshot_.dim = domain_.number();
  snapshot_.points = coordinates;
  snapshot_.weights = weight_;
  snapshot_.neighbours = diagram_.neighbour_counts();
  snapshot_.mass = mass_;
  snapshot_.name = diagram_.name();
  snapshot_.save(prefix_ + "-snapshot-iter-" + std::to_string(iteration_) + ".json" );
}

template<typename type>
struct nlopt_data_optimal_transport
{
	SemiDiscreteOptimalTransport<type>& transport;
	index_t eval_count;
	real_t objective;
  int mode; // 0 for coordinates , 1 for weights
};

template<typename type>
void
SemiDiscreteOptimalTransport<type>::set_delaunay( const real_t* x , coord_t dim )
{
  for (index_t k=0;k<delaunay_.nb();k++)
  for (coord_t d=0;d<dim;d++)
    delaunay_[k][d] = x[k*dim+d];
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::set_weights( const real_t* w )
{
  for (index_t k=0;k<delaunay_.nb();k++)
    weight_[k] = w[k];
  coord_t dim = delaunay_.dim() -1;
  real_t wm = * std::max_element( weight_.begin() , weight_.end() );
  //printf("setting dim = %u, wm = %g\n",dim,wm);
  for (index_t k=0;k<weight_.size();k++)
    delaunay_[k][dim] = std::sqrt( std::max( wm - weight_[k] , 0.0 ) );
}

template<typename type>
real_t
SemiDiscreteOptimalTransport<type>::transport_objective( index_t n , const double* x , double* grad )
{
  // set the coordinates of the delaunay points
  if (mode_ == 0)
    set_delaunay(x,domain_.number());
  else if (mode_ == 1)
    set_weights(x);
  else
    avro_assert_not_reached;

  // recompute the voronoi diagram
  compute_laguerre();

  index_t nb_points = delaunay_.nb();
  coord_t dim = domain_.number();

  // compute the energy and gradients
  std::vector<real_t> dE_dZ(nb_points*dim),dE_dW(nb_points);
  real_t energy = evaluate( dE_dZ.data() , dE_dW.data() );
  if (grad != nullptr)
  {
    if (mode_ == 0)
    {
      // cvt mode
      for (index_t k=0;k<dE_dZ.size();k++)
        grad[k] = dE_dZ[k];
    }
    else if (mode_ == 1)
    {
      // otm mode: flip the sign to make it convex
      energy *= -1.0;
      for (index_t k=0;k<dE_dW.size();k++)
        grad[k] = -dE_dW[k];
    }
    else
      avro_assert_not_reached;
  }
  return energy;
}

template<typename type>
double
nlopt_transport_objective( unsigned n , const double* x , double* grad, void* data0 )
{
  nlopt_data_optimal_transport<type>* data = (nlopt_data_optimal_transport<type>*)(data0);
  SemiDiscreteOptimalTransport<type>& transport = data->transport;
  transport.iteration()++;
  return transport.transport_objective( n , x , grad );
}

#define OPTIMIZER_NLOPT 1
#define OPTIMIZER_HLBFGS 0
OptimalTransportBase* __transport__ = nullptr;

#if OPTIMIZER_HLBFGS
static void
newiteration(int iter, int call_iter, double *x, double* f, double *g,  double* gnorm)
{}

void
hlbfgs_transport_objective(int N, double* x, double *prev_x, double* f, double* g)
{
  avro_assert( __transport__ != nullptr );
  OptimalTransportBase* transport = static_cast<OptimalTransportBase*>(__transport__);
  transport->iteration()++;
  *f = transport->transport_objective(index_t(N),x,g);
}
#endif

template<typename type>
void
SemiDiscreteOptimalTransport<type>::start()
{
  printf("--------------------------------------------------------------------------------------\n");
  printf("%5s|%10s|%10s|%10s|%10s|%10s|%10s|%10s\n","iter  "," energy "," gradient "," t_vor (s) "," t_knn (s) " , " t_tri (s) " , " t_int (s) ", " min(m) ");
  printf("--------------------------------------------------------------------------------------\n");
  iteration_ = 0;
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::optimize_points( index_t nb_iter )
{
  properties_.clear();

  mode_ = 0;
  coord_t dim = domain_.number();
  index_t n = delaunay_.nb()*dim;

  std::vector<real_t> x(n);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_(k,d);

  // setup the optimizer
  #if 1 // OPTIMIZER_NLOPT
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_optimal_transport<type> data = {*this,0,1,0};

	// set the objective function
	opt.set_min_objective( &nlopt_transport_objective<type> , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-12);
  opt.set_ftol_rel(1e-12);
  opt.set_maxeval(nb_iter);

  // set the lower and upper bounds on the weights
  std::vector<real_t> lower_bound( n , 0.0 );
  std::vector<real_t> upper_bound( n ,  1.0 );
  opt.set_lower_bounds(lower_bound);
  opt.set_upper_bounds(upper_bound);

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  start();
  result = opt.optimize(x,f_opt);
  std::string desc = nloptResultDescription(result);
  printf("nlopt result: %s\n",desc.c_str());

  #else

  double parameter[20];
  int info[20];
  int T = 0;
  int M = 7;
  bool with_hessian = false;

  // initialize
  INIT_HLBFGS(parameter, info);
  parameter[5] = 1e-16;
  info[3] = 1;
  info[4] = nb_iter;
  info[6] = T;
  info[5] = 1; // verbose = 1
  info[7] = with_hessian ? 1 : 0;
  info[10] = 0;
  info[11] = 1;

  __transport__ = this;
  start();
  HLBFGS(n, M, x.data() , hlbfgs_transport_objective , 0 , HLBFGS_UPDATE_Hessian, newiteration, parameter, info);

  #endif

  save_snapshot();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::optimize_points_lloyd( index_t nb_iter )
{
  properties_.clear();

  mode_ = 0;
  coord_t dim = domain_.number();
  index_t n = delaunay_.nb()*dim;

  centroid_.clear();
  centroid_.resize(n);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    centroid_[k*dim+d] = delaunay_(k,d);

  start();
  for (iteration_ = 0; iteration_ < nb_iter ; iteration_++)
  {
    // set the current delaunay vertices
    set_delaunay( centroid_.data() , dim );

    // evaluate the voronoi diagram for these points
    std::vector<real_t> dc_dx(n); // not used
    evaluate( dc_dx.data() , nullptr );
  }

  save_snapshot();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::optimize_weights( index_t nb_iter )
{
  properties_.clear();

  mode_ = 1;
  index_t n = delaunay_.nb();
  std::vector<real_t> x(n);
  for (index_t k=0;k<delaunay_.nb();k++)
    x[k] = 0.0;//mass_[k];//weight_[k];

  #if OPTIMIZER_NLOPT

  // setup the optimizer
  nlopt::opt opt( nlopt::LD_LBFGS , n );
  nlopt_data_optimal_transport<type> data = {*this,0,1,1};

	// set the objective function
	opt.set_min_objective( &nlopt_transport_objective<type> , static_cast<void*>(&data) );

  // set some optimization parameters
  opt.set_xtol_rel(1e-16);
  opt.set_ftol_rel(1e-16);
  opt.set_maxeval(nb_iter);

  // set the lower and upper bounds on the weights
  std::vector<real_t> lower_bound( n , 0.0 );
  std::vector<real_t> upper_bound( n ,  weight_max_ );
  opt.set_lower_bounds(lower_bound);
  opt.set_upper_bounds(upper_bound);

  real_t f_opt;
  nlopt::result result = nlopt::result::SUCCESS;
  start();
  result = opt.optimize(x,f_opt);
  std::string desc = nloptResultDescription(result);
  printf("nlopt result: %s\n",desc.c_str());

  #elif OPTIMIZER_HLBFGS

  double parameter[20];
  int info[20];
  int T = 0;
  int M = 7;
  bool with_hessian = false;

  // initialize
  INIT_HLBFGS(parameter, info);
  parameter[5] = 1e-12;
  info[3] = 1;
  info[4] = nb_iter;
  info[6] = T;
  info[5] = 1; // verbose = 1
  info[7] = with_hessian ? 1 : 0;
  info[10] = 0;
  info[11] = 1;

  __transport__ = this;
  start();
  HLBFGS(n, M, x.data() , hlbfgs_transport_objective , 0 , HLBFGS_UPDATE_Hessian, newiteration, parameter, info);

  #elif OPTIMIZER_LIBLBFGS

  lbfgsfloatval_t fx;
  lbfgsfloatval_t *x0 = lbfgs_malloc(n);
  lbfgs_parameter_t param;

  for (index_t i = 0; i < n; i++)
    x0[i] = x[i];

  // initialize
  lbfgs_parameter_init(&param);
  param.linesearch = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
  param.max_iterations = nb_iter;

  // optimize!
  __transport__ = this;
  start();
  int ret = lbfgs(n, x0, &fx, liblbfgs_evaluate, liblbfgs_progress, NULL, &param);

  // report the result
  printf("L-BFGS optimization terminated with status code = %d\n", ret);
  lbfgs_free(x0);

  #else

  start();
  std::vector<real_t> gradient( n , 0. );
  std::vector<real_t> gradient0( n , 0. );
  std::vector<real_t> direction( n , 0 );
  bool bfgs = false;

  for (iteration_ = 0; iteration_ < nb_iter; iteration_++)
  {
    // compute the gradient
    std::vector<real_t> dc_dw(n);
    real_t f0 = transport_objective( n , x.data() , dc_dw.data() );

    // determine the descent direction
    for (index_t i = 0; i < n; i++)
    {
      gradient[i] = dc_dw[i];
      direction[i] = -gradient[i];
    }

    // option for bfgs udate
    if (bfgs)
    {
      avro_implement;
    }
    real_t mass_min0 = * std::min_element( mass_.begin() , mass_.end() );

    real_t alpha = 1.0;
    std::vector<real_t> x0( x.begin() , x.end() );
    real_t gnorm0 = calculate_norm( gradient.data() , n );
    real_t c1 = 1e-4;

    // perform a line search on the optimization direction
    index_t nb_linesearch = 0;
    while (true)
    {
      // update the weights in the descent direction
      for (index_t i = 0; i < n; i++)
        x[i] = x0[i] + alpha*direction[i];

      print_ = false;
      real_t f = transport_objective( n , x.data() , dc_dw.data() );
      real_t gnorm = calculate_norm( dc_dw.data() , n );
      print_ = true;

      real_t pg = 0.0;
      for (index_t i = 0; i < n; i++)
        pg += direction[i]*dc_dw[i];

      // determine the minimum mass
      real_t mass_min = * std::min_element( mass_.begin() , mass_.end() );
      if (mass_min > 1e-12 && f < f0 + c1*alpha*pg) break;
      //if (mass_min > 1e-12 && gnorm <= (1-.5*alpha)*gnorm0) break;
      alpha *= 0.5;

      //printf("\t -> linesearch %lu, alpha = %g, gnorm = %g, gnorm0 = %g, mass_min = %g\n",nb_linesearch,alpha,gnorm,gnorm0,mass_min);

      nb_linesearch++;
    }
  }

  #endif
  save_snapshot();
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::generate_bluenoise()
{
  // optimize the points
  optimize_points(100);

  // compute the target mass
  real_t m_total = 0.0;
  for (index_t k = 0; k < mass_.size(); k++)
    m_total += mass_[k];
  nu_.clear();
  nu_.resize( delaunay_.nb() , m_total/delaunay_.nb() );

  printf("total mass = %g, target mass = %g\n",m_total,m_total/delaunay_.nb());

  // generate blue noise
  for (index_t iter = 0; iter < 5; iter++)
  {
    // compute the transport map
    //optimize_weights(10);

    // optimize the points
    //optimize_points_lloyd(1);
  }

  // generate blue noise
  for (index_t iter = 0; iter < 5; iter++)
  {
    // compute the transport map
    //optimize_weights(10);

    // optimize the points
  //  optimize_points(2);
  }
}

template<typename type>
void
SemiDiscreteOptimalTransport<type>::stochastic_gradient_descent( index_t nb_iter )
{

  start();
  mode_ = 0;

  index_t n = delaunay_.nb();
  coord_t dim = domain_.number();

  // initialize the gradient
  std::vector<real_t> f(n, 0.0);
  std::vector<real_t> dE_dw( n ,0.0 );
  std::vector<real_t> a( n , 1.0 );
  for (iteration_ = 0; iteration_ < nb_iter; iteration_++)
  {
    // sample the density measure
    std::vector<real_t> y(dim,0.0);

    // determine which laguerre cell contains y
    index_t i = 0;

    // update the gradient
    dE_dw = a;
    dE_dw[i] -= 1.0;

    real_t l0 = 10.0;
    real_t tau = 0.1/( 1.0 + iteration_/l0 );

    f[i] += tau*dE_dw[i];

    // set the weights and re-evaluate
    set_weights( a.data() );
    evaluate( nullptr , nullptr );
  }

}

template<typename type>
std::vector<real_t>
SemiDiscreteOptimalTransport<type>::get_sites() const {
  index_t nb_points = delaunay_.nb();
  index_t dim = domain_.number(); // maybe make dim a parameter to this function
  std::vector<real_t> sites( nb_points*dim );
  index_t i = 0;
  for (index_t k = 0; k < nb_points; k++)
  for (coord_t d = 0; d < dim; d++)
    sites[i++] = delaunay_[k][d];
  return sites;
}


template class LaguerreCellBase<Simplex>;
template class LaguerreCellBase<Polytope>;
template class LaguerreDiagram<Simplex>;
template class LaguerreDiagram<Polytope>;
template class SemiDiscreteOptimalTransport<Simplex>;
template class SemiDiscreteOptimalTransport<Polytope>;

} // delaunay

} // avro
