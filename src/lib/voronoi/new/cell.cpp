#include "mesh/points.h"
#include "numerics/geometry.h"
#include "voronoi/new/cell.h"

#include <time.h>

namespace avro
{

namespace voronoi
{

Cell::Cell( index_t site , const Points& delaunay ,
      const Topology<Simplex>& domain , GEO::NearestNeighborSearch& searcher ) :
  Topology<Polytope>(points_,domain.number()),
  site_(site),
  delaunay_(delaunay),
  domain_(domain),
  search_(searcher),
  points_(delaunay.dim()),
  exact_(true)
{}

void
Cell::compute( index_t elem ) {

  // do the clipping
  clip(elem);

  // option to decompose the simplices
  //if (number_ < 5) generate_simplices();
}

void
Cell::clip( index_t elem ) {

  if (delaunay_.nb() < 10)
    neighbours_.resize(delaunay_.nb());
  else
    neighbours_.resize(10);
  std::vector<double> dist2(neighbours_.size(),0.0);
  search_.get_nearest_neighbors( neighbours_.size() , site_ , neighbours_.data() , dist2.data() );

  // clip by the closest simplex
  // this will march into simplex neighbours until no more clipping is needed
  clip_simplex(elem);
}

void
Cell::clip_simplex( index_t elem ) {

  // check if this element has already been visited
  if (visited_.find(elem) != visited_.end()) return;

  // initialize the edges and bisectors
  visited_.insert( elem );
  initialize( elem );

  // start clipping with the nearest neighbour and continue until the security radius is reached
  index_t j = 1; // first nearest neighbour
  index_t nb_clip = 0;
  while (true) {

    if (j == neighbours_.size()) {
      // we reached the last nearest neighbour, so we need to enlarge them
      enlarge_neighbours();
    }

    index_t bj = neighbours_[j];
    clip_by_bisector( j , bj );
    nb_clip++;

    if (security_radius_reached(bj)) break;
    j++;
    if (j == delaunay_.nb()) break;
  }

  if (polytope_.size() == 0) return;

  // add the points
  index_t nb_points = points_.nb();
  for (index_t j = 0; j < polytope_.size(); j++) {
    // save the coordinates and bisector information
    points_.create( vertex_[polytope_[j]].X() );
    const std::vector<int>& bisectors = vertex_[ polytope_[j]].bisectors();
    points_.incidence().add( bisectors.data() , bisectors.size() );
    polytope_[j] = nb_points + j;
  }

  // add the polytope
  add( polytope_.data() , polytope_.size() );

  return;
  // move to the neighbours of the current element
  for (index_t j = 0; j < domain_.number()+1; j++) {
    int n = domain_.neighbours()(elem,j);
    if (n < 0) continue; // do not step into boundaries
    clip_simplex( index_t(n) );
  }
}

void
Cell::initialize( index_t elem ) {

  const coord_t dim = delaunay_.dim();
  const coord_t n = domain_.number();

  // get the list of facets of this element
  std::vector<int> facets(n+1);
  for (index_t j = 0; j < n+1; j++)
    facets[j] = -j - 1; // TODO use the mesh facets which have global numbering (may be needed for boundary conditions)

  // create the initial points
  vertex_.clear();
  vertex_.resize( domain_.nv(elem) );
  for (index_t j = 0; j < domain_.nv(elem); j++) {

    vertex_[j] = Vertex(dim,n);
    vertex_[j].set_coordinates( domain_.points()[ domain_(elem)[j] ] , dim );
    vertex_[j].set_number(n);
    vertex_[j].add_simplex_vertex( domain_.points()[ domain_(elem,j) ] );

    // get the mesh facet this vertex is opposite
    for (index_t i = 0; i < facets.size(); i++) {
      if (i == j) continue;
      vertex_[j].add_bisector(facets[i]);
    }
  }

  printf("--> intersecting with elem %lu\n",elem);

  // initialize the polytope
  polytope_ = linspace(domain_.nv(elem));

  // initialize the polytope edges
  pedges_.clear();
  #if 0
  for (index_t ii = 0; ii < polytope_.size(); ii++)
  for (index_t jj = ii+1; jj < polytope_.size(); jj++) {
    index_t e0 = polytope_[ii];
    index_t e1 = polytope_[jj];
    if (element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) ) {
      pedges_.push_back(e0);
      pedges_.push_back(e1);
    }
  }

  #else
  for (index_t ii = 0; ii < n+1; ii++)
  for (index_t jj = ii+1; jj < n+1; jj++) {
    pedges_.push_back(ii);
    pedges_.push_back(jj);
  }
  #endif
}

void
Cell::enlarge_neighbours() {

  index_t nb_nns = neighbours_.size();
  nb_nns += 10;
  if (nb_nns > delaunay_.nb()) nb_nns = delaunay_.nb();
  std::vector<real_t> dist2(nb_nns); // not used
  neighbours_.resize( nb_nns );
  search_.get_nearest_neighbors( nb_nns , delaunay_[site_] , neighbours_.data() , dist2.data()  );
}

void
Cell::generate_simplices() {
  avro_implement;
}

int
Cell::add_bisector( index_t p0 , index_t p1 ) {
  Bisector b(p0,p1);
  if (bisector_.find(b) == bisector_.end()) {
    int id = bisector_.size();
    bisector_.insert( {b,id});
    ids_.insert({id,b});
  }
  return bisector_[b];
}

void
Cell::get_bisector( int b , index_t& p0 , index_t& p1 ) const {
  std::map<int,Bisector>::const_iterator it;
  it = ids_.find(b);
  avro_assert_msg( it != ids_.end() , "bisector %d not found" , b );
  p0 = it->second.p0;
  p1 = it->second.p1;
}

void
Cell::clip_by_bisector( index_t j , index_t bj ) {

  // retrieve the bisector
  int b = add_bisector( site_ , bj );

  printf("  ==> clipping by bisector %d\n",b);
  delaunay_.print(  bj );

  // initialize the clipped polytope
  qpolytope_.clear();
  qedges_.clear();
  qplane_.clear();

  // retrieve the edges
  std::set< std::pair<index_t,index_t> > E;
  for (index_t i = 0; i < pedges_.size()/2; i++){

    // retrieve edge indices
    index_t e0 = pedges_[2*i];
    index_t e1 = pedges_[2*i+1];

    // can probably turn this off since it adds computation
    avro_assert( element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) );

    printf("\tclipping with edge (%lu,%lu)\n",e0,e1);

    // clip the edge by the bisector
    int q0 = -1, q1 = -1;
    int qs = clip_edge(e0,e1,b,qpolytope_,q0,q1);
    if (q0 < 0 || q1 < 0) continue;

    if (qs < 0) {
      printf("\t --> edge (%lu,%lu) is on the same side\n",e0,e1);
    }
    else {
      printf("\t--> edge( %lu,%lu) is intersected!\n",e0,e1);
      vertex_[e0].print("\t--> e0",true);
      vertex_[e1].print("\t--> e1",true);
      vertex_[qs].print("\t--> intersection",true);
    }

    // add the vertices to the new polytope
    qedges_.push_back( q0 );
    qedges_.push_back( q1 );

    E.insert( {q0,q1} );

    // add the vertex to the list on the bisector
    if (qs >= 0) qplane_.push_back(qs);
  }

  // go back through the new vertices on the bisector and determine the edges
  for (index_t ii = 0; ii < qplane_.size(); ii++)
  for (index_t jj = ii+1; jj< qplane_.size(); jj++) {
    index_t e0 = qplane_[ii];
    index_t e1 = qplane_[jj];
    if (element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) ) {
      qedges_.push_back(e0);
      qedges_.push_back(e1);
      avro_assert( E.find( {e0,e1} ) == E.end() );
    }
  }

  // the current polytope becomes the clipped one
  polytope_.assign( qpolytope_.begin() , qpolytope_.end() );
  pedges_.assign( qedges_.begin() , qedges_.end() );
  uniquify(polytope_);
  for (index_t ii = 0; ii < polytope_.size(); ii++)
    vertex_[ polytope_[ii] ].set_base_site( delaunay_[site_] );
}

int
Cell::clip_edge( const index_t e0 , const index_t e1 , const int b , std::vector<index_t>& q , int& q0, int& q1 ) {

  int on_bisector = -1;
  Vertex& v0 = vertex_[e0];
  Vertex& v1 = vertex_[e1];

  index_t pi,pj;
  get_bisector( b , pi , pj );

  const real_t* zi = delaunay_[site_];
  const real_t* zj;
  if (pi != site_) {
    avro_assert_msg( pj == site_, "bisector = %d, pj = %lu, site = %lu, int limits = (%d,%d)",
                    b,pj,site_,std::numeric_limits<int>::min(),std::numeric_limits<int>::max());
    zj = delaunay_[pi];
  }
  else {
    zj = delaunay_[pj];
  }

  GEO::Sign side1 = v0.side( zi , zj , exact_ );
  GEO::Sign side2 = v1.side( zi , zj , exact_ );

  if (side1 == GEO::ZERO || side2 == GEO::ZERO) {
    avro_assert_not_reached;
    return on_bisector;
  }
  avro_assert( side1 != GEO::ZERO && side2 != GEO::ZERO );

  if (side1 != side2) {

    Vertex v2(delaunay_.dim(),domain_.number());

    // perform the symbolic intersection of the bisectors, simplices, meshes
    v2.set_base_site(delaunay_[site_]);
    v2.intersect_symbolic( &v0 , &v1 );
    v2.set_sites( delaunay_ , ids_ );

    // this vertex lies exactly on the slicing bisector
    v2.add_bisector( b );
    v2.add_site( zj );

    // compute the geometric intersection
    v2.intersect_geometric( v0.X() , v1.X() , zi , zj );

    avro_assert_msg( v0.nb_bisectors()==domain_.number() ,
                "nb_bisectors = %lu, number = %u",v0.nb_bisectors(),number_ );

    // add the vertex after saving its location
    index_t vs = vertex_.size();
    vertex_.emplace_back(v2);

    if (side1 == GEO::POSITIVE) {
      q.push_back( e0 );
      q0 = e0;
    }
    else {
      q.push_back( e1 );
      q0 = e1;
    }

    q.push_back(vs);
    q1 = vs;
    on_bisector = vs;

    avro_assert( v2.number() == domain_.number() );
    avro_assert( vertex_[vs].number() == domain_.number() );
  }
  else {
    // both points are on the same side
    if (side1 == GEO::POSITIVE) {
      // both points are in the voronoi cell
      q.push_back( e0 );
      q.push_back( e1 );

      q0 = e0;
      q1 = e1;
    }
    else {
      // both points are outside the voronoi cell
      q0 = -1;
      q1 = -1;
    }
  }
  return on_bisector;
}

bool
Cell::security_radius_reached( index_t bj ) const
{
  const coord_t dim = delaunay_.dim();
  real_t R = -1.;
  real_t d;
  for (index_t k = 0; k < polytope_.size(); k++) {
    d = numerics::distance2( vertex_[polytope_[k]].X() , delaunay_[site_] , dim );
    if (d > R) R = d;
  }
  // a little more than 2^2 (4.1) because we need to make sure that we clip
  // against all possible bisectors (this just means we might try to clip
  // with a non-contributing cell which will terminate anyway)
  if (numerics::distance2( delaunay_[site_] , delaunay_[bj] , dim ) > 4.1*R) {
    // the security radius has been reached, signal we need to stop clipping
    return true;
  }
  return false;
}

}

} // avro
