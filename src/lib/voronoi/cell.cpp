#include "mesh/points.h"
#include "numerics/geometry.h"
#include "numerics/mat.hpp"
#include "element/quadrature.h"
#include "numerics/vec.hpp"
#include "voronoi/cell.h"

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
  exact_(false),
  ambient_dim_(domain.number())
{}

void
Cell::compute( const std::vector<index_t>& candidates ) {

  clear();

  // compute nearest neighbours
  index_t nb_nns = 50;
  if (delaunay_.nb() < nb_nns)
    neighbours_.resize(delaunay_.nb());
  else
    neighbours_.resize(nb_nns);
  std::vector<double> dist2(neighbours_.size(),0.0);
  search_.get_nearest_neighbors( neighbours_.size() , site_ , neighbours_.data() , dist2.data() );

  // do the clipping
  for (index_t k = 0; k < candidates.size(); k++) {
    // break as soon as we found an element (from the candidates) to clip with
    // it's possible that there is no intersection when weights are used
    if (clip_simplex(candidates[k])) {
      break;
    }
  }

  // finish off the facet calculation
  finish_facets();
}

void
Cell::clear() {
  Topology<Polytope>::clear();
  points_.clear();
  visited_.clear();
  triangles_.clear();
  edges_.clear();
  bisector_.clear();
  ids_.clear();
  moment_.resize( ambient_dim_ );
  for (coord_t d = 0; d < moment_.size(); d++)
    moment_[d] = 0.0;
  volume_ = 0.0;
  boundary_area_ = 0.0;
  facets_.clear();
  energy_ = 0.0;
}

bool
Cell::clip_simplex( index_t elem ) {

  // check if this element has already been visited
  if (visited_.find(elem) != visited_.end()) return false;

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

  if (polytope_.size() == 0) return false;

  // decompose the simplices (integration and visualization)
  generate_simplices();

  // add the points
  index_t nb_points = points_.nb();
  for (index_t j = 0; j < polytope_.size(); j++) {
    // save the coordinates and bisector information
    points_.create( vertex_[polytope_[j]].X() );
    polytope_[j] = nb_points + j;
  }

  // add the polytope
  add( polytope_.data() , polytope_.size() );

  // move to the neighbours of the current element
  for (int j = 0; j < domain_.number()+1; j++) {
    int n = domain_.neighbours()(elem,j);
    if (n < 0) continue; // do not step into boundaries
    clip_simplex( index_t(n) );
  }

  return true;
}

void
Cell::initialize( index_t elem ) {

  const coord_t dim = delaunay_.dim();
  const coord_t n = domain_.number();
  elem_ = elem;

  // get the list of facets of this element
  std::vector<int> facets(n+1);
  for (int j = 0; j < n+1; j++) {
    int id = encode_mesh_facet(j);
    facets[j] = id;
  }

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

  // initialize the polytope to the simplex
  polytope_.resize(n+1);
  for (int j = 0; j < n+1; j++)
    polytope_[j] = j;

  // initialize the polytope edges
  pedges_.clear();
  for (int ii = 0; ii < n+1; ii++)
  for (int jj = ii+1; jj < n+1; jj++) {
    pedges_.push_back(ii);
    pedges_.push_back(jj);
  }
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
fan_triangulation( const std::vector<index_t>& edges , std::vector<index_t>& triangles ) {

  const index_t nb_edges = edges.size()/2;

  // make enough room for the new triangles
  index_t n = triangles.size();
  triangles.resize( n + (nb_edges-2)*3 );

  // connect everything to the first point
  index_t p = edges[0];
  index_t t = 0;
  for (index_t k = 1; k < nb_edges; k++) {

    index_t e0 = edges[2*k  ];
    index_t e1 = edges[2*k+1];

    if (e0 == p || e1 == p) continue;

    triangles[n + 3*t    ] = p;
    triangles[n + 3*t + 1] = e0;
    triangles[n + 3*t + 2] = e1;
    t++;
  }
}

real_t
Cell::compute_energy( const std::vector<const real_t*>& X ) const {

  coord_t n = domain_.number();
  avro_assert( X.size() == index_t(n+1) );

  int quad_order = 2;
  const matd<real_t>& B = __store_simplex_lagrange__.get_basis(n,1,quad_order);
  const Quadrature& Q = __store_simplex_lagrange__.quadrature(n,quad_order);
  avro_assert( B.m() == X.size() );

  real_t e = 0.0;
  for (index_t k = 0; k < B.n(); k++) {

    // evaluate the coordinate at this quadrature point
    std::vector<real_t> p(ambient_dim_,0.0);
    for (index_t i = 0; i < B.m(); i++)
    for (coord_t d = 0; d < ambient_dim_; d++)
      p[d] += B(i,k)*X[i][d];

    // evaluate the energy integrand
    real_t w = Q.w(k); // quadrature weight
    e += w*numerics::distance2( p.data() , delaunay_[site_] , ambient_dim_ );
  }

  return e;
}

void
Cell::generate_simplices() {

  // compute the map from the points to what they will be in the final saved polytope
  index_t nb_points = points_.nb();
  std::map<index_t,index_t> point_map;
  for (index_t k = 0; k < polytope_.size(); k++)
    point_map.insert( {polytope_[k] , nb_points+k} );

  if (number_ == 2) {

    // produce a fan triangulation
    std::vector<index_t> ptriangles;
    fan_triangulation( pedges_ , ptriangles );

    // add the triangles
    for (index_t j = 0; j < ptriangles.size(); j++)
      triangles_.push_back( point_map.at(ptriangles[j]) );

    // add the edges
    for (index_t j = 0; j < pedges_.size(); j++)
      edges_.push_back( point_map.at(pedges_[j]) );

    // add the volume terms
    // if we just want the volume and first moment, then calculate these on the fly and save the result
    std::vector<const real_t*> X(3);
    for (index_t k = 0; k < ptriangles.size()/3; k++) {

      // compute the volume of the triangle
      for (coord_t i = 0; i < 3; i++)
        X[i] = vertex_[ptriangles[3*k+i]].X();

      // add the volume term
      real_t vk = fabs(numerics::simplex_volume( X , ambient_dim_ , false ));
      volume_ += vk;

      // add the energy
      energy_ += compute_energy(X)*vk;

      // compute the centroid of this triangle
      std::vector<real_t> ck( ambient_dim_ , 0.0 );
      for (coord_t i = 0; i < X.size(); i++)
      for (coord_t d = 0; d < ambient_dim_; d++)
        ck[d] += X[i][d]/X.size();

      // compute the first moment
      for (index_t d = 0; d < ambient_dim_; d++)
        moment_[d] += ck[d]*vk;
    }

    // add the face terms
    // i.e. compute dij: the distance between this site and the bisector
    //              lij: the distance between this site and the opposite site
    //              Aij: the area of the polygon on this bisector
    //              bij: the centroid of the polygon on the bisector
    for (index_t k = 0; k < pedges_.size()/2; k++) {

      index_t e0 = pedges_[2*k];
      index_t e1 = pedges_[2*k+1];

      // bisectors of each vertex
      const std::vector<int>& b0 = vertex_[e0].bisectors();
      const std::vector<int>& b1 = vertex_[e1].bisectors();
      avro_assert( b0.size() == 2 && b1.size() == 2 );

      // determine which bisector this edge is on
      int h;
      if      (b0[0] == b1[0]) h = b0[0];
      else if (b0[0] == b1[1]) h = b0[0];
      else if (b0[1] == b1[0]) h = b0[1];
      else if (b0[1] == b1[1]) h = b0[1];
      else avro_assert_not_reached;

      // we only save power facets if they are voronoi bisectors, or if they are boundary mesh facets
      bool interior = false;
      if (h < 0) {
        index_t f = decode_mesh_facet(elem_,h);
        if (domain_.neighbours()(elem_,f) >= 0) interior = true;
      }

      // this should be computed in the ambient dimension (either 2d or 3d)
      avro_assert( ambient_dim_ == 2 || ambient_dim_ == 3);
      if (!interior) {

        // compute the length of this edge
        real_t Aij = numerics::distance( vertex_[e0].X() , vertex_[e1].X() , ambient_dim_ );

        // the last coordinate will be 0 if the ambient dimension is 2d
        vecs<3,real_t> x0,x1,x2;
        for (coord_t d = 0; d < ambient_dim_; d++) {
          x0(d) = delaunay_[site_][d];
          x1(d) = vertex_[e0][d];
          x2(d) = vertex_[e1][d];
        }
        real_t dij = norm( numerics::cross(x0-x1,x0-x2) ) / norm(x2-x1);
        vecs<3,real_t> abij = Aij * 0.5 * ( x1 + x2 );

        if (facets_.find(h) == facets_.end()) {

          PowerFacet facet;

          facet.Aij = Aij;
          facet.dij = dij;
          facet.Abij.resize(ambient_dim_);
          facet.bij.resize(ambient_dim_);
          facet.normal.resize(ambient_dim_);
          facet.elem = elem_;

          for (coord_t d = 0; d < ambient_dim_; d++) {
            facet.Abij[d]   = abij(d);
            facet.normal[d] = (x2(d) - x1(d))/Aij; // i'm not really sure about this one
          }

          facets_.insert( {h,facet} );
        }
        else {

          PowerFacet& facet = facets_.at(h);

          facet.Aij += Aij;
          for (coord_t d = 0; d < ambient_dim_; d++)
            facet.Abij[d] += abij(d);

          // should we check that dij is close to the previous value?
          // this may not be accurate for general 3d meshes
          //avro_assert_msg( fabs(facet.dij - dij) < 1e-8 , "|dij - f.dij| = %g" , fabs(facet.dij-dij) );
        }
      }
    }
  }
  else if (number_ == 3) {

    // compute the centroid of the polyhedrons
    std::vector<real_t> centroid(delaunay_.dim());
    for (index_t k = 0; k < polytope_.size(); k++)
    for (coord_t d = 0; d < delaunay_.dim(); d++)
      centroid[d] += vertex_[polytope_[k]][d] / polytope_.size();

    // produce a fan triangulation of the polygon on each bisector
    std::set<int> bisectors;
    for (index_t k = 0; k < polytope_.size(); k++) {
      const std::vector<int>& b = vertex_[polytope_[k]].bisectors();
      for (index_t j = 0; j < b.size(); j++) {
        bisectors.insert( b[j] );
      }
    }

    for (std::set<int>::iterator it = bisectors.begin(); it != bisectors.end(); ++it) {

      int h = *it;

      // retrieve the edges along this bisector
      std::vector<index_t> bedges;
      for (index_t ii = 0; ii < pedges_.size()/2; ii++) {

        // both vertices must have this bisector for this to be an edge on the bisector
        index_t e0 = pedges_[2*ii];
        index_t e1 = pedges_[2*ii+1];

        bool on_bisector = false;
        const std::vector<int>& b0 = vertex_[e0].bisectors();
        for (index_t jj = 0; jj < b0.size(); jj++) {
          if (b0[jj] == h) {
            on_bisector = true;
            break;
          }
        }
        if (!on_bisector) continue;

        on_bisector = false;
        const std::vector<int>& b1 = vertex_[e1].bisectors();
        for (index_t jj = 0; jj < b1.size(); jj++) {
          if (b1[jj] == h) {
            on_bisector = true;
            break;
          }
        }
        if (!on_bisector) continue;

        bedges.push_back(e0);
        bedges.push_back(e1);
      }

      // we only save visualization triangles/edges if they are voronoi bisectors, or if they are boundary mesh facets
      bool interior = false;
      if (h < 0) {
        index_t f = decode_mesh_facet(elem_,h);
        if (domain_.neighbours()(elem_,f) >= 0) interior = true;
      }

      // compute the fan triangulation of the polygon on this bisector
      std::vector<index_t> ptriangles;
      fan_triangulation( bedges , ptriangles );

      if (!interior) {
        // save the visualization triangles
        for (index_t j = 0; j < ptriangles.size(); j++)
          triangles_.push_back( point_map.at(ptriangles[j]) );

        // save the visualization edges
        for (index_t j = 0; j < pedges_.size(); j++)
          edges_.push_back( point_map.at(pedges_[j]) );
      }

      // if we just want the volume and first moment, then calculate these on the fly and save the result
      std::vector<const real_t*> X(4);
      std::vector<const real_t*> Y(3);
      X[0] = centroid.data();
      real_t Aij = 0.0;
      int k_max = -1;
      real_t aij_max = -1;
      vecs<3,real_t> abij;
      for (index_t k = 0; k < ptriangles.size()/3; k++) {

        // compute the volume of the tetrahedron
        for (coord_t i = 0; i < 3; i++) {
          X[i+1] = vertex_[ptriangles[3*k+i]].X();
          Y[i]   = vertex_[ptriangles[3*k+i]].X();
        }

        // add the volume term
        real_t vk = fabs(numerics::simplex_volume( X , 3 , false ));
        volume_ += vk;
        energy_ += compute_energy(X)*vk;

        // compute the area in 3d
        real_t aij = fabs(numerics::simplex_volume( Y , 3 , false ) );
        Aij += aij;

        if (aij > aij_max) {
          // compute a basis using this triangle
          k_max   = k;
          aij_max = aij;
        }

        // compute the centroid of this tetrahedron
        std::vector<real_t> ck( ambient_dim_ , 0.0 );
        for (coord_t i = 0; i < 4; i++)
        for (coord_t d = 0; d < ambient_dim_; d++)
          ck[d] += 0.25*X[i][d];

        // compute the first moment
        for (index_t d = 0; d < ambient_dim_; d++)
          moment_[d] += ck[d]*vk;

        // compute the centroid of this triangle
        for (index_t d = 0; d < 3; d++)
          abij(d) += aij*(Y[0][d] + Y[1][d] + Y[2][d])/3.0;
      }
      if (aij_max == 0.0) {
        for (index_t j = 0; j < polytope_.size(); j++)
          vertex_[polytope_[j]].print("v");
      }
      // it may still be possible to have a zero area due to precision issues
      avro_assert_msg( aij_max > 0.0 , "aij_max = %g, nb_triangles = %lu" , aij_max , ptriangles.size()/3 );

      avro_assert( ambient_dim_ == 3 );
      vecs<3,real_t> u,v,n,c,qi;
      index_t i0 = ptriangles[3*k_max];
      index_t i1 = ptriangles[3*k_max+1];
      index_t i2 = ptriangles[3*k_max+2];

      for (coord_t d = 0; d < ambient_dim_; d++) {
        u(d)  = vertex_[i1][d] - vertex_[i0][d];
        v(d)  = vertex_[i2][d] - vertex_[i0][d];
        c(d)  = vertex_[i0][d];
        qi(d) = delaunay_[site_][d];
      }
      n = numerics::cross(u,v);
      numerics::normalize(n);

      // distance from the site to the plane
      real_t dij = fabs(dot( n , qi - c ));

      // compute the face terms in the actual ambient dimension that we will need to compute the discrete operators
      // i.e. compute dij: the distance between this site and the bisector
      //              lij: the distance between this site and the opposite site
      //              Aij: the area of the polygon on this bisector
      //              bij: the centroid of the polygon on the bisector
      //              normal: normal to the facet

      // this should be computed in 3d
      if (!interior) {
        if (facets_.find(h) == facets_.end()) {
          PowerFacet facet;
          facet.elem = elem_;

          facet.normal.resize(ambient_dim_);
          for (coord_t d = 0; d < ambient_dim_; d++)
            facet.normal[d] = n(d);

          facet.Aij = Aij;
          facet.dij = dij;
          facet.Abij.resize(ambient_dim_);
          for (coord_t d = 0; d < ambient_dim_; d++)
            facet.Abij[d] = abij(d);
          facets_.insert( {h,facet} );
        }
        else {
          PowerFacet& facet = facets_.at(h);

          // add to the Aij term
          facet.Aij += Aij;

          // and to the moment term
          for (coord_t d = 0; d < ambient_dim_; d++)
            facet.Abij[d] += abij(d);

          // should we check that dij is close to the previous value?
          avro_assert_msg( fabs(facet.dij - dij) < 1e-8 , "|dij - f.dij| = %g" , fabs(facet.dij-dij) );
        }
      }
    }
  }
  else if (number_ == 4) {
    // we can decompose the polytope into integration pentatopes
    // but retrieving the visualization triangles will be a little more involved
    // i don't plan on needing this in the near future - for higher-dimensional domains, use the implementation in optimal_transport.cpp
    avro_implement;
  }
  else
    avro_assert_not_reached;
}

void
Cell::finish_facets() {

  boundary_area_ = 0.0;
  for (std::map<int,PowerFacet>::iterator it = facets_.begin(); it != facets_.end(); ++it) {

    // retrieve the facet
    int h = it->first;
    PowerFacet& facet = it->second;

    // save the site information
    facet.pi = site_;

    // save the opposite site information, as well as the distance lij
    if (h >= 0) {
      index_t pi, pj;
      get_bisector(h,pi,pj);
      facet.lij = numerics::distance( delaunay_[pi] , delaunay_[pj] , ambient_dim_ );

      if (pi == site_) facet.pj = pj;
      else {
        avro_assert( pj == site_ );
        facet.pj = pi;
      }
    }
    else {
      facet.lij = facet.dij;
      facet.pj  = h;
    }

    // finish off the centroid calculation by dividing by the area
    facet.bij.resize(ambient_dim_);
    for (coord_t d = 0; d < ambient_dim_; d++)
      facet.bij[d] = facet.Abij[d] / facet.Aij;

    if (h < 0) {
      // this is a mesh element facet
      index_t f = decode_mesh_facet(facet.elem,h);
      int n = domain_.neighbours()(facet.elem,f);

      // only add the contribution to the boundary area if the mesh element at this facet has no neighbour (n < 0)
      if (n < 0) boundary_area_ += facet.Aij;
    }
  }
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

  // initialize the clipped polytope
  qpolytope_.clear();
  qedges_.clear();
  qplane_.clear();

  // retrieve the edges
  for (index_t i = 0; i < pedges_.size()/2; i++){

    // retrieve edge indices
    index_t e0 = pedges_[2*i];
    index_t e1 = pedges_[2*i+1];

    // can probably turn this off since it adds computation
    avro_assert( element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) );

    // clip the edge by the bisector
    int q0 = -1, q1 = -1;
    int qs = clip_edge(e0,e1,b,qpolytope_,q0,q1);
    if (q0 < 0 || q1 < 0) continue;

    // add the vertices to the new polytope
    qedges_.push_back( q0 );
    qedges_.push_back( q1 );

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

    // perform the symbolic intersection of the bisectors and simplex vertices
    v2.set_base_site(delaunay_[site_]);
    v2.intersect_symbolic( &v0 , &v1 );
    v2.set_sites( delaunay_ , ids_ );

    // this vertex lies exactly on the slicing bisector
    v2.add_bisector( b );
    v2.add_site( zj );

    // compute the geometric intersection (i.e. the coordinates of the intersection point)
    v2.intersect_geometric( v0.X() , v1.X() , zi , zj );

    avro_assert_msg( v0.nb_bisectors() == domain_.number() ,
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

} // voronoi

} // avro
