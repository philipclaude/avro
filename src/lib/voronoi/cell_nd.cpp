#include "numerics/geometry.h"

#include "voronoi/cell_nd.h"

namespace avro
{

namespace voronoi
{

template<typename type>
LaguerreCellBase<type>::LaguerreCellBase( const Points& delaunay, GEO::NearestNeighborSearch& nns,
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

LaguerreCell<Polytope>::LaguerreCell(index_t site , const Points& delaunay ,
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
    vertex_[j].set_coordinates( domain_.points()[ domain_(elem_)[j] ] , domain_.points().dim() );
    vertex_[j].set_number( domain_.number() );
    vertex_[j].add_simplex_vertex( domain_.points()[ domain_(elem_,j) ] );
    std::vector<int> b = domain_.points().incidence().get( domain_(elem_,j) );
    for (index_t i=0;i<b.size();i++)
      vertex_[j].add_bisector(b[i]);
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
LaguerreCell<Simplex>::LaguerreCell(index_t elem , const Points& delaunay ,
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
    vertex_[j].set_coordinates( domain_.points()[ domain_(elem_)[j] ] ,
                domain_.points().dim() );
    vertex_[j].set_number( domain_.number() );
    vertex_[j].add_simplex_vertex( domain_.points()[ domain_(elem_,j) ] );
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
      vertex_[i].add_bisector(b);
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

static index_t
closest( const Points& points , const real_t* x )
{
	avro_assert( points.nb() > 0 );
	int k0 = 0;
	real_t d,dmin = numerics::distance2( x , points(0) , points.dim() );
	for (index_t k = 1; k < points.nb(); k++)
	{
		d = numerics::distance2( x , points(k) , points.dim() );
		if (d<dmin)
		{
			dmin = d;
			k0   = k;
		}
	}
	return index_t(k0);
}


void
LaguerreCell<Simplex>::clip()
{
  clock_t t0;

  // reset everything
  reset();
  clipped_.resize(delaunay_.nb(),false);

  // find the nearest delaunay site
  index_t s = closest( delaunay_ , domain_.points()[domain_(elem_,0)] );
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
    vertex_[ polytope_[ii] ].set_base_site( delaunay_[site_] );
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

template class LaguerreCellBase<Simplex>;
template class LaguerreCellBase<Polytope>;

} // voronoi

} // avro
