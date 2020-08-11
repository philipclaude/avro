#include "common/parallel_for.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi_cell.h"

#define SAFE_HASH 1

namespace avro
{

namespace delaunay
{

// needed to create a set/map of elements
bool
operator==( const Bisector& bx , const Bisector& by )
{
  if (bx.p0 != by.p0) return false;
  if (bx.p1 != by.p1) return false;
  return true;
}

// needed to create a map of elements
bool
operator<( const Bisector& f , const Bisector& g )
{
  if (f.p0 < g.p0) return true;
  if (f.p0 == g.p0)
  {
    if (f.p1 < g.p1) return true;
  }
  return false;
}

VoronoiCell::VoronoiCell( index_t site , const Delaunay& delaunay , const NearestNeighbours& neighbours ,
                          const TopologyBase& domain , bool exact , bool simplex ,
                          GEO::NearestNeighborSearch* nns , index_t nb_nns ) :
  Topology<Polytope>(points_,domain.number()),
  site_(site),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  neighbours_(neighbours),
  domain_(domain),
  exact_(exact),
  facets_(nullptr),
  simplex_(simplex),
  nns_(nns),
  nn_(nb_nns,0),
  recycle_neighbours_(false),
  simplex_points_(delaunay_.dim()),
  simplices_(simplex_points_,domain.number())
{}

void
VoronoiCell::initialize()
{
  clear();
  points_.clear();
  simplices_.clear();
  simplex_points_.clear();

  if (simplex_) initialize_simplex();
  else initialize_polytope();
}

void
VoronoiCell::initialize_polytope()
{
  // create the initial points
  vertex_.resize( domain_.nv(0) );
  for (index_t j=0;j<domain_.nv(0);j++)
  {
    vertex_[j] = Vertex(delaunay_.dim(),domain_.number());
    vertex_[j].setCoordinates( domain_.points()[ domain_(0)[j] ] , domain_.points().dim() );
    vertex_[j].setNumber( domain_.number() );
    vertex_[j].addSimplexVertex( domain_.points()[ domain_(0,j) ] );
    std::vector<int> b = domain_.points().incidence().get( domain_(0,j) );
    for (index_t i=0;i<b.size();i++)
      vertex_[j].addBisector(b[i]);
  }

  avro_assert( domain_.nb() == 1 );
  polytope_ = linspace(domain_.nv(0));
}

void
VoronoiCell::initialize_simplex()
{
}

void
VoronoiCell::enlarge_neighbours()
{
  avro_assert( nns_ != nullptr );
  index_t nb_nns = nn_.size();
  nb_nns += 10;
  std::vector<real_t> dist2(nb_nns); // not used
  nn_.resize( nb_nns );
  clock_t t0 = clock();
  nns_->get_nearest_neighbors( nb_nns , delaunay_[site_] , nn_.data() , dist2.data()  );
  time_neighbours_ += real_t(clock()-t0)/real_t(CLOCKS_PER_SEC);
}

void
VoronoiCell::compute()
{
  time_neighbours_ = 0;
  time_clip_       = 0;
  time_decompose_  = 0;

  initialize();
  if (nns_ != nullptr)
  {
    std::vector<double> dist2(nn_.size(),0.0);
    clock_t t0 = clock();
    if (!recycle_neighbours_)
    {
      nns_->get_nearest_neighbors( nn_.size() , site_ , nn_.data() , dist2.data() );
      //recycle_neighbours_ = true; // seems to be a bug when keeping initial values..
    }
    else
    {
      avro_implement; // there's a bug
      nns_->get_nearest_neighbors( nn_.size() , delaunay_[site_] , nn_.data() , dist2.data() , GEO::NearestNeighborSearch::KeepInitialValues() );
    }
    avro_assert( nn_[0] == site_ );
    time_neighbours_ += real_t(clock() -t0) / real_t(CLOCKS_PER_SEC);
  }
  else
  {
    avro_assert_msg( nn_.size() == neighbours_.knear() , "|nn| = %lu, knear = %lu\n",nn_.size(),neighbours_.knear() );
    for (index_t k=0;k<nn_.size();k++)
      nn_[k] = neighbours_(site_,k);
  }

  clock_t t0 = clock();
  if (simplex_) compute_simplex();
  else compute_polytope();
  time_clip_ += real_t(clock() - t0)/real_t(CLOCKS_PER_SEC);

  t0 = clock();
  generate_simplices();
  time_decompose_ = real_t(clock() - t0)/real_t(CLOCKS_PER_SEC);
}

void
VoronoiCell::compute_polytope()
{
  index_t j = 1;
  while (true)
  {
    //index_t bj = neighbours_(site_,j);
    index_t bj = nn_[j];
    clip_by_bisector( j , bj );

    if (security_radius_reached(bj)) break;
    j++;
    if (j == delaunay_.nb()) break;
  }

  // add the points
  for (index_t j=0;j<polytope_.size();j++)
  {
    // save the geometry and bisector information
    points_.create( vertex_[polytope_[j]].X() );
    const std::vector<int>& bisectors = vertex_[ polytope_[j]].bisectors();
    points_.incidence().add( bisectors.data() , bisectors.size() );
  }
  polytope_ = linspace( points_.nb() );
  add( polytope_.data() , polytope_.size() );
}

void
VoronoiCell::compute_simplex()
{
  avro_assert( facets_ != nullptr );

  for (index_t k=0;k<domain_.nb();k++)
  {
    vertex_.clear();
    bisector_.clear();
    ids_.clear();

    // create the initial points
    polytope_ = linspace(domain_.nv(k));
    vertex_.resize( domain_.nv(k) );
    for (index_t j=0;j<domain_.nv(k);j++)
    {
      vertex_[j] = Vertex(delaunay_.dim(),number_);
      vertex_[j].setCoordinates( domain_.points()[ domain_(k)[j] ] ,
                  domain_.points().dim() );
      vertex_[j].setNumber( domain_.number() );
      vertex_[j].addSimplexVertex( domain_.points()[ domain_(k,j) ] );
    }

    // create the initial simplex
    std::vector<index_t> simplex_k = domain_.get(k);

    // add the bounding facets as bisectors
    index_t nf = domain_.number()+1;
    for (index_t j=0;j<nf;j++)
    {
      std::vector<index_t> facet = simplex_k;
      facet.erase( facet.begin() +j );
      std::sort( facet.begin() , facet.end() );

      // get the facet label
      int b = facets_->facet(facet);

      for (index_t i=0;i<simplex_k.size();i++)
      {
        if (i==j) continue; // skip the vertex opposite the facet
        vertex_[i].addBisector(b);
      }
    }

    // clip the simplex
    index_t j = 1;
    while (true)
    {
      //index_t bj = neighbours_(site_,j);
      index_t bj = nn_[j];
      clip_by_bisector( j , bj );

      if (security_radius_reached(bj)) break;
      j++;
      if (j == delaunay_.nb()) break;
      if (j == nn_.size())
      {
        // max nearest neighbours
        enlarge_neighbours();
      }
    }
    if (polytope_.size() == 0 ) continue;

    // add the points
    index_t nb_points = points_.nb();
    for (index_t j=0;j<polytope_.size();j++)
    {
      // save the geometry and bisector information
      points_.create( vertex_[polytope_[j]].X() );
      const std::vector<int>& bisectors = vertex_[ polytope_[j]].bisectors();
      points_.incidence().add( bisectors.data() , bisectors.size() );
    }
    polytope_ = linspace( polytope_.size() );
    for (index_t j=0;j<polytope_.size();j++)
      polytope_[j] += nb_points;
    add( polytope_.data() , polytope_.size() );
  }
}

int
VoronoiCell::add_bisector( index_t p0 , index_t p1 )
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

void
VoronoiCell::get_bisector( int b , index_t& p0 , index_t& p1 ) const
{
  std::map<int,Bisector>::const_iterator it;
  it = ids_.find(b);
  avro_assert( it != ids_.end() );
  p0 = it->second.p0;
  p1 = it->second.p1;
}

void
VoronoiCell::clip_by_bisector( index_t j , index_t bj )
{
  // retrieve the bisector
  #if SAFE_HASH
  //int b = delaunay_.bisector( 0 , j ); // local retrieval
  int b = add_bisector( site_ , bj );
  #else
  int b = delaunay_.bisector( site_ , bj );
  #endif

  // initialize the clipped polytope
  qpolytope_.clear();

  // retrieve the edges
  for (index_t ii=0;ii<polytope_.size();ii++)
  for (index_t jj=ii+1;jj<polytope_.size();jj++)
  {
    index_t e0 = polytope_[ii];
    index_t e1 = polytope_[jj];
    if (element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) )
    {
      // clip the edge and save the result into q
      clip_edge(e0,e1,b,qpolytope_);
    }
  }

  // the current polytope becomes the clipped one
  polytope_.assign( qpolytope_.begin() , qpolytope_.end() );
  //polytope_ = q;
  uniquify(polytope_);
  for (index_t ii=0;ii<polytope_.size();ii++)
    vertex_[ polytope_[ii] ].setBaseSite( delaunay_[site_] );
}

void
VoronoiCell::clip_edge( const index_t e0 , const index_t e1 , const int b , std::vector<index_t>& q )
{
  Vertex& v0 = vertex_[e0];
  Vertex& v1 = vertex_[e1];

  index_t pi,pj;
  delaunay_.seeds(b,pi,pj);

  #if SAFE_HASH
  //pi = neighbours_(site_,pi);
  //pj = neighbours_(site_,pj);
  get_bisector( b , pi , pj );
  #endif

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

  avro_assert( side1!=GEO::ZERO && side2!=GEO::ZERO );

  if (side1!=side2)
  {
    Vertex v2(delaunay_.dim(),domain_.number());

    // perform the symbolic intersection of the bisectors, simplices, meshes
    v2.setBaseSite(delaunay_[site_]);
    v2.intersectSymbolic( &v0 , &v1 , delaunay_ );
    v2.setSites( delaunay_ , ids_ );

    const std::vector<const real_t*>& sites = v2.sites();
    for (index_t s=0;s<sites.size();s++)
    {
      //avro_assert( sites[s] != zj );
      if (sites[s] == zj) printf("conflict %p\n",(void*)zj);
    }

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
      q.push_back( e0 );
    else
      q.push_back( e1 );

    q.push_back(vs);

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
    }
    else
    {
      // both points are outside the voronoi cell
    }
  }
}

bool
VoronoiCell::security_radius_reached( index_t bj ) const
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
  if (numerics::distance2( delaunay_[site_] , delaunay_[bj] , dim )>4.1*R) return true;
  return false;
}

void
VoronoiCell::generate_simplices()
{
  simplices_.clear();
  simplex_points_.clear();

  coord_t dim = delaunay_.dim();

  points_.copy( simplex_points_ );

  std::vector<real_t> xc(dim,0.);
  std::vector<int> facets;
  std::vector<index_t> simplex(number_+1);

  if (number_ == 2)
  {
    // add a point for the centroid of each polygon and triangulate
    for (index_t k=0;k<nb();k++)
    {
      index_t idc = simplex_points_.nb();
      numerics::centroid( (*this)(k) , nv(k) , points_ , xc );
      simplex_points_.create( xc.data() );

      // get the hrep of this polygon
      element().hrep( (*this)(k) ,  nv(k) , facets );

      avro_assert( facets.size() == nv(k) ); // for polygons
      for (index_t j=0;j<facets.size();j++)
      {
        // get the points with this bisector
        std::vector<index_t> vf;
        element().vrep( (*this)(k) , nv(k) , facets[j] , vf );
        avro_assert( vf.size() == 2 ); // for polygons

        simplex[0] = vf[0];
        simplex[1] = vf[1];
        simplex[2] = idc;

        simplices_.add( simplex.data() , simplex.size() );
      }
    }
  }
  else if (number_ == 3)
  {
    // add a point for the centroid of each polygon and triangulate
    std::vector<real_t> xf(dim,0.);
    std::vector<index_t> edges;
    for (index_t k=0;k<nb();k++)
    {
      index_t idc = simplex_points_.nb();
      numerics::centroid( (*this)(k) , nv(k) , points_ , xc );
      simplex_points_.create( xc.data() );

      // get the hrep of this polygon
      element().hrep( (*this)(k) ,  nv(k) , facets );

      for (index_t j=0;j<facets.size();j++)
      {
        // get the points with this bisector
        std::vector<index_t> vf;
        element().vrep( (*this)(k) , nv(k) , facets[j] , vf );

        // compute the centroid of the facet
        index_t idf = simplex_points_.nb();
        numerics::centroid( vf.data() , vf.size() , points_ , xf );
        simplex_points_.create( xf.data() );

        // retrieve the edges
        std::vector<int> edges;
        element().hrep( vf.data() , vf.size() , edges );
        for (index_t e=0;e<edges.size();e++)
        {
          if (edges[e] == facets[j]) continue;

          std::vector<index_t> ve;
          element().vrep( vf.data() , vf.size() , edges[e] , ve );
          avro_assert( ve.size() == 2 );

          simplex[0] = ve[0];
          simplex[1] = ve[1];
          simplex[2] = idf;
          simplex[3] = idc;
          simplices_.add( simplex.data() , simplex.size() );
        }
      }
    }
  }
  else
    avro_implement;
}

VoronoiDiagram::VoronoiDiagram( Delaunay& delaunay , const TopologyBase& domain , bool simplex ) :
  Topology<Polytope>(points_,domain.number()),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  domain_(domain),
  cells_(delaunay.nb(),nullptr),
  simplex_(simplex),
  nns_(nullptr)
{}

void
VoronoiDiagram::compute( bool exact )
{
  std::shared_ptr<RVDFacets> facets;
  if (simplex_)
    facets = std::make_shared<RVDFacets>( static_cast<const Topology<Simplex>&>(domain_) );

  real_t t0 = clock();
  index_t nb_nns = 50;
  if (delaunay_.nb() < nb_nns) nb_nns = delaunay_.nb();
  NearestNeighbours neighbours(delaunay_,nb_nns);

  #if 0
  neighbours.compute();
  #else
  const coord_t dim = delaunay_.dim();
  if (nns_ == nullptr)
    nns_ = GEO::NearestNeighborSearch::create(dim,"ANN");
  std::vector<real_t> x(delaunay_.nb()*dim);
  for (index_t k=0;k<delaunay_.nb();k++)
  for (index_t d=0;d<dim;d++)
    x[k*dim+d] = delaunay_[k][d];
  t0 = clock();
  nns_->set_points( delaunay_.nb() , x.data() );
  #endif
  time_neighbours_ = real_t( clock() - t0 )/real_t(CLOCKS_PER_SEC);

  for (index_t k=0;k<delaunay_.nb();k++)
  {
    if (cells_[k] == nullptr || nns_ == nullptr)
      cells_[k] = std::make_shared<VoronoiCell>(k,delaunay_,neighbours,domain_,exact,simplex_,nns_,nb_nns );
    cells_[k]->set_facets( facets.get() );
  }

  //t0 = clock();
  #if 1
  ProcessCPU::parallel_for(
    parallel_for_member_callback( this , &thisclass::clip ),
    0,cells_.size()
  );
  #else
  for (index_t k=0;k<cells_.size();k++)
    cells_[k]->compute();
  #endif
  //time_voronoi_ = real_t(clock() - t0 )/real_t(CLOCKS_PER_SEC);

  time_voronoi_ = 0;
  real_t time_decompose = 0;

  // accumulate the result
  sites_.clear();
  for (index_t k=0;k<cells_.size();k++)
  {
    const VoronoiCell& cell = *cells_[k].get();

    time_neighbours_ += cell.time_neighbours() / ProcessCPU::maximum_concurrent_threads();
    time_voronoi_    += cell.time_clip() / ProcessCPU::maximum_concurrent_threads();
    time_decompose   += cell.time_decompose() / ProcessCPU::maximum_concurrent_threads();

    // cells could be submerged
    if (cell.nb() == 0) continue;

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
      points_.create( cell.points()[j] );

    // add the vertex-to-facet information
    const Table<int>& vf = cell.element().incidence();
    for (index_t j=0;j<vf.nb();j++)
    {
      std::vector<int> f = vf.get(j);
      points_.incidence().add( f.data() , f.size() );
    }

  }

  printf("--> timing: neighbours (%3.4g sec.), clipping (%3.4g sec.), decompose (%3.4g sec.)\n",time_neighbours_,time_voronoi_,time_decompose);
}

} // delaunay


} // avro
