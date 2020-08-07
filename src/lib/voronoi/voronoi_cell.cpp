#include "common/parallel_for.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi_cell.h"


namespace avro
{

namespace delaunay
{

VoronoiCell::VoronoiCell( index_t site , const Delaunay& delaunay , const NearestNeighbours& neighbours , const TopologyBase& domain , bool exact , bool simplex ) :
  Topology<Polytope>(points_,domain.number()),
  site_(site),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  neighbours_(neighbours),
  domain_(domain),
  exact_(exact),
  facets_(nullptr),
  simplex_(simplex)
{
  avro_assert( domain.number() == points_.dim() );
  initialize();
}

void
VoronoiCell::initialize()
{
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
VoronoiCell::compute()
{
  if (simplex_) compute_simplex();
  else compute_polytope();
}

void
VoronoiCell::compute_polytope()
{
  index_t j = 1;
  while (true)
  {
    index_t bj = neighbours_(site_,j);
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
    printf("clipping with simplex %lu\n",k);

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

      vertex_[j].print("v",true);
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
      index_t bj = neighbours_(site_,j);
      clip_by_bisector( j , bj );

      if (security_radius_reached(bj)) break;
      j++;
      if (j == delaunay_.nb()) break;
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
    print_inline(polytope_);

    add( polytope_.data() , polytope_.size() );

  }

}

void
VoronoiCell::clip_by_bisector( index_t j , index_t bj )
{
  // retrieve the bisector
  //int b = delaunay_.bisector( site_ , bj );
  int b = delaunay_.bisector( 0 , j ); // local retrieval

  // initialize the clipped polytope
  std::vector<index_t> q;

  // retrieve the edges
  for (index_t ii=0;ii<polytope_.size();ii++)
  for (index_t jj=ii+1;jj<polytope_.size();jj++)
  {
    index_t e0 = polytope_[ii];
    index_t e1 = polytope_[jj];
    if (element().is_edge( vertex_[e0].bisectors() , vertex_[e1].bisectors() ) )
    {
      // clip the edge and save the result into q
      clip_edge(e0,e1,b,q);
    }
  }

  // the current polytope becomes the clipped one
  polytope_ = q;
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

  pi = neighbours_(site_,pi);
  pj = neighbours_(site_,pj);

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

VoronoiDiagram::VoronoiDiagram( Delaunay& delaunay , const TopologyBase& domain , bool simplex ) :
  Topology<Polytope>(points_,domain.number()),
  points_(delaunay.dim()),
  delaunay_(delaunay),
  domain_(domain),
  simplex_(simplex)
{}

void
VoronoiDiagram::compute( bool exact )
{
  std::shared_ptr<RVDFacets> facets;
  if (simplex_)
    facets = std::make_shared<RVDFacets>( static_cast<const Topology<Simplex>&>(domain_) );

  NearestNeighbours neighbours(delaunay_,100);
  neighbours.compute();
  printf("computed neighbours!\n");

  for (index_t k=0;k<delaunay_.nb();k++)
  {
    cells_.push_back( std::make_shared<VoronoiCell>(k,delaunay_,neighbours,domain_,exact,simplex_) );
    cells_[k]->set_facets( facets.get() );
  }

  #if 0
  ProcessCPU::parallel_for(
    parallel_for_member_callback( this , &thisclass::clip ),
    0,cells_.size()
  );
  #else
  for (index_t k=0;k<cells_.size();k++)
    cells_[k]->compute();
  #endif

  // accumulate the result
  for (index_t k=0;k<cells_.size();k++)
  {
    const VoronoiCell& cell = *cells_[k].get();

    // add the cells
    avro_assert( cell.nb() >= 1 );
    for (index_t j=0;j<cell.nb();j++)
    {
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
}

} // delaunay


} // avro
