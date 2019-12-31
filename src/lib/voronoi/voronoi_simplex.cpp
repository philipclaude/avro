// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/array.h"

#include "mesh/geometrics.h"

#include "mesh/delaunay/voronoi.h"
#include "mesh/delaunay/delaunay.h"

#include "graph/neighbours.h"

namespace avro
{

namespace delaunay
{

RestrictedVoronoiSimplex::RestrictedVoronoiSimplex(
   const index_t k , const Topology<Simplex>& topology , RVDFacets& facets ,
   Delaunay& _delaunay , graph::Neighbours& _neighbours , bool _exact ) :
   Mesh<ConvexPolytope>(topology.vertices().dim(),topology.number()),
   facets_(facets),
   delaunay_(_delaunay),
   neighbours_(_neighbours),
   exact_(_exact),
   topology_(vertices_,topology.number())
{
  // mark all the delaunay seeds as unclipped
  clipped_.resize( delaunay_.nb() );
  std::fill( clipped_.begin() , clipped_.end() , false );

  // create the initial vertices
  vertex_.resize( topology.nv(k) );
  for (index_t j=0;j<topology.nv(k);j++)
  {
    vertex_[j] = Vertex(delaunay_.dim(),number_);
    vertex_[j].setCoordinates( topology.vertices()[ topology(k)[j] ] ,
                topology.vertices().dim() );
    vertex_[j].setNumber( topology.number() );
    vertex_[j].addSimplexVertex( topology.vertices()[ topology(k,j) ] );
  }

  // create the initial simplex
  simplex_ = topology.get(k);

  // add the bounding facets as bisectors
  index_t nf = topology.number()+1;
  for (index_t j=0;j<nf;j++)
  {
    std::vector<index_t> facet = simplex_;
    facet.erase( facet.begin() +j );
    std::sort( facet.begin() , facet.end() );

    // get the facet label
    int b = facets_.facet(facet);

    for (index_t i=0;i<simplex_.size();i++)
    {
      if (i==j) continue; // skip the vertex opposite the facet
      vertex_[i].addBisector(b);
    }
  }
}

index_t
RestrictedVoronoiSimplex::nb_sites() const
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
RestrictedVoronoiSimplex::nextSite()
{
  index_t k;
  for (k=0;k<sites_.size();k++)
  {
    if (clipped_[sites_[k]]) continue;
      break;
  }
  clipped_[sites_[k]] = true;
  site_ = sites_[k];
}

void
RestrictedVoronoiSimplex::addSite( const index_t zj )
{
  // we may have already clipped by this cell
  if (clipped_[zj]) return;
  sites_.push_back(zj);
}

void
RestrictedVoronoiSimplex::clip()
{
  // get the delaunay seed closest to some arbitrary vertex of the simplex
  site_ = delaunay_.closest( vertex_[0].X() );
  sites_.push_back( site_ );

  while ( nb_sites()>0 )
  {
    // reset back to the original simplex
    reset();

    // move to the next site
    nextSite();
    clip( site_ );
  }
}

void
RestrictedVoronoiSimplex::reset()
{
  // the original simplex vertices are stored at the beginning
  polytope_ = linspace( number_+1 );
}

void
RestrictedVoronoiSimplex::clip( const index_t i )
{
  clipPolytope( 0 );

  // only polytopes of an admissible size are created
  if (polytope_.size()>number_)
  {
    // we create the vertices in order so the created cell is also in order
    std::vector<index_t> cell = linspace( polytope_.size() );
    for (index_t k=0;k<cell.size();k++)
      cell[k] += vertices_.nb();

    // add the vertices
    for (index_t j=0;j<polytope_.size();j++)
    {
      // save the geometry and bisector information
      vertices_.create( vertex_[polytope_[j]].X() );
      topology_.master().vertexFacetMatrix().add( vertex_[ polytope_[j]].bisectors() );
    }
    // save the cell and which delaunay seed it belongs to
    topology_.add( cell );
    seed_.push_back(site_);
  }
}

bool
RestrictedVoronoiSimplex::securityRadiusReached( const real* zi ,
                                                 const real* zj ) const
{
  const coord_t dim = delaunay_.dim();
  real R = -1.;
  real d;
  for (index_t k=0;k<polytope_.size();k++)
  {
    d = geometrics::distance2( vertex_[polytope_[k]].X() , zi , dim );
    if (d>R) R = d;
  }
  // a little more than 2^2 (4.1) because we need to make sure that we clip
  // against all possible bisectors (this just means we might try to clip
  // with a non-contributing cell which will terminate anyway)
  if (geometrics::distance2(zi,zj,dim)>4.1*R) return true;
  return false;
}

void
RestrictedVoronoiSimplex::clipPolytope( index_t j )
{

  if (polytope_.size()==0) return;

  // get the next neighbour
  index_t zj;
  if (neighbours_(site_,j)==site_) j++;
  zj = neighbours_(site_,j);

  // TODO save the region zi-zj

  // add the cell zj as site that we will clip
  addSite( zj );

  // retrieve the bisector
  int b = delaunay_.bisector( site_ , zj );

  // initialize the clipped polytope
  std::vector<index_t> q;

  // retrieve the edges using the master
  for (index_t ii=0;ii<polytope_.size();ii++)
  for (index_t jj=ii+1;jj<polytope_.size();jj++)
  {
    index_t e0 = polytope_[ii];
    index_t e1 = polytope_[jj];
    if (topology_.master().isEdge( vertex_[e0].bisectors() ,
                        vertex_[e1].bisectors() ) )
    {
      // clip the edge and save the result into q
      clipEdge(e0,e1,b,q);
    }
  }

  // the current polytope becomes the clipped one
  polytope_ = q;
  uniquify(polytope_);

  for (index_t ii=0;ii<polytope_.size();ii++)
    vertex_[ polytope_[ii] ].setBaseSite( delaunay_[site_] );

  // check if this is the last neighbour
  if (j==neighbours_.nb(site_)-1) return;

  // check if the security readius has been reached
  if (securityRadiusReached(delaunay_[site_],delaunay_[zj]))
    return;

  // keep clipping by the next bisector
  clipPolytope( j+1 );

}

void
RestrictedVoronoiSimplex::clipEdge( const index_t e0 , const index_t e1 ,
  const int b , std::vector<index_t>& q )
{

  Vertex& v0 = vertex_[e0];
  Vertex& v1 = vertex_[e1];

  index_t pi,pj;
  delaunay_.seeds(b,pi,pj);

  const real* zi = delaunay_[site_];
  const real* zj;
  if (pi!=site_)
  {
    avro_assert(pj==site_);
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
    Vertex v2(delaunay_.dim(),number_);

    // perform the symbolic intersection of the bisectors, simplices, meshes
    v2.setBaseSite(delaunay_[site_]);
    v2.intersectSymbolic( &v0 , &v1 , delaunay_ );

    // this vertex lies exactly on the slicing bisector
    v2.addBisector( b );
    v2.addSite( zj );

    for (index_t j=0;j<v2.nb_sites();j++)
    {
      if (delaunay_[site_]==v2.site(j))
      {
        v2.print("vs",true);
        printf("zi = %p\n",(void*)delaunay_[site_]);
        avro_implement;
      }
    }

    // compute the geometric intersection
    v2.intersectGeometric( v0.X() , v1.X() , zi , zj );

    avro_assert_msg( v0.nb_bisectors()==number_ ,
                "nb_bisectors = %lu, number = %u",v0.nb_bisectors(),number_ );

    // add the vertex after saving its location
    index_t vs = vertex_.size();
    vertex_.emplace_back(v2);

    if (side1==GEO::POSITIVE)
      q.push_back( e0 );
    else
      q.push_back( e1 );

    q.push_back(vs);

  }
  else
  {
    // both vertices are on the same side
    if (side1==GEO::POSITIVE)
    {
      // both vertices are in the voronoi cell
      q.push_back( e0 );
      q.push_back( e1 );
    }
    else
    {
      // both vertices are outside the voronoi cell
    }
  }

}

} // delaunay

} // avro
