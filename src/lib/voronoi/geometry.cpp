#include "geometry/entity.h"

#include "mesh/boundary.h"

#include "voronoi/delaunay.h"
#include "voronoi/geometry.h"
#include "voronoi/voronoi.h"

namespace avro
{

GeometryConformingRVD::GeometryConformingRVD( const Topology<Simplex>& topology , Delaunay& sites ) :
  Topology<Polytope>( vertices_ , topology.number() ),
  topology_(topology),
  sites_(sites),
  vertices_( sites.dim() ),
  triangulation_( sites , topology.number() )
{
  initialize();
}

void
GeometryConformingRVD::initialize()
{
  boundary_ = std::make_shared<Boundary<Simplex>>( topology_ );
  boundary_->extract();
  boundary_->print();

  // create an rvd and rdt for all geometry entities
  entity_.assign( boundary_->entities().begin() , boundary_->entities().end() );
  uniquify(entity_);

  // check if we have a null entity on interior vertices
  for (index_t k=0;k<sites_.nb();k++)
  {
    if (k<sites_.nb_ghost()) continue;
    if (sites_.entity(k)==nullptr)
    {
      entity_.push_back(nullptr);
      break;
    }
  }

  printf("total number of entities = %lu\n",entity_.size());
  for (index_t k=0;k<entity_.size();k++)
  {
    std::shared_ptr< delaunay::RestrictedVoronoiDiagram > rvdk;
    std::shared_ptr< RestrictedDelaunayTriangulation > rdtk;
    if (entity_[k]==nullptr)
    {
      printf("adding interior!\n");
      rvdk = std::make_shared< delaunay::RestrictedVoronoiDiagram >(topology_,sites_);
      rdtk = std::make_shared< RestrictedDelaunayTriangulation >(sites_,topology_.number());
    }
    else
    {
      entity_[k]->print_header();
      Topology<Simplex>& tk = boundary_->child( boundary_->indexof(entity_[k]) );
      avro_assert_msg( entity_[k]->number() == tk.number() , "tk number = %u" , tk.number() );
      rvdk = std::make_shared< delaunay::RestrictedVoronoiDiagram >(tk,sites_,entity_[k]);
      rdtk = std::make_shared< RestrictedDelaunayTriangulation >(sites_, tk.number() );
    }

    // create the rvd and add it as a child to this topology
    rvd_.push_back(rvdk);
    add_child(rvdk);

    // create the rdt and add it as a child to the main triangulation topology
    rdt_.push_back(rdtk);
    triangulation_.add_child(rdtk);
  }

  printf("total number of rvd's = %lu\n",rvd_.size());
}

void
GeometryConformingRVD::compute()
{
  for (index_t k=0;k<rvd_.size();k++)
  {
    rvd_[k]->compute();

    // copy the site for the interior topology
    if (entity_[k]==nullptr)
    {
      this->fields().make( "sites" , rvd_[k]->sites_ptr() );
    }
  }
}

void
GeometryConformingRVD::extract_triangulations()
{
  for (index_t k=0;k<rvd_.size();k++)
  {
    rvd_[k]->extract(*rdt_[k].get());
    printf("rdt on entity %lu has %lu simplices\n",k,rdt_[k]->nb());
  }
}

} // avro
