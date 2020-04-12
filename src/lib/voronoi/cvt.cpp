#include "geometry/entity.h"

#include "mesh/boundary.h"

#include "voronoi/cvt.h"
#include "voronoi/voronoi.h"

namespace avro
{

namespace delaunay
{

CentroidalVoronoiTessellation::CentroidalVoronoiTessellation( const Topology<Simplex>& topology , Points& sites , bool hierarchical ) :
  Topology<Polytope>(points_,topology.number()),
  points_(sites.dim()),
  sites_(sites),
  hierarchical_(hierarchical),
  exact_(true)
{
  if (hierarchical_)
  {
    Boundary<Simplex> boundary( topology );
    boundary.extractall();

    // create an rvd for each boundary topology
    for (index_t k=0;k<boundary.nb();k++)
    {
      entities_.push_back( boundary.entity(k) );
      std::shared_ptr<Topology<Simplex>> bnd_topology = std::make_shared<Topology<Simplex>>(topology.points(),boundary.child(k).number() );
      topologies_.push_back( bnd_topology );
    }
  }

  std::shared_ptr<Topology<Simplex>> main_topology = std::make_shared<Topology<Simplex>>(topology.points(),topology.number() );
  main_topology->TopologyBase::copy(topology);
  topologies_.push_back( main_topology );
  entities_.push_back(nullptr);

}

void
CentroidalVoronoiTessellation::compute( index_t nb_iter )
{
  avro_assert( topologies_.size() == entities_.size() );

  // optimize the lower-dimensional rvd's first
  for (index_t k=0;k<topologies_.size();k++)
  {

    // get the geometry entity this rvd corresponds to
    Entity* entity = entities_[k];

    // extract the topology to use as the background
    const Topology<Simplex>& topology = *topologies_[k].get();

    // extract the sites to be included inthe computatation of this rvd
    Delaunay z( sites_.dim() );
    for (index_t j=0;j<sites_.nb();j++)
    {
      Entity* ej = sites_.entity(j);
      if (entity==NULL)
      {
        z.create( sites_[j] );
        continue;
      }

      if (ej==NULL) continue; // entity canot be nonnull with ek null

      if (entity->above(ej))
        z.create(sites_[j]);
    }

    printf("nb sites on entity %lu = %lu\n",k,z.nb());

    // create the restricted voronoi diagram structure
    RestrictedVoronoiDiagram rvd( topology , z ); // TODO add entity information to know which points to keep fixed

    rvd.parallel() = true;

    // compute the rvd
    rvd.compute(exact_);

    // optimize this rvd
    for (index_t iter=0;iter<nb_iter;iter++)
    {
      // TODO move these sites closer to the centroids

      // TODO track the cvt energy
    }

    // add the rvd to this topology
    index_t offset = points_.nb();
    for (index_t j=0;j<rvd.points().nb();j++)
    {
      points_.create( rvd.points()[j] );
      std::vector<int> facets = rvd.points().incidence().get(j);
      points_.incidence().add( facets.data() , facets.size() );
    }

    if (k==topologies_.size()-1)
    {
      for (index_t j=0;j<rvd.nb();j++)
      {
        std::vector<index_t> pj = rvd.get(j);
        for (index_t i=0;i<pj.size();i++)
          pj[i] += offset;
        this->add(pj.data(),pj.size());
      }

      std::shared_ptr<VoronoiSites> s = std::make_shared<VoronoiSites>(*this);
      s->build();

      for (index_t k=0;k<rvd.sites().nb_data();k++)
        s->value(k) = rvd.sites().value(k);

      sites_fields_.push_back(s);
      fields_.make("sites",s);
    }

  }

}

} // delaunay

} // avro
