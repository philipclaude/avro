//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
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

    // optimize the rvd
    rvd.optimise( nb_iter , exact_ );

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

void
random_point_in_simplex( const Points& points , const index_t* v , index_t nv , std::vector<real_t>& p )
{
  coord_t dim = points.dim();

  // calculate random barycentric coordinates
  std::vector<real_t> alpha(nv);
  alpha[nv-1] = 1.0;
  for (coord_t j=0;j<nv-1;j++)
  {
    alpha[j]    = real_t(rand())/real_t(RAND_MAX);
    alpha[nv-1] -= alpha[j];
  }

  // use barycentric coordinates to compute p
  std::fill(p.begin(),p.end(),0.);
  for (coord_t i=0;i<nv;i++)
  {
    for (coord_t d=0;d<dim;d++)
      p[d] += alpha[i] * points[v[i]][d];
  }
}

void
CentroidalVoronoiTessellation::sample_geometry( Entity* entity , index_t nb_samples )
{
  int idx = -1;
  for (index_t k=0;k<entities_.size();k++)
  {
    if (entities_[k] == entity)
    {
      idx = k;
      break;
    }
  }
  avro_assert( idx>=0 );
  const Topology<Simplex>& topology = *topologies_[idx].get();

  coord_t np = topology.number() +1;
  index_t simplex_begin = 0;
  index_t simplex_end = topology.nb();
  coord_t dim = sites_.dim();

  avro_assert( dim == topology.points().dim() );
  avro_assert( topology.nb()>0 );

  if (np==1)
  {
    avro_assert(topology.nb()==1);
    index_t id = sites_.nb();
    sites_.create(topology.points()[ topology(0)[0] ]);
    sites_.set_entity(id,entity);
    return;
  }

  // generate random numbers for every sample
  std::vector<real_t> s(nb_samples);
  for (index_t i=0;i<nb_samples;i++)
    s[i] = real_t( double(rand())/double(RAND_MAX) );
  std::sort(s.begin(),s.end());

  // measure the volume of each element and the total volume
  std::vector<real_t> vt(topology.nb());
  real_t vtot = 0.;
  topology.get_volumes(vt);
  for (index_t k=0;k<vt.size();k++)
    vtot += vt[k];

  int first_s = -1;
  int last_s  = 0;

  index_t cur_t = simplex_begin;
  real_t cur_s = vt[0]/vtot;
  for (index_t i=0;i<nb_samples;i++)
  {
    avro_assert( i<s.size() );
    while (s[i] > cur_s && cur_t < simplex_end-1 )
    {
      cur_t++;
      avro_assert( cur_t < simplex_end );
      cur_s += vt[cur_t]/vtot;
    }
    if (first_s==-1)
      first_s = int(cur_t);
    last_s = ( last_s > int(cur_t) ) ? last_s : int(cur_t);

    // generate random point in simplex
    std::vector<real_t> cur_p(dim);
    random_point_in_simplex( topology.points() , topology(cur_t) , topology.nv(cur_t) , cur_p );

    // do something with the point
    index_t id = sites_.nb();
    sites_.create(cur_p.data());
    sites_.set_entity(id,entity);
  }
}

void
CentroidalVoronoiTessellation::generate_sites( int nb_samples )
{
  if (nb_samples<0)
  {
    // TODO calculate the total volume and make an assumption about
    // valency to compute the total number of simplices and hence points
    avro_implement;
  }

  // TODO calculate how many samples are needed on each geometry entity
  // using the valency assumptions

  for (index_t k=0;k<entities_.size();k++)
  {
    sample_geometry( entities_[k] , nb_samples );
  }
}

} // delaunay

} // avro
