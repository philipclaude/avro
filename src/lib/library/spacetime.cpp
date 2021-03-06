//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/entity.h"
#include "geometry/psc/facet.h"

#include "library/meshb.h"
#include "library/spacetime.h"

#include "mesh/boundary.h"

#include "numerics/linear_algebra.h"

namespace avro
{

template<typename type>
Topology_Spacetime<type>::Topology_Spacetime( const Topology<type>& topology ) :
  Topology<type>(points_,3),
  points_(3),
  topology_(topology)
{
  // setup the child topologies
  // first retrieve any possible boundary entity
  std::set<Entity*> entities;
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    Entity* entity = topology_.points().entity(k);
    if (entity==nullptr) continue;
    if (entity->number()!=3) continue;
    entities.insert( entity );
  }

  std::set<Entity*>::iterator it;
  for (it=entities.begin();it!=entities.end();it++)
  {
    std::shared_ptr<Topology<type>> t = std::make_shared<Topology<type>>(points_,3);
    entity2topology_.insert( {*it,t.get()} );
    topology2entity_.insert( {t.get(),*it} );
    this->add_child(t);
  }
}

template<typename type>
void
Topology_Spacetime<type>::extract()
{
  clear();
  Boundary<type> boundary(topology_);
  boundary.extract();

  const coord_t dim = topology_.points().dim();
  avro_assert( dim==4 ); // spacetime mesh but just checking

  // loop through the entities on the boundary
  const std::vector<Entity*>& entities = boundary.entities();

  for (index_t k=0;k<entities.size();k++)
  {
    Entity* entity = entities[k];
    Topology<type>& tk = *entity2topology_.at(entity);
    Topology<type>& bk = boundary.child( boundary.indexof(entity) );

    //const matd<real_t>& B = static_cast<PSC::Facet*>(entity)->basis();
    //const vecd<real_t>& x0 = static_cast<PSC::Facet*>(entity)->x0();

    // create a vertex for every vertex on the boundary
    std::set<index_t> vertices;
    for (index_t k=0;k<bk.nb();k++)
    {
      for (index_t j=0;j<bk.nv(k);j++)
        vertices.insert( bk(k,j) );
    }

    std::set<index_t>::iterator it;
    std::map<index_t,index_t> identifier;
    index_t n0 = points_.nb();
    for (it=vertices.begin();it!=vertices.end();it++)
    {
      // create a vertex for this one
      index_t idx = points_.nb();

      Entity* ep = topology_.points().entity(*it);
      avro_assert( ep==entity || entity->above(ep) );

      std::vector<real_t> X( topology_.points()[*it] , topology_.points()[*it]+dim );
      std::vector<real_t> U( 3 );
      entity->inverse( X , U );

      points_.create( U.data() );
      points_.set_entity( idx , entity );
      identifier.insert( {*it,idx} );
    }

    #if 1
    std::vector<real_t> xmin(3, 1e20);
    std::vector<real_t> xmax(3,-1e20);

    for (index_t j = n0; j < points_.nb(); j++) {
      for (coord_t d = 0; d < 3; d++) {
        real_t x = points_(j,d);
        if (x < xmin[d]) xmin[d] = x;
        if (x > xmax[d]) xmax[d] = x;
      }
    }
    //printf("xmin = %g, xmax = %g\n",xmin,xmax);

    for (index_t j = n0; j < points_.nb(); j++) {
      for (coord_t d = 0; d < 3; d++) {
        points_(j,d) = (points_(j,d) - xmin[d]) / (xmax[d] - xmin[d]);
      }
    }
    #endif

    // create the elements in the topology
    std::vector<index_t> simplex(4);
    for (index_t k=0;k<bk.nb();k++)
    {
      for (index_t j=0;j<bk.nv(k);j++)
        simplex[j] = identifier[ bk(k,j) ];
      tk.add( simplex.data() , simplex.size() );
    }
    tk.orient();
  }


}

template<typename type>
void
Topology_Spacetime<type>::write( const std::string& filename )
{
  library::meshb writer;

  writer.open( 3 , filename );
  writer.write( points_ );

  Topology<type> topology( points_ , 3 );
  std::vector<index_t> refs;
  for (index_t j=0;j<this->nb_children();j++)
  {
    const Topology<type>& tk = this->child(j);
    for (index_t k=0;k<tk.nb();k++)
    {
      topology.add( tk(k) , tk.nv(k) );
      refs.push_back( j );
    }
  }
  writer.write( topology , refs );

  writer.close();
}

template<typename type>
void
Topology_Spacetime<type>::clear()
{
  points_.clear();
  for (index_t k=0;k<this->nb_children();k++)
      this->child(k).clear();
}

template class Topology_Spacetime<Simplex>;

} // avro
