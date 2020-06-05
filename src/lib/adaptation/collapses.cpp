//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "graphics/application.h"
#include "graphics/user_interface.h"
#include "graphics/window.h"

#include "library/meshb.h"

namespace avro
{

template<typename type>
Collapse<type>::Collapse( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("collapser");
  nb_accepted_.resize(_topology.number()+1);
  nb_rejected_.resize(_topology.number()+1);
}

template<typename type>
bool
Collapse<type>::visible_geometry( index_t p , index_t q , Entity* g , bool edge )
{
  avro_assert( g->number()==2 );
  if (!this->curved_) return true;
  if (g->interior()) return true;

  this->set_entity(g);
  this->geometry_cavity_.set_entity(g);

  // extract the geometry cavity
  this->extract_geometry( g, {p} );
  avro_assert( this->geometry_topology_.nb()>0 );

  // based on the type of face, we may need to flip the sign for the volume calculation
  this->geometry_cavity_.sign() = g->sign();
  if (this->topology_.master().parameter()) this->geometry_cavity_.sign() = 1;

  if (this->topology_.master().parameter())
  {
    this->convert_to_parameter(g);
  }

  // compute the cavity about the original configuration
  bool accept = this->geometry_cavity_.compute( this->v2u_[p] , this->u_[this->v2u_[p]] , this->S_ );
  avro_assert( accept );

  if (this->topology_.master().parameter())
  {
    this->convert_to_physical();
  }

  // reset the geometry checker to this face (which also computes the normals at the vertices of the geometry topology)
  this->geometry_inspector_.reset(g);
  int s = this->geometry_inspector_.signof( this->geometry_cavity_ );
  avro_assert( s > 0 );
  avro_assert( this->geometry_inspector_.positive_volumes(this->geometry_cavity_,this->geometry_cavity_.sign()));

  if (this->topology_.master().parameter())
  {
    this->convert_to_parameter(g);
  }

  // check visibility in the new configuration
  nb_parameter_tests_++;
  accept = this->geometry_cavity_.compute( this->v2u_[q] , this->u_[this->v2u_[q]] , this->S_ );
  if (!accept)
  {
    if (edge) nb_rej_vis_Edge_++;
    nb_parameter_rejections_++;
    return false;
  }

  if (this->topology_.master().parameter())
  {
    this->convert_to_physical();
  }

  // check the normal orientations
  s = this->geometry_inspector_.signof( this->geometry_cavity_ );
  if (s<0)
  {
    if (edge) nb_rej_sgn_Edge_++;
    return false;
  }

  // check the geometry topology can be extracted correctly
  if (this->geometry_inspector_.invalidates_geometry(this->geometry_cavity_) ||
      this->geometry_inspector_.invalidates_geometry(this->boundary()))
  {
    if (edge) nb_rej_geo_Edge_++;
    nb_invalid_geometry_++;
    return false;
  }

  // ensure we have all positive volumes
  avro_assert( this->geometry_inspector_.positive_volumes(this->geometry_cavity_,this->geometry_cavity_.sign()));

  return true;
}

template<typename type>
bool
Collapse<type>::valid( const index_t p , const index_t q )
{
  if (this->topology_.points().fixed(p)) return false;

  // check if a lower dimensional geometry entity is collapsed to a higher
  // dimensional one
  Entity* e0 = this->topology_.points().entity(p);
  Entity* e1 = this->topology_.points().entity(q);

  // check if the removed vertex is a volume one (always valid)
  if (e0==NULL) return true;

  if (e1==NULL)
  {
    // e0 is has a non-null geometry, so e1 cannot be a volume point
    return false;
  }

  // e0 must be higher in the geometry hierarchy than e1
  if (e0->above(e1))
  {
    if (this->geometry(p,q)==NULL) return false;
    return true;
  }

  // if the entities are not the same then we cannot collapse the edge
  // this also accounts for node-node collapses
  if (e0!=e1) return false;

  // check that the edge has a ghost
  if (this->geometry(p,q)==NULL) return false;

  if (e0!=NULL && e1!=NULL)
  {
    // check if the edge between the points contains a ghost
    bool contains = false;
    std::vector<index_t> shell;

    this->topology_.intersect( {p,q} , shell );
    for (index_t k=0;k<shell.size();k++)
    {
      if (this->topology_.ghost(shell[k]))
      {
        contains = true;
        break;
      }
    }
    Entity* g = this->geometry(p,q);
    if (!contains && !g->interior()) return false;
  }

  return true;
}

template<typename type>
bool
Collapse<type>::apply( const index_t p , const index_t q , bool delay )
{
  this->check_visibility(true);

  // attempt to collapse p onto q
  // assign the cavity
  this->topology_.intersect( {p} , this->C_ );

  this->info_ = "trying to collapse " + stringify<index_t>(p) +
                " onto " + stringify<index_t>(q);

  // count initial number of ghosts
  index_t nb_ghost0 = 0;
  for (index_t k=0;k<this->C_.size();k++)
  {
    if (this->topology_.ghost(this->C_[k]))
      nb_ghost0++;
  }

  bool accept;
  Entity* g = this->geometry(p,q);
  this->set_entity(g);

  if (this->topology_.master().parameter())
  {
    avro_assert( g!=nullptr );

    if (g->number()!=2)
    {
      avro_assert(g->number()==1);

      // loop through the parents of this geometry Edge
      for (index_t k=0;k<g->nb_parents();k++)
      {
        if (g->parents(k)->number()!=2 || !g->parents(k)->tessellatable())
          continue;

        // determine all elements on this parent Face
        this->Cavity<type>::clear();
        for (index_t j=0;j<this->C_.size();j++)
        {
          Entity* entityp = BoundaryUtils::geometryFacet( this->topology_.points() , this->topology_(this->C_[j]) , this->topology_.nv(this->C_[j]) );
          avro_assert( entityp!=nullptr );
          avro_assert( entityp->number()==2 );
          if (entityp != g->parents(k)) continue;
          this->add_cavity( this->C_[j] );
        }

        Entity* parent = g->parents(k);
        if (!visible_geometry(p,q,parent,true)) // flag this is an edge for the counters
        {
          nb_rejected_[g->number()]++;
          return false;
        }
      }

      // we need to insert the correct topology
      // so recompute the cavity (visibility check off) for all
      // original cavity elements
      this->check_visibility(false);
      bool accept = this->compute(q,this->topology_.points()[q],this->C_);
      avro_assert(accept);

      // apply the operator if no delay was requested
      if (!delay)
        this->topology_.apply(*this);

      return true;
    }
    else
    {
      this->Cavity<type>::clear();
      for (index_t j=0;j<this->C_.size();j++)
        this->add_cavity( this->C_[j] );
      if (!visible_geometry(p,q,g,false))
      {
        nb_rejected_[g->number()]++;
        return false;
      }

      // we need to insert the correct topology
      // so recompute the cavity (visibility check off) for all
      // original cavity elements
      this->check_visibility(false);
      bool accept = this->compute(q,this->topology_.points()[q],this->C_);
      avro_assert(accept);

      // apply the operator if no delay was requested
      if (!delay)
        this->topology_.apply(*this);

      return true;
    }

    return true;
  }

  // turn off enlarging
  this->enlarge_ = false;
  accept = this->compute( q , this->topology_.points()[q] , this->C_ );
  if (!accept)
  {
    // the point is not visible
    return false;
  }

  if (!this->has_unique_elems())
  {
    // we don't want duplicate elements! can happen with ghosted topologies..
    return false;
  }

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positive_implied_metrics())
    return false;

  // determine if the receiving point is visible in the parameter space
  if (g!=NULL && g->number()==2 && !visible_geometry(p,q,g))
  {
    nb_rejected_[g->number()]++;
    return false;
  }

  // determine if the point is visible on every face parenting the Edge
  if (g!=NULL && this->topology_.number()==3 && g->number()==1 && this->curved_)
  {
    // we need to check for visibility on all Faces that parent this Edge
    for (index_t k=0;k<g->nb_parents();k++)
    {
      Entity* parent = g->parents(k);
      if (parent->number()!=2 || !parent->tessellatable())
      {
        continue;
      }
      if (!visible_geometry(p,q,parent,true)) // flag this is an edge for the counters
      {
        nb_rejected_[g->number()]++;
        return false;
      }
    }
  }

  // check if only ghosts are created
  index_t only_ghost = true;
  for (index_t k=0;k<this->nb();k++)
  {
    if (!this->ghost(k))
    {
      only_ghost = false;
      break;
    }
  }
  if (only_ghost) return false;

  // apply the operator if no delay was requested
  if (!delay)
    this->topology_.apply(*this);

  if (g==NULL) nb_accepted_[this->topology_.number()]++;
  else nb_accepted_[g->number()]++;
  return true;
}

template<typename type>
void
AdaptThread<type>::collapse_edges( bool limitLength , bool swapout )
{
  index_t pass = 0;
  index_t nb_candidates,nb_collapsed,nb_collapsed_total = 0;
  index_t nb_swaps;
  real_t lmin0,lmax0,lmin1,lmax1;

  printf("-> performing collapses:\n");

  // get the current list of edges
  std::vector<index_t> edges;
  topology_.get_edges(edges);
  std::vector<real_t> lengths;

  collapser_.nb_parameter_rejections() = 0;
  collapser_.nb_parameter_tests() = 0;
  collapser_.nb_invalid_geometry() = 0;

  // evaluate the current worst quality
  // this is used as the improvement criterion when swapping out of
  // uncollapsable configuration
  real_t Q0 = worst_quality(topology_,metric_);

  while (true)
  {

    pass++;
    bool collapsed = false;
    nb_collapsed = 0;
    nb_swaps = 0;

    // compute all the lengths
    index_t ne = edges.size()/2;
    lengths.resize( ne );
    for (index_t k=0;k<ne;k++)
      lengths[k] = metric_.length( topology_.points() ,
                                  edges[2*k] , edges[2*k+1] );

    // sort the edges by length
    std::vector<index_t> idx = linspace( ne );
    std::sort( idx.begin() , idx.end() , SortBy<real_t>(lengths) );

    // compute the original length bounds
    lmin0 = * std::min_element( lengths.begin() , lengths.end() );
    lmax0 = * std::max_element( lengths.begin() , lengths.end() );

    // list the collapse candidates
    std::vector<index_t> node0;
    std::vector<index_t> node1;
    std::vector<bool> swapok;
    for (index_t k=0;k<ne;k++)
    {
      index_t i = idx[k];
      if (lengths[i]<sqrt(.5))
      {
        // add the nodes of the edge
        node0.push_back( edges[2*i] );
        node1.push_back( edges[2*i+1] );
        swapok.push_back( false );

        node0.push_back( edges[2*i+1] );
        node1.push_back( edges[2*i] );
        swapok.push_back( true );
      }
      else
        break; // the next edge is long enough that it doesn't need collapsing
    }

    nb_candidates = node0.size()/2;

    for (coord_t d=0;d<topology_.number()+1;d++)
    {
      collapser_.nb_rejected(d) = 0;
      collapser_.nb_accepted(d) = 0;
    }
    collapser_.nb_rej_vis_Edge() = 0;
    collapser_.nb_rej_sgn_Edge() = 0;
    collapser_.nb_rej_geo_Edge() = 0;

    printf("\tpass %lu: ne = %lu, short = %lu, l = [%3.4f,%3.4f]\n",
                      pass,ne,nb_candidates,lmin0,lmax0);

    // remove the nodes
    std::vector<bool> removed( topology_.points().nb() , false );
    index_t edge = 0;
    while (true)
    {
      if (node0.size()==0) break;
      if (edge==node0.size()-1) break;

      // get the next edge to collapse
      index_t n = node0[edge];
      index_t p = node1[edge];

      if (removed[n] || removed[p])
      {
        edge++;
        continue;
      }

      bool result = false;
      if (!collapser_.valid(n,p))
      {
        // try swapping out of this configuration
        if (swapout && swapok[edge])
        {
          result = swap_edge( n , p , Q0 , lmin0 , lmax0 );
          if (result) nb_swaps++;
        }

        // go to the next edge, we already added the reverse nodes so
        // those will be checked
        edge++;
        continue;
      }

      if (removed[n] || removed[p])
      {
        // nodes were removed, skip this collapse
        edge++;
        continue;
      }

      // attempt the collapse, delay application
      result = collapser_.apply( n , p , true );
      if (!result)
      {

        // try swapping out of this configuration
        if (swapout && swapok[edge])
        {
          result = swap_edge( n , p , Q0 , lmin0 , lmax0 );
          if (result) nb_swaps++;
        }

        // the operator didn't like this due to visibility
        edge++;
        continue;
      }

      // check if we want to bound the maximum length
      if (limitLength)
      {
        std::vector<real_t> lens;
        metric_.lengths( collapser_ , lens );
        if (lens.size()==0)
        {
          collapser_.print();
          avro_assert_not_reached;
        }
        real_t lmax = * std::max_element( lens.begin() , lens.end() );
        if (lmax>lmax0)
        {
          // try swapping out of this configuration
          if (swapout && swapok[edge])
          {
            result = swap_edge( n , p , Q0 , lmin0 , lmax0 );
            if (result) nb_swaps++;
          }

          edge++;
          continue;
        }
      }

      // make sure the quality does not globally degrade
      real_t qwc = worst_quality( collapser_ , metric_ );
      if (qwc<Q0)
      {
        edge++;
        continue;
      }

      // the collapse was finally accepted! apply the topology change
      topology_.apply(collapser_);

      collapsed = true;
      nb_collapsed++;

      // mark the removed nodes
      removed[n] = true;

      // go to the next edge
      edge++;

    } // loop over current edges to collapse

    nb_collapsed_total += nb_collapsed;

    // if no collapses were performed, we're done
    if (!collapsed)
    {
      printf("-> done collapses. total collapses = %lu. nb_param_rej = (%lu/%lu), nb_geom_rej = %lu.\n",
                nb_collapsed_total,collapser_.nb_parameter_rejections(),
                collapser_.nb_parameter_tests(),collapser_.nb_invalid_geometry());
      break;
    }

    // now actually delete the points
    index_t count = 0;
    for (index_t k=0;k<removed.size();k++)
    {
      if (removed[k])
      {
        topology_.remove_point( k-count );
        metric_.remove(k-count);
        topology_.inverse().remove( k-count );

        count++;
      }
    }

    // analyze the resulting edges (and recompute for the next iteration)
    edges.clear();
    topology_.get_edges(edges);

    // compute all the lengths
    ne = edges.size()/2;
    lengths.resize( ne );
    for (index_t k=0;k<ne;k++)
      lengths[k] = metric_.length( topology_.points() ,
                                  edges[2*k] , edges[2*k+1] );

    lmin1 = * std::min_element( lengths.begin() , lengths.end() );
    lmax1 = * std::max_element( lengths.begin() , lengths.end() );
    printf("\t\tcol = %lu, swap = %lu l = [%3.4f,%3.4f]\n",
                nb_collapsed,nb_swaps,lmin1,lmax1);
    if (topology_.number()>=1) printf("\t\t--> Edges:   accepted %lu, rejected %lu\n",collapser_.nb_accepted(1),collapser_.nb_rejected(1));
    if (topology_.number()>=2) printf("\t\t--> Faces:   accepted %lu, rejected %lu\n",collapser_.nb_accepted(2),collapser_.nb_rejected(2));
    if (topology_.number()>=3) printf("\t\t--> Volumes: accepted %lu, rejected %lu\n",collapser_.nb_accepted(3),collapser_.nb_rejected(3));
    if (topology_.number()>=4) printf("\t\t--> HypVols: accepted %lu, rejected %lu\n",collapser_.nb_accepted(4),collapser_.nb_rejected(4));
    printf("\t\tEdge rejections: vis = %lu, sgn = %lu, geo = %lu\n",collapser_.nb_rej_vis_Edge(),collapser_.nb_rej_sgn_Edge(),collapser_.nb_rej_geo_Edge());

  } // loop over recomputation of edges
}

template class Collapse<Simplex>;
template class AdaptThread<Simplex>;

} // avro
