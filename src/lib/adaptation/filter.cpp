//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/filter.h"
#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include "geometry/entity.h"

#include "mesh/topology.h"

#include "common/kdtree.h"

#include <algorithm>

namespace avro
{

static inline void
prioritizeGeometry( std::vector<index_t>& edges0 , const Points& points , Insert<Simplex>& prim )
{
  index_t ne = edges0.size()/2;
  std::vector<bool> traversed( ne , false );

  std::vector<index_t> edges;
  for (index_t k=0;k<ne;k++)
  {
    // check if this is a geometry edge
    if (prim.geometry(edges0[2*k],edges0[2*k+1]))
      traversed[k] = true;
    else
      continue;
    edges.push_back( edges0[2*k] );
    edges.push_back( edges0[2*k+1] );
  }

  // now add the remaining edges
  for (index_t k=0;k<ne;k++)
  {
    // check if this is a geometry edge
    if (traversed[k])
      continue;
    edges.push_back( edges0[2*k] );
    edges.push_back( edges0[2*k+1] );
  }

  edges0 = edges;
}

static inline void
prioritizeLongest( std::vector<index_t>& edges0 , const Points& points , std::vector<real_t>& lengths0 , Insert<Simplex>& prim )
{
  index_t ne = edges0.size()/2;
	std::vector<bool> traversed( ne , false );

  // sort by length, the longest are at the end
  std::vector<index_t> idx = linspace( ne );
  std::sort( idx.begin() , idx.end() , SortBy<real_t>(lengths0) );

	// add the geometry edges first
  std::vector<index_t> edges;
	std::vector<real_t> lengths;
  for (index_t k=0;k<ne;k++)
	{
		// longest edge
		index_t e = idx[ne-k-1];

		// we want to keep geometry edges first
		if (!prim.geometry( edges0[2*e] , edges0[2*e+1] ))
			continue;

		traversed[e] = true;

		edges.push_back( edges0[2*e] );
		edges.push_back( edges0[2*e+1] );
		lengths.push_back( lengths0[e] );
	}

	// go back through the remaining edges
	for (index_t k=0;k<ne;k++)
	{
		index_t e = idx[ne-k-1];
		if (traversed[e])
			continue;

		edges.push_back( edges0[2*e] );
		edges.push_back( edges0[2*e+1] );
		lengths.push_back( lengths0[e] );
	}
	edges0 = edges;
	lengths0 = lengths;
}

Filter::Filter( const coord_t dim ) :
  Points(dim),
  lmin_(0.),
  lmax_(0.)
{
}

void
Filter::createPermanent( const real_t* x )
{
  node0_.push_back(-1);
  node1_.push_back(-1);
  s_.push_back( -1. );
  idx_.push_back( nb() );
  Points::create(x);
}

void
Filter::createCandidate( const index_t n0 , const index_t n1 ,
                         const real_t s , const real_t* x,
                         Entity * e, const real_t* params )
{
  avro_assert( n0 < n1 );
  node0_.push_back(n0);
  node1_.push_back(n1);
  s_.push_back(s);
  candidates_.push_back(nb());
  idx_.push_back(-1);
  Points::create(x);
  Points::set_entity(nb()-1,e);
  Points::set_param(nb()-1,params);
}

void
Filter::generateCandidates( Topology<Simplex>& topology ,
                            MetricField<Simplex>& metric , real_t lmax , Insert<Simplex>& inserter )
{
  // get all the edges
  std::vector<index_t> edges;
  topology.get_edges( edges );

	// prioritize geometry edges
  prioritizeGeometry(edges,topology.points(),inserter);

  nb_long_ = 0;

  index_t ne = edges.size()/2;

  std::vector<real_t> lengths( ne , -1. );
  for (index_t k=0;k<ne;k++)
  {
    index_t n0 = edges[2*k];
    index_t n1 = edges[2*k+1];

    // do not create candidates on edges touching a ghost
    if (n0<topology.points().nb_ghost() || n1<topology.points().nb_ghost() )
      continue;

    lengths[k] = metric.length(topology.points(),n0,n1);
  }

	// prioritize by longest length
	prioritizeLongest( edges , topology.points() , lengths , inserter );

  real_t Lmax = *std::max_element(lengths.begin(),lengths.end());

  lmin_ = * std::min_element( lengths.begin() , lengths.end() );
  lmax_ = Lmax;

  const coord_t udim = topology.points().udim();
  index_t v[2];
  std::vector<real_t> params(4,0);
  std::vector<real_t> xs(std::max(coord_t(3),dim_)); // EGADS evaluations are always 3

  // loop through all edges and create candidates along long edges
  for (index_t k=0;k<ne;k++)
  {
    index_t n0 = edges[2*k];
    index_t n1 = edges[2*k+1];

    // do not create candidates on edges touching a ghost
    if (n0<topology.points().nb_ghost() || n1<topology.points().nb_ghost() )
      continue;

    real_t* x0 = topology.points()[n0];
    real_t* x1 = topology.points()[n1];

    real_t length = lengths[k];
    if (length<lmax) continue;
    nb_long_++;

    // divide the edge
    index_t ns;
    if (Lmax<2)
      ns = ceil(length);
    else
      ns = floor(length);
    ns = 2.;
    real_t ds = 1./real_t(ns);

    // create candidates that will subdivide the edge
    for (index_t j=1;j<ns;j++)
    {
      for (coord_t d=0;d<dim_;d++)
        xs[d] = x0[d] +j*ds*( x1[d] -x0[d] );

      // determine if the coordinates need to be projected
      Entity* eg = inserter.geometry(n0,n1);
      if (eg!=NULL && inserter.curved())
      {
        if (eg->egads())
        //if (eg->tesseractGeometry())
        //  eg->project( xs );
        //else
        {
          v[0] = n0;
          v[1] = n1;
          geometry_params( eg , topology.points() , v , 2 , params.data() );

          // interpolate in parameter space
          params[0] = params[0] +j*ds*( params[udim  ] -params[0] );
          params[1] = params[1] +j*ds*( params[udim+1] -params[1] );

          // evaluate on the geometry
          eg->evaluate( params, xs );
        }
      }

      // create the insertion candidate
      createCandidate( n0 , n1 , j*ds , xs.data() , eg , params.data() );

      // reset parameters to 'none'
      params[0] = Points::INFTY;
      params[1] = Points::INFTY;
    }
  }
}

void
Filter::buildTree()
{
  // setup the kdtree
  cloud_ = std::make_shared<PointCloud>(*this);
  tree_ = initializeKdTree(*cloud_);
  tree_->build();
}

void
Filter::accept( const index_t k , const index_t idx )
{
  avro_assert_msg( node0_[k]>=0 && node1_[k]>=0 ,
      "node %lu is already active",k);

  // save the original values
  int n0 = node0_[k];
  int n1 = node1_[k];

  // look for any other candidate which touches this edge
  for (index_t i=0;i<nb();i++)
  {
    // permanent points do not get modified
    if (permanent(i)) continue;

    if (node0_[i]!=n0) continue;
    if (node1_[i]!=n1) continue;

    // check the s value for which node gets modified
    if (s_[i]<s_[k])
    {
      // node1 gets modified such that the edges is now [node0,idx]
      node1_[i] = idx;
    }
    else
    {
      // node1 gets modified such that the edge is now [idx,node1]
      node0_[i] = idx;
    }
  }

  // turn off the edges, which informs us that the node is active
  node0_[k] = -1;
  node1_[k] = -1;
  s_[k] = -1;
  idx_[k] = idx;

}

void
Filter::clearCandidates()
{
  for (int k=nb()-1;k>=0;k--)
  {
    if (permanent(k)) continue;

    node0_.erase( node0_.begin() +k );
    node1_.erase( node1_.begin() +k );
    s_.erase( s_.begin() +k );
    Points::remove( k );
  }
  candidates_.clear();
}

void
Filter::print() const
{
  for (index_t k=0;k<nb_candidates();k++)
  {
    index_t idx = candidate(k);
    printf("candidate[%lu] = %lu, node0 = %d, node1 = %d, s = %g\n",
                              k,idx,node0_[idx],node1_[idx],s_[idx]);
  }
}

} // avro
