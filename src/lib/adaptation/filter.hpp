// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/kdtree.h"
#include "common/tools.h"

#include "adaptation/cavity.h"
#include "adaptation/filter.h"
#include "adaptation/metric.h"
#include "adaptation/primitive.h"

#include "mesh/topology.h"

#include <egads.h>

namespace luna
{

bool
Filter::tooclose( Points& v , const index_t k ,
                  MetricField<Simplex>& metric , const real_t lmin , index_t nn )
{
  // search for nn neighbours which are permanent (i.e. skip candidates)
  index_t nu = nn;
  while (true)
  {
    // query the tree for nu nearest neighbours
    std::vector<index_t> idx(nu);
    std::vector<real_t> dist(nu);
    tree_->getNearestNeighbours( operator[](k) , idx , dist , nu );

    // loop through the returned indices and count number of permanents
    index_t np = 0;
    for (index_t j=0;j<nu;j++)
    {
      // other candidates are not considere yet
      // in any case, the other candidates do not have a metric evaluated
      // them so the metric cannot be evaluated
      if (!permanent( idx[j] )) continue;
      np++;

      // if the point is permanent, it has an associated index in the provided
      // points (or at least it should)
      index_t pt = idx_[idx[j]];

      luna_assert( pt < v.nb() );
      luna_assert( k < v.nb() );

      // check the distance, if it's too short then we can obviously skip the rest
      real_t length = metric.length( v , k , pt );
      if (length<lmin)
      {
        return true;
      }
    }

    // increase the number of neighbours searched for the next iteration
    nu = 2*nu;
    if (nu>nb()) nu = nb();

    if (np>=nn) break;
    if (nu>=nb()) break;
  }

  return false;
}

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


template<typename Metric>
void
Filter::generateCandidates( Topology<Simplex>& topology ,
                            Metric& metric , real_t lmax , Insert<Simplex>& inserter )
{
  // get all the edges
  std::vector<index_t> edges;
  topology.getEdges( edges );

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
  real_t params[4] = {0,0,0,0};
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
    real_t ds = 1./real(ns);

    // create candidates that will subdivide the edge
    for (index_t j=1;j<ns;j++)
    {
      for (coord_t d=0;d<dim_;d++)
        xs[d] = x0[d] +j*ds*( x1[d] -x0[d] );

      // determine if the coordinates need to be projected
      Entity* eg = inserter.geometry(n0,n1);
      if (eg!=NULL && inserter.curved())
      {
        if (eg->tesseractGeometry())
          eg->project( xs );
        else
        {
          v[0] = n0;
          v[1] = n1;
          geometryParams( eg , topology.points() , v , 2 , params );

          // interpolate in parameter space
          params[0] = params[0] +j*ds*( params[udim  ] -params[0] );
          params[1] = params[1] +j*ds*( params[udim+1] -params[1] );

          // evaluate on the geometry
          eg->evaluate( params, xs.data() );
        }
      }

      // create the insertion candidate
      createCandidate( n0 , n1 , j*ds , xs.data() , eg , params );

      // reset parameters to 'none'
      params[0] = Points::INFTY;
      params[1] = Points::INFTY;
    }
  }
}

} // local

} // avro
