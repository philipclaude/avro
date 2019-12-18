// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/sort.h"
#include "common/types.h"

#include "mesh/detective.h"

#include "mesh/local/filter.h"

#include <unordered_set>

#include "mesh/boundary.h"
#include "library/mesh.h"
#include "graphics/plotter.h"

namespace avro
{

namespace local
{

template<typename type>
template<typename field>
bool
SerialAdaptation<type>::swapedge( numerics::Metric<field>& metric , index_t p , index_t q , real Q0 , real lmin0 , real lmax0 )
{
  // try edge swapping around the edge (p,q) such that the
  // global worst quality does not get worse
  std::vector<index_t> shell;
  topology_.intersect( {p,q} , shell );
  if (shell.size()==0) return false;

  // list all the vertices in the shell
  std::vector<index_t> candidates;
  for (index_t j=0;j<shell.size();j++)
  for (index_t i=0;i<topology_.nv(shell[j]);i++)
    candidates.push_back(topology_(shell[j],i));
  uniquify(candidates);

  // assign the cavity used by the swapper since this is the same
  // for all candidates
  edge_swapper_.setCavity( shell );

  // loop through all candidates
  int m = -1;
  real qw = Q0;
  for (index_t j=0;j<candidates.size();j++)
  {

    if (candidates[j]==p || candidates[j]==q) continue;

    // check if the swap is valid in terms of geometry topology
    bool accept = edge_swapper_.valid( candidates[j] , p , q );
    if (!accept)
    {
      continue;
    }

    // check if the swap is valid in terms of visibility
    accept = edge_swapper_.apply( candidates[j] , p , q );
    if (!accept)
    {
      continue;
    }

    // option to check the produced lengths
    std::vector<real> lens;
    metric.lengths( edge_swapper_ , lens );
    if (lens.size()==0) return false;
    real lmin = * std::min_element( lens.begin() , lens.end() );
    real lmax = * std::max_element( lens.begin() , lens.end() );
    if (lmin0>0 && lmin<lmin0) continue;
    if (lmax0>0 && lmax>lmax0) continue;

    edge_swapper_.evaluate(metric);
    if (edge_swapper_.worst()>qw)
    {
      if (!edge_swapper_.philipcondition()) continue;
      m  = candidates[j];
      qw = edge_swapper_.worst();
    }

  }

  if (m<0) return false;

  edge_swapper_.apply( index_t(m) , p , q );
  edge_swapper_.evaluate(metric);
  topology_.apply(edge_swapper_);
  topology_.update(edge_swapper_,metric);

  return true;

}

template<typename type>
template<typename field>
void
SerialAdaptation<type>::collapses( numerics::Metric<field>& metric , bool limitLength , bool swapout )
{
  index_t pass = 0;
  index_t nb_candidates,nb_collapsed,nb_collapsed_total = 0;
  index_t nb_swaps;
  real lmin0,lmax0,lmin1,lmax1;

  printf("-> performing collapses:\n");

  // get the current list of edges
  std::vector<index_t> edges;
  topology_.getEdges(edges);
  std::vector<real> lengths;

  collapser_.nb_parameter_rejections() = 0;
  collapser_.nb_parameter_tests() = 0;
  collapser_.nb_invalid_geometry() = 0;

  // evaluate the current worst quality
  // this is used as the improvement criterion when swapping out of
  // uncollapsable configuration
  topology_.evaluate(metric);
  real Q0 = topology_.worst();

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
      lengths[k] = metric.length( topology_.vertices() ,
                                  edges[2*k] , edges[2*k+1] );

    // sort the edges by length
    std::vector<index_t> idx = linspace( ne );
    std::sort( idx.begin() , idx.end() , SortBy<real>(lengths) );

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
    std::vector<bool> removed( topology_.vertices().nb() , false );
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
          result = swapedge( metric , n , p , Q0 , lmin0 , lmax0 );
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
          result = swapedge( metric , n , p , Q0 , lmin0 , lmax0 );
          if (result) nb_swaps++;
        }

        // the operator didn't like this due to visibility
        edge++;
        continue;
      }

      // check if we want to bound the maximum length
      if (limitLength)
      {
        std::vector<real> lens;
        metric.lengths( collapser_ , lens );
        if (lens.size()==0)
        {
          collapser_.print();
          avro_assert_not_reached;
        }
        real lmax = * std::max_element( lens.begin() , lens.end() );
        if (lmax>lmax0)
        {
          // try swapping out of this configuration
          if (swapout && swapok[edge])
          {
            result = swapedge( metric , n , p , Q0 , lmin0 , lmax0 );
            if (result) nb_swaps++;
          }

          edge++;
          continue;
        }
      }

      // make sure the quality does not globally degrade
      collapser_.evaluate(metric);
      if (collapser_.worst()<Q0)
      {
        edge++;
        continue;
      }

      // the collapse was finally accepted! apply the topology change
      topology_.apply(collapser_);

      // the quality needs to be updated if we are allowing swaps
      if (swapout)
      {
        collapser_.evaluate(metric);
        topology_.update(collapser_,metric);
      }

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

    // now actually delete the vertices
    index_t count = 0;
    for (index_t k=0;k<removed.size();k++)
    {
      if (removed[k])
      {
        topology_.removeVertex( k-count );
        metric.remove(k-count);
        topology_.inverse().remove( k-count );

        count++;
      }
    }

    // analyze the resulting edges (and recompute for the next iteration)
    edges.clear();
    topology_.getEdges(edges);

    // compute all the lengths
    ne = edges.size()/2;
    lengths.resize( ne );
    for (index_t k=0;k<ne;k++)
      lengths[k] = metric.length( topology_.vertices() ,
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

    //check("collapse pass "+std::to_string(pass));

  } // loop over recomputation of edges
}

template<typename type>
template<typename field>
void
SerialAdaptation<type>::splits( numerics::Metric<field>& metric , real lt, bool limitlength , bool swapout )
{
  avro_assert( metric.check(topology_) );

  real dof_factor = params_.insertion_volume_factor();

  index_t nb_swaps;
  index_t nb_inserted,nb_inserted_total = 0;
  index_t nb_length_rejected = 0;
  index_t nb_quality_rejected = 0;
  index_t nb_visiblity_rejected = 0;
  index_t nb_count_rejected = 0;
  real lmin1,lmax1;

  inserter_.nb_parameter_rejections() = 0;
  inserter_.nb_parameter_tests() = 0;

  topology_.evaluate(metric);
  real Q0 = topology_.worst();

  // don't be too restrictive with insertions when the quality is already good
  if (Q0>0.4) Q0 = 0.1;

  std::vector<real> lens;
  metric.lengths(topology_,lens);
  real lmin0 = *std::min_element( lens.begin() , lens.end() );
  real lmax0 = *std::max_element( lens.begin() , lens.end() );

  index_t pass = 0;
  bool done = false;
  printf("-> performing edge splits on edges with lt > %1.3f and dof_factor %g:\n",lt,dof_factor);
  while (!done)
  {
    // anything beyond 20 passes is borderline ridiculous
    // the metric is probably way off from the current mesh
    if (pass>20)
    {
      printf("warning: too many insertions, metric too far from current mesh.\n");
      break;
    }

    nb_inserted = 0;
    nb_swaps = 0;
    nb_length_rejected = 0;
    nb_quality_rejected = 0;
    nb_visiblity_rejected = 0;
    nb_count_rejected = 0;

    done = true;

    // setup the insertion filter
    Filter filter( topology_.vertices().dim() );

    // add the current vertices in the topology
    for (index_t k=0;k<topology_.vertices().nb();k++)
      filter.createPermanent( topology_.vertices()[k] );

    // create the candidates on edges longer than lt in the target space
    filter.generateCandidates( topology_ , metric , lt , inserter_ );

    printf("\t pass %lu: long = %lu, l = [%3.4f,%3.4f] -> insert %lu\n",
                pass,filter.nb_long(),filter.minlength(),filter.maxlength(),
                filter.nb_candidates());

    // build the kd-tree search structure only if we need to check the
    // distance to the nearest neighbours
    filter.buildTree();

    // insert the points
    std::unordered_set<index_t> removed;
    std::unordered_set<index_t> flagged;
    std::vector<index_t> shell;
    std::vector<index_t> N;
    for (index_t k=0;k<filter.nb_candidates();k++)
    {
      // index of the vertex stored in the filter
      index_t idx = filter.candidate(k);

      // get the edge along which this insertion takes place
      index_t n0 = filter.edge( idx , 0 );
      index_t n1 = filter.edge( idx , 1 );

      real lk = metric.length( topology_.vertices() , n0 , n1 );

      // insertions on the edges with fixed nodes are not allowed
      // as these are partition boundaries
      if (topology_.vertices().fixed(n0) && topology_.vertices().fixed(n1))
        continue;

      // do not insert on ghost edges
      if (n0<topology_.vertices().nb_ghost() ||
          n1<topology_.vertices().nb_ghost())
        continue;

      // check if the vertices were removed or flagged
      if (removed.find(n0)!=removed.end() || removed.find(n1)!=removed.end())
        continue;

      if (flagged.find(n0)!=flagged.end() || flagged.find(n1)!=flagged.end())
        continue;

      // the metric needs to be interpolated for the filter to evaluate lengths
      // so we need to add the vertex and also add the interpolated metric
      index_t ns = topology_.vertices().nb();
      topology_.vertices().create(filter[idx]);
      metric.add(n0,n1,filter[idx]);

      // notify the inverse topology that we want to store extra data
      topology_.inverse().create(1);

      // check the distance to nearby vertices
      shell.clear();
      topology_.intersect( {n0,n1} , shell );
      N.clear();
      for (index_t j=0;j<shell.size();j++)
      {
        for (index_t i=0;i<topology_.nv(shell[j]);i++)
          N.push_back( topology_(shell[j],i) );
      }
      uniquify(N);
      bool bad = false;

      // set the minimum length any insertion can create
      real Lmin = sqrt(0.5);

      // if the current length is greater than 4.0, we need to be more flexible
      if (lk>4.0) Lmin = 0.0;

      // also relax the insertion criterion when we insert on geometry Edges
      Entity* ge = inserter_.geometry(n0,n1);
      if (ge!=NULL && ge->number()==1)
        Lmin = 0.25;

      for (index_t j=0;j<N.size();j++)
      {
        real lj = metric.length( topology_.vertices() , N[j] , ns );
        if (lj<Lmin)
        {
          bad = true;
          break;
        }
      }

      bool swapped;
      if (bad && limitlength)
      {
        // remove the interpolated tensor with its associated vertex
        topology_.vertices().remove(ns);
        metric.remove(ns);
        topology_.inverse().remove(ns);
        nb_length_rejected++;

        // option to swap out of the rejected configuration
        if (swapout)
        {
          swapped = swapedge(metric,n0,n1,Q0,lmin0,lmax0);
          if (swapped) nb_swaps++;
        }
        continue;
      }
      avro_assert( metric.check(topology_) );

      // the vertex needs to be removed because the inserter will add it again
      // TODO clean up this inefficiency
      topology_.vertices().remove(ns);

      // apply the insertion
      inserter_.clear();
      inserter_.restart();
      inserter_.delay() = true;
      bool result = inserter_.apply( n0 , n1 , filter[idx] , filter.u(idx) );
      if (!result)
      {
        // the metric needs to be removed because the vertex was rejected.
        // for curverd geometries, the cavity might need to be enlarged.
        // the vertex was already removed (above)
        // so only the metric requires removing
        metric.remove(ns);
        nb_visiblity_rejected++;
        continue;
      }

      // if the inserter was enlarged, don't be too restrictive with quality
      inserter_.evaluate( metric );
      if (inserter_.worst()<Q0)
      {
        topology_.vertices().remove(ns);
        metric.remove(ns);
        topology_.inverse().remove(ns);
        nb_quality_rejected++;
        continue;
      }

      // check if the metric volume is respected
      real vol = 0.0;
      for (index_t k=0;k<inserter_.nb();k++)
        vol += metric.volume( inserter_ , k );
      index_t count = index_t(vol/topology_.master().unitVolume());
      if (dof_factor>0 && dof_factor*count<inserter_.nb_real())
      {
        if (ge==NULL || ge->number()>2)
        {
          // it's going to be really hard to come back from something like this
          topology_.vertices().remove(ns);
          metric.remove(ns);
          topology_.inverse().remove(ns);
          nb_count_rejected++;
          continue;
        }
      }

      // apply the insertion into the topology
      topology_.apply(inserter_);
      avro_assert( metric.check(topology_) );

      if (swapout)
      {
        // if we're doing swaps then we need to update the quality information
        inserter_.evaluate(metric);
        topology_.update( inserter_ , metric );
      }

      // determine if any vertices were removed
      for (index_t j=0;j<inserter_.nb_removedNodes();j++)
      {
        index_t removed_node = inserter_.removedNode(j);
        if (removed.find(removed_node)==removed.end())
        {
          printf("vertex %lu was removed!\n",removed_node);
          removed.insert( removed_node );
        }
      }

      // check if the cavity was enlarged to turn off some of the existing edges
      if (inserter_.enlarged())
      {
        // turn off vertices in the cavity, we can attempt them on the next pass
        const std::vector<index_t>& nodes = inserter_.nodes();
        for (index_t j=0;j<nodes.size();j++)
          flagged.insert( nodes[j] );
      }

      // the filter is happy
      filter.accept( idx , ns );
      done = false; // since an insertion was performed

      nb_inserted++;
    }

    // delete the vertices that were removed during cavity enlargements
    std::vector<index_t> removed_indices;
    for (std::unordered_set<index_t>::iterator it=removed.begin();it!=removed.end();it++)
      removed_indices.push_back( *it );

    std::sort( removed_indices.begin() , removed_indices.end() );
    std::reverse( removed_indices.begin() , removed_indices.end() );
    for (index_t j=0;j<removed_indices.size();j++)
    {
      topology_.removeVertex( j );
      metric.remove(j);
      topology_.inverse().remove(j);
    }

    printf("\t\tinserted %lu, swapped %lu, lrej = %lu, qrej = %lu, vrej = %lu, prej = %lu/%lu, dof_rej = %lu\n",
                nb_inserted,nb_swaps,nb_length_rejected,nb_quality_rejected,nb_visiblity_rejected,
                inserter_.nb_parameter_rejections(),inserter_.nb_parameter_tests(),nb_count_rejected);
    nb_inserted_total += nb_inserted;

    pass++;
  }

  // analyze the resulting edge lengths
  std::vector<index_t> edges;
  topology_.getEdges(edges);
  std::vector<real> lengths(edges.size()/2);
  for (index_t k=0;k<edges.size()/2;k++)
    lengths[k] = metric.length( topology_.vertices() , edges[2*k] , edges[2*k+1] );
  lmin1 = *std::min_element( lengths.begin() , lengths.end() );
  lmax1 = *std::max_element( lengths.begin() , lengths.end() );
  printf("done insertions: total insert = %lu, l = [%3.4g,%3.4g]\n",
    nb_inserted_total,lmin1,lmax1);
  printf("\tnb_parameter_rejections = (%lu/%lu)\n",
            inserter_.nb_parameter_rejections(),inserter_.nb_parameter_tests());

}

template<typename type>
template<typename field>
void
SerialAdaptation<type>::unitise( numerics::Metric<field>& metric )
{
  // collapse all short edges
  // find all edges in the mesh shorter that 1/sqrt(2)
  // accumulate the list of nodes on these edges
  // the cavity for each node is the ball of the node
  // the star_ will be any node in the ball of the node
  // perform the operation if that short edge is purged
  collapses(metric);

  // split all long edges
  // first propose candidates by subdividing the long edges
  // then check the distances are not too short by checking a subset
  // of nearest neighbours of the inserted points in the actual mesh
  splits(metric,2.0,false);

}

} // local

} // avro
