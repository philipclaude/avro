#include "adaptation/adapt.h"
#include "adaptation/metric.h"

#include "mesh/topology.h"

namespace luna
{

template<typename type>
void
edgesInLengthBounds( MetricField<type>& metric_ ,
                     Topology<type>& topology ,
                     std::vector<index_t>& edges ,
                     real_t alpha=-1 )
{
  edges.clear();
  std::vector<index_t> edges0;
  topology.get_edges(edges0);
  std::vector<real_t> lengths( edges0.size()/2 );
  for (index_t k=0;k<edges0.size()/2;k++)
  {
    lengths[k] = metric_.length( edges0[2*k] , edges0[2*k+1] );
    if (alpha<0)
    {
      edges.push_back( edges0[2*k] );
      edges.push_back( edges0[2*k+1] );
    }
    if (lengths[k]<alpha || lengths[k]>1./alpha)
    {
      edges.push_back( edges0[2*k] );
      edges.push_back( edges0[2*k+1] );
    }
  }
}

template<typename type>
void
AdaptThread<type>::swap_edges( real_t qt , index_t npass , bool lcheck )
{
  index_t pass = 0;

  real_t lmin0 = -1;
  real_t lmax0 = -1;
  if (lcheck)
  {
    // compute the min and max lengths coming in
    std::vector<real_t> lens;
    metric_.lengths( topology_ , lens );
    lmin0 = * std::min_element( lens.begin() , lens.end() );
    lmax0 = * std::max_element( lens.begin() , lens.end() );
  }

  printf("-> performing edge swaps with target qt < %g:\n",qt);
  while (true)
  {
    if (pass>npass) break;

    // evaluate the quality on the current topology
    //topology_.evaluate( metric_ );
    //real_t qmin = topology_.worst();
    real_t qmin = worst_quality(topology_,metric_);
    index_t count = 0;
    for (index_t k=0;k<topology_.nb();k++)
    {
      if (topology_.ghost(k)) continue;
      //if (topology_.quality(k)<qt) count++;
      if (metric_.quality(topology_,k)<qt) count++;
    }
    printf("\tpass %lu: qmin = %g, nb_simplices < %g = %lu / %lu\n",
            pass,qmin,qt,count,topology_.nb_real());

    index_t nb_swaps = 0;
    index_t nb_length_rejected = 0;
    index_t nb_edges_considered = 0;

    // now try edge swaps
    nb_swaps = 0;
    nb_length_rejected = 0;
    std::vector<index_t> edges;
    edgesInLengthBounds( metric_ , topology_ , edges , 0.8 );

    printf("\tnb_elem = %lu, nb_edges = %lu, nb_vert = %lu\n",
            topology_.nb(),edges.size()/2,topology_.points().nb());

    edge_swapper_.nb_parameter_rejections() = 0;
    edge_swapper_.nb_parameter_tests() = 0;
    for (coord_t d=0;d<topology_.number();d++)
      edge_swapper_.nb_geometry_rejections(d) = 0;
    edge_swapper_.nb_wake() = 0;
    edge_swapper_.nb_invalid_geometry() = 0;

    std::vector<index_t> elems;
    std::vector<index_t> candidates;
    std::vector<real_t> lens;
    for (index_t k=0;k<edges.size()/2;k++)
    {

      elems.clear();
      candidates.clear();
      edge_swapper_.restart();

      // swapping an edge doesn't change the other edges we might want to swap
      // it will affect the visibility of each swap so they might not be accepted
      index_t e0 = edges[2*k];
      index_t e1 = edges[2*k+1];
      std::vector<index_t> edge = {e0,e1};

      topology_.intersect(edge,elems);
      if (elems.size()==0) continue;
      luna_assert( elems.size()>0 );

      // initial worst quality
      real_t q0 = -1;
      for (index_t j=0;j<elems.size();j++)
      {
        if (topology_.ghost( elems[j] ) )
          continue;
        //if (q0<0 || topology_.quality(elems[j])<q0)
        //  q0 = topology_.quality(elems[j]);
        if (q0<0 || metric_.quality(topology_,elems[j])<q0)
          q0 = metric_.quality(topology_,elems[j]);
      }
      if (q0>qt) continue;

      nb_edges_considered++;

      // list all the vertices in elems
      for (index_t j=0;j<elems.size();j++)
      for (index_t i=0;i<topology_.nv(elems[j]);i++)
        candidates.push_back(topology_(elems[j],i));
      uniquify(candidates);

      // assign the cavity used by the swapper since this is the same
      // for all candidates
      edge_swapper_.setCavity( elems );

      // loop through all candidates
      int m = -1;
      real_t qw = q0;
      for (index_t j=0;j<candidates.size();j++)
      {

        // skip candidates that are endpoints of the initial edge
        if (candidates[j]==e0 || candidates[j]==e1) continue;

        // try the swap
        bool accept = edge_swapper_.apply( candidates[j] , e0 , e1 );
        if (!accept)
        {
          // it was rejected because of geometry topology or visiblity
          continue;
        }

        // option to check the produced lengths
        if (lcheck)
        {
          lens.clear();
          metric_.lengths( edge_swapper_ , lens );
          if (lens.size()==0)
          {
            continue;
          }
          real_t lmin = * std::min_element( lens.begin() , lens.end() );
          real_t lmax = * std::max_element( lens.begin() , lens.end() );
          if (lmin<lmin0 || lmax>lmax0)
          {
            nb_length_rejected++;
            continue;
          }
        }

        //edge_swapper_.evaluate(metric_);
        real_t qw_swap = worst_quality(edge_swapper_,metric_);
        if (qw_swap>qw)
        {
          // check that no elements (ghosts) become duplicated
          if (!edge_swapper_.philipcondition()) continue;
          m  = candidates[j];
          qw = qw_swap;
        }

      }

      if (m<0) continue;

      nb_swaps++;
      edge_swapper_.apply( index_t(m) , e0 , e1 );
      //edge_swapper_.evaluate(metric_);
      topology_.apply(edge_swapper_);
      //topology_.update(edge_swapper_,metric_);

    } // loop over edges

    printf("\t--> performed %lu edge swaps (out of %lu), lrej = %lu\n",nb_swaps,nb_edges_considered,nb_length_rejected);
    printf("\t\tgeometry rejections: ");
    printf("Invalid (%lu)",edge_swapper_.nb_invalid_geometry());
    if (topology_.number()>1) printf(", Edges (%lu)",edge_swapper_.nb_geometry_rejections(1));
    if (topology_.number()>2) printf(", Faces (%lu)",edge_swapper_.nb_geometry_rejections(2));
    if (topology_.number()>3) printf(", Cubes (%lu)",edge_swapper_.nb_geometry_rejections(3));
    printf(", Wake (%lu)",edge_swapper_.nb_wake());
    printf("\n");
    if (nb_swaps==0) break;
    pass++;

  } // loop over passes

}

template class AdaptThread<Simplex>;

} // luna
