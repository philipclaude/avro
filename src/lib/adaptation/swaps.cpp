//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/adapt.h"
#include "adaptation/geometry.h"
#include "adaptation/metric.h"

#include "mesh/topology.h"

#include <time.h>

namespace avro
{

template<typename type>
void
lengths_in_bounds( MetricField<type>& metric_ ,
                   Topology<type>& topology ,
                   std::vector<index_t>& edges ,
                   real_t alpha=-1 ) {
  edges.clear();
  std::vector<index_t> edges0;
  topology.get_edges(edges0);
  std::vector<real_t> lengths( edges0.size()/2 );
  for (index_t k = 0; k < edges0.size()/2; k++) {
    lengths[k] = metric_.length( edges0[2*k] , edges0[2*k+1] );
    if (alpha < 0) {
      edges.push_back( edges0[2*k] );
      edges.push_back( edges0[2*k+1] );
    }
    if (lengths[k]<alpha || lengths[k]>1./alpha) {
      edges.push_back( edges0[2*k] );
      edges.push_back( edges0[2*k+1] );
    }
  }
}

template<typename type>
bool
AdaptThread<type>::swap_edge( index_t p , index_t q , real_t Q0 , real_t lmin0 , real_t lmax0 ) {

  if (topology_.points().fixed(p) && topology_.points().fixed(q)) return false;

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
  edge_swapper_.set_cavity( shell );

  // loop through all candidates
  int m = -1;
  real_t qw = Q0;
  for (index_t j = 0; j < candidates.size(); j++) {

    if (candidates[j]==p || candidates[j]==q) continue;

    // check if the swap is valid in terms of geometry topology
    bool accept = edge_swapper_.valid( candidates[j] , p , q );
    if (!accept) {
      continue;
    }

    // check if the swap is valid in terms of visibility
    accept = edge_swapper_.apply( candidates[j] , p , q );
    if (!accept) {
      continue;
    }

    // option to check the produced lengths
    #if 0
    std::vector<real_t> lens;
    metric_.lengths( edge_swapper_ , lens );
    if (lens.size()==0) return false;
    real_t lmin = * std::min_element( lens.begin() , lens.end() );
    real_t lmax = * std::max_element( lens.begin() , lens.end() );
    if (lmin0>0 && lmin<lmin0) continue;
    if (lmax0>0 && lmax>lmax0) continue;
    #endif

    edge_swapper_.cavity_quality().resize( edge_swapper_.nb() );
    real_t qws = worst_quality(edge_swapper_,metric_ , edge_swapper_.cavity_quality().data() );
    if (qws > qw) {
      if (!edge_swapper_.has_unique_elems()) continue;
      m  = candidates[j];
      qw = qws;
    }
  }

  if (m < 0) return false;

  edge_swapper_.apply( index_t(m) , p , q );
  topology_.apply(edge_swapper_);

  return true;
}

template<typename type>
void
AdaptThread<type>::swap_edges( real_t qt , index_t npass , bool lcheck )
{
  index_t pass = 0;

  printf("-> performing edge swaps with target qt < %g:\n",qt);
  while (true) {

    if (pass > npass) break;

    clock_t TIME0,TIME1;
    clock_t PASS_T0 = clock();

    real_t candidate_time = 0.0;
    real_t metric_time = 0.0;
    real_t shell_time = 0.0;
    real_t apply_time = 0.0;
    real_t candidate_apply_time = 0.0;
    real_t candidate_metric_time = 0.0;

    // evaluate the quality on the current topology
    TIME0 = clock();
    real_t qmin = worst_quality(topology_,metric_);
    index_t count = 0;
    for (index_t k=0;k<topology_.nb();k++) {
      if (topology_.ghost(k)) continue;
      if (metric_.quality(topology_,k)<qt) count++;
    }
    TIME1 = clock();
    metric_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);
    printf("\tpass %lu: qmin = %g, nb_simplices < %g = %lu / %lu\n",
            pass,qmin,qt,count,topology_.nb_real());

    index_t nb_swaps = 0;
    index_t nb_edges_considered = 0;

    // now try edge swaps
    nb_swaps = 0;
    std::vector<index_t> edges;
    TIME0 = clock();
    lengths_in_bounds( metric_ , topology_ , edges , 0.8 );
    TIME1 = clock();
    metric_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);

    printf("\tnb_elem = %lu, nb_edges = %lu, nb_vert = %lu\n",
            topology_.nb(),edges.size()/2,topology_.points().nb());

    edge_swapper_.nb_parameter_rejections() = 0;
    edge_swapper_.nb_parameter_tests() = 0;
    for (coord_t d=0;d<topology_.number();d++)
      edge_swapper_.nb_geometry_rejections(d) = 0;
    edge_swapper_.nb_interior() = 0;
    edge_swapper_.nb_invalid_geometry() = 0;

    std::vector<index_t> elems;
    std::vector<index_t> candidates;
    std::unordered_set<index_t> candidate_set;
    std::vector<real_t> lens;
    for (index_t k = 0; k < edges.size()/2; k++) {

      elems.clear();
      candidates.clear();
      candidate_set.clear();
      edge_swapper_.restart();

      // swapping an edge doesn't change the other edges we might want to swap
      // it will affect the visibility of each swap so they might not be accepted
      index_t e0 = edges[2*k];
      index_t e1 = edges[2*k+1];

      // skip fixed edges (when working in parallel)
      if (topology_.points().fixed(e0) && topology_.points().fixed(e1)) {
        topology_.points().age( e0 )++;
        topology_.points().age( e1 )++;
        continue;
      }

      // we can skip a shell calculation since this swap is never valid
      Entity* g = edge_swapper_.geometry(e0,e1);
      if (g != nullptr && g->number() == 1) continue;

      TIME0 = clock();
      std::vector<index_t> edge = {e0,e1};
      topology_.intersect(edge,elems);
      TIME1 = clock();
      shell_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);
      if (elems.size()==0) continue;

      // initial worst quality
      TIME0 = clock();
      real_t q0 = -1;
      index_t worst_elem = 0;
      std::unordered_set<index_t> bad_elems;
      for (index_t j=0;j<elems.size();j++) {
        if (topology_.ghost( elems[j] ) )
          continue;
        real_t q = metric_.quality(topology_,elems[j]);
        if (q0 < 0 || q  < q0) {
          q0 = q;
          worst_elem = elems[j];
        }
        if (q < qt) bad_elems.insert(elems[j]);
      }
      if (q0 > qt) continue;
      TIME1 = clock();
      metric_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);

      nb_edges_considered++;

      TIME0 = clock();

      // list all the swap candidates
      for (index_t i = 0; i < topology_.nv(worst_elem); i++) {
        candidates.push_back( topology_(worst_elem,i) );
        candidate_set.insert( topology_(worst_elem,i) );
      }

      #if 0
      for (index_t j = 0; j < elems.size(); j++) {
        if (bad_elems.find(elems[j]) == bad_elems.end()) continue; // only consider bad elements
        for (index_t i = 0; i < topology_.nv(elems[j]); i++) {
          index_t m = topology_(elems[j],i);
          if (candidate_set.find(m) != candidate_set.end()) continue;
          candidate_set.insert(m);
          candidates.push_back(m);
        }
      }
      #endif

      // assign the cavity used by the swapper since this is the same
      // for all candidates
      edge_swapper_.set_cavity( elems );

      // the boundary of the swapper is the same for all candidates, so compute it once
      // and let the swapper know that the boundary is saved
      edge_swapper_.save_boundary() = false;
      edge_swapper_.boundary().clear();
      edge_swapper_.clear();
      for (index_t j = 0; j < elems.size(); j++)
        edge_swapper_.add_cavity(elems[j]);
      edge_swapper_.compute_boundary();
      edge_swapper_.save_boundary() = true;

      // loop through all candidates
      int m = -1;
      real_t qw = q0;
      for (index_t j = 0; j < candidates.size(); j++) {

        // skip candidates that are endpoints of the initial edge
        if (candidates[j] == e0 || candidates[j] == e1) continue;

        // try the swap
        clock_t A0 = clock();
        bool accept = edge_swapper_.apply( candidates[j] , e0 , e1 );
        clock_t A1 = clock();
        candidate_apply_time += real_t(A1 - A0)/real_t(CLOCKS_PER_SEC);
        if (!accept) {
          // it was rejected because of geometry topology or visiblity
          continue;
        }

        A0 = clock();
        edge_swapper_.cavity_quality().resize( edge_swapper_.nb() );
        real_t qw_swap = worst_quality(edge_swapper_,metric_,edge_swapper_.cavity_quality().data() );
        if (qw_swap > qw) {
          // check that no elements (ghosts) become duplicated
          // this might only be needed for curved geometries
          if (!edge_swapper_.has_unique_elems()) continue;
          if (!edge_swapper_.closed_boundary()) continue;
          m  = candidates[j];
          qw = qw_swap;
        }
        A1 = clock();
        candidate_metric_time += real_t(A1-A0)/real_t(CLOCKS_PER_SEC);
        if (m >= 0) break; // if we found a swap candidate, we're done (save extra computation)
      }
      TIME1 = clock();
      candidate_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);

      // check if there was no candidate
      if (m < 0) continue;

      // an acceptable swap was found, reset the age of the vertices involved
      topology_.points().age(e0) = 0;
      topology_.points().age(e1) = 0;

      TIME0 = clock();
      nb_swaps++;
      //edge_swapper_.apply( index_t(m) , e0 , e1 ); // not needed if we immediately break when a swap is found
      topology_.apply(edge_swapper_);
      TIME1 = clock();
      apply_time += real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC);

    } // loop over edges

    printf("\t--> performed %lu edge swaps (out of %lu)\n",nb_swaps,nb_edges_considered);
    printf("\t\tgeometry rejections: ");
    printf("Invalid (%lu)",edge_swapper_.nb_invalid_geometry());
    if (topology_.number()>1) printf(", Edges (%lu)",edge_swapper_.nb_geometry_rejections(1));
    if (topology_.number()>2) printf(", Faces (%lu)",edge_swapper_.nb_geometry_rejections(2));
    if (topology_.number()>3) printf(", Cubes (%lu)",edge_swapper_.nb_geometry_rejections(3));
    printf(", Interior (%lu)",edge_swapper_.nb_interior());
    printf("\n");
    if (nb_swaps==0) break;
    pass++;

    printf("\t\tcandidate time = %g\n",candidate_time);
    printf("\t\tmetric time = %g\n",metric_time);
    printf("\t\tshell time = %g\n",shell_time);
    printf("\t\tapply time = %g\n",apply_time);
    printf("\t\tcandidate apply time = %g\n",candidate_apply_time);
    printf("\t\tcandidate metric time = %g\n",candidate_metric_time);

    printf("\t\tpass time = %g\n",real_t(clock()-PASS_T0)/real_t(CLOCKS_PER_SEC));

  } // loop over passes

}

template<typename type>
void
AdaptThread<type>::swap_cells( real_t qt , index_t npass )
{

  printf("-> performing edge swaps with target qt < %g:\n",qt);

  index_t nb_swaps_total = 0;
  index_t pass = 0;

  std::unordered_set<index_t> candidate_set;
  std::vector<index_t> candidates;

  index_t nb_candidates_avg = 0;
  index_t counter = 0;

  std::vector<index_t> shell;
  index_t nb_swaps0 = 0;
  while (true) {

    index_t nb_swaps = 0;

    std::unordered_set<std::string> edges_considered;

    for (index_t k = 0; k < topology_.nb(); k++) {

      if (topology_.ghost(k)) continue;

      // evaluate the quality of this element
      real_t q0 = metric_.quality( topology_ , k );
      if (q0 > qt) continue;

      // loop through the edges of this element
      std::vector<index_t> edges;
      topology_.element().get_edges( topology_(k) , topology_.nv(k) , edges );

      for (index_t j = 0; j < edges.size()/2; j++) {

        index_t e0 = edges[2*j];
        index_t e1 = edges[2*j+1];

        if (e0 < topology_.points().nb_ghost()) continue;
        if (e1 < topology_.points().nb_ghost()) continue;

        if (e0 > e1) std::swap(e0,e1);

        // only loop through edges that have yet to be considered
        std::string label = std::to_string(e0) + "-" + std::to_string(e1);
        if (edges_considered.find(label) != edges_considered.end()) continue;
        edges_considered.insert(label);

        //real_t lj = metric_.length( topology_.points() , e0 , e1 );
        //if (lj < 1.2 && lj > 0.8) continue;

        // get the shell
        shell.clear();
        topology_.intersect({e0,e1},shell);
        if (shell.size() == 0) continue; // this probably isn't possible, but just being safe

        // list all the swap candidates
        candidate_set.clear();
        candidates.clear();
        for (index_t i = 0; i < shell.size(); i++) {
          //if (shell[i] != k) continue;
          if (topology_.ghost(shell[i])) continue;
          for (index_t p = 0; p < topology_.nv(shell[i]); p++) {
            index_t m = topology_(shell[i],p);
            if (m < topology_.points().nb_ghost()) continue;
            if (candidate_set.find(m) != candidate_set.end()) continue;
            candidate_set.insert(m);
            candidates.push_back(m);
          }
        }
        nb_candidates_avg += candidates.size();
        counter++;

        // assign the cavity used by the swapper since this is the same
        // for all candidates
        edge_swapper_.set_cavity( shell );

        // the boundary of the swapper is the same for all candidates, so compute it once
        // and let the swapper know that the boundary is saved
        edge_swapper_.save_boundary() = false;
        edge_swapper_.boundary().clear();
        edge_swapper_.clear();
        for (index_t i = 0; i < shell.size(); i++)
          edge_swapper_.add_cavity(shell[i]);
        edge_swapper_.compute_boundary();
        edge_swapper_.save_boundary() = true;

        // loop through all candidates
        int m = -1;
        real_t qw = q0;
        for (index_t i = 0; i < candidates.size(); i++) {

          // skip candidates that are endpoints of the initial edge
          if (candidates[i] == e0 || candidates[i] == e1) continue;

          // try the swap
          bool accept = edge_swapper_.apply( candidates[i] , e0 , e1 );
          if (!accept) {
            // it was rejected because of geometry topology or visiblity
            continue;
          }

          edge_swapper_.cavity_quality().resize( edge_swapper_.nb() );
          real_t qw_swap = worst_quality(edge_swapper_,metric_,edge_swapper_.cavity_quality().data() );
          if (qw_swap > qw) {
            // check that no elements (ghosts) become duplicated
            // this might only be needed for curved geometries
            //if (!edge_swapper_.has_unique_elems()) continue;
            //if (!edge_swapper_.closed_boundary()) continue;
            m  = candidates[i];
            qw = qw_swap;
          }
          if (m >= 0) break; // if we found a swap candidate, we're done (save some computation)
        } // loop over candidates

        // check if there was no candidate
        if (m < 0) continue;

        // an acceptable swap was found, reset the age of the vertices involved
        topology_.points().age(e0) = 0;
        topology_.points().age(e1) = 0;

        nb_swaps++;
        topology_.apply(edge_swapper_);
        break;
      } // loop over edges
    }

    if (pass == 0) nb_swaps0 = nb_swaps;

    printf("--> pass %lu, performed %lu swaps\n",pass,nb_swaps);
    if (nb_swaps == 0) break;
    if (100*nb_swaps < nb_swaps0) break;
    nb_swaps_total += nb_swaps;
    pass++;
    if (pass > npass) break;
  }

  printf("performed %lu swaps in total\n",nb_swaps_total);
  printf("average number of candidates = %g\n",real_t(nb_candidates_avg)/counter);

}

template<typename type>
EdgeSwap<type>::EdgeSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("edge swapper");
  nb_geometry_rejections_.resize( _topology.number() );
}

template<typename type>
bool
EdgeSwap<type>::visible_geometry( index_t p , index_t e0 , index_t e1 , Entity* face )
{
  // check if the edge swap is valid in the parameter space
  // 2d swaps area already checks, 4d is limited to linear geometries for now
  //if (this->topology_.number()!=3) return true;
  if (!this->curved_) return true;
  if (face->interior()) return true;

  this->geometry_cavity_.set_entity(face);

  // based on the type of face, we may need to flip the sign for the volume calculation
  this->geometry_cavity_.sign() = face->sign();

  // extract the geometry cavity
  this->extract_geometry( face, {e0,e1} );
  avro_assert( this->geometry_topology_.nb()>0 );

  if (!this->geometry_topology_.closed())
    avro_assert_msg( this->geometry_topology_.nb()==2 ,
                    "|G| = %lu, |swap ghosts| = %lu" , this->geometry_topology_.nb() , this->nb_ghost() );

  if (this->geometry_topology_.closed())
    avro_assert( this->geometry_topology_.nb_real()==2 );

  if (this->topology_.element().parameter())
  {
    this->convert_to_parameter(face);
    this->geometry_cavity_.sign() = 1;
  }

  // ensure the e0 is visible to the cavity boundary
  bool accept = this->geometry_cavity_.compute( this->v2u_[e0] , this->u_[this->v2u_[e0]] , this->S_ );
  avro_assert( accept );

  if (this->topology_.element().parameter())
  {
    this->convert_to_physical();
  }

  // reset the geometry checker to this face (which also computes the normals at the vertices of the geometry topology)
  this->geometry_inspector_.reset(face);
  int s = this->geometry_inspector_.signof( this->geometry_cavity_ );
  avro_assert_msg( s > 0 , "negative orientation for edge (%lu,%lu) with vertex %lu" , e0,e1,e0 );
  avro_assert( this->geometry_inspector_.positive_volumes(this->geometry_cavity_,this->geometry_cavity_.sign()));

  if (this->topology_.element().parameter())
  {
    this->convert_to_parameter(face);
  }

  // ensure e1 is visible to the cavity boundary
  accept = this->geometry_cavity_.compute( this->v2u_[e1] , this->u_[this->v2u_[e1]] , this->S_ );
  avro_assert( accept );

  if (this->topology_.element().parameter())
  {
    this->convert_to_physical();
  }

  s = this->geometry_inspector_.signof( this->geometry_cavity_ );
  avro_assert_msg( s > 0 , "negative orientation for edge (%lu,%lu) with vertex %lu" , e0,e1,e1 );
  avro_assert( this->geometry_inspector_.positive_volumes(this->geometry_cavity_,this->geometry_cavity_.sign()));

  if (this->topology_.element().parameter())
  {
    this->convert_to_parameter(face);
  }

  // check for visibility of p
  nb_parameter_tests_++;
  accept = this->geometry_cavity_.compute( this->v2u_[p] , this->u_[this->v2u_[p]] , this->S_ );
  avro_assert( this->geometry_cavity_.nb_real()==2 );
  if (!accept)
  {
    nb_parameter_rejections_++;
    return false;
  }

  if (this->topology_.element().parameter())
  {
    this->convert_to_physical();
  }

  s = this->geometry_inspector_.signof( this->geometry_cavity_ );
  if (s<0)
  {
    return false;
  }

  if (this->geometry_inspector_.invalidates_geometry(this->geometry_cavity_))
  {
    nb_invalid_geometry_++;
    return false;
  }

  // the geometry cavity should have positive volumes
  avro_assert( this->geometry_inspector_.positive_volumes(this->geometry_cavity_,this->geometry_cavity_.sign()));

  return true;
}

template<typename type>
bool
EdgeSwap<type>::valid( const index_t p , const index_t e0 , const index_t e1 )
{
  // topology checks
  if (this->topology_.points().fixed(e0) &&
      this->topology_.points().fixed(e1))
    return false;
  if (p<this->topology_.points().nb_ghost()) return false;
  if (e0<this->topology_.points().nb_ghost()) return false;
  if (e1<this->topology_.points().nb_ghost()) return false;

  std::vector<index_t> edge = {e0,e1};
  if (this->C_.empty()) {
    this->topology_.intersect(edge,this->C_);
  }
  if (this->C_.size()==0) return false;

  // check if the two points lie on an Edge entity
  Entity* ge = this->geometry(e0,e1);

  // check if this edge is in the volume
  if (ge==NULL)
  {
    // edge is in the volume, always valid
    return true;
  }

  // check if this edge is on a geometry Edge (can't swap these)
  if (ge!=NULL && ge->number()<2)
  {
    nb_geometry_rejections_[ge->number()]++;
    return false;
  }

  // check the entity on the re-inserted vertex
  if (ge!=NULL)
  {
    Entity* ep = this->topology_.points().entity(p);
    if (ep==NULL) return false;
    if (ep!=ge)
    {
      if (ge->above(ep))
      {} // the swap is okay
      else
      {
        nb_geometry_rejections_[ge->number()]++;
        return false;
      }
    }
  }

  if (ge->interior())
  {
    nb_interior_++;
  }

  return true;
}

template<typename type>
bool
EdgeSwap<type>::apply( const index_t p , const index_t e0 , const index_t e1 ) {

  // check if the swap is valid in terms of geometry
  if (!valid(p,e0,e1)) return false; // computes cavity C_

  this->info_ = "trying to swap edge (" + stringify<index_t>(e0) + "/"
                + stringify<index_t>(e1) +")"+
                " with reinsertion " + stringify<index_t>(p);

  // compute the original number of ghost elements
  index_t nb_ghost0 = 0;
  for (index_t j = 0; j < this->C_.size(); j++)
    if (this->topology_.ghost(this->C_[j]))
      nb_ghost0++;

  Entity* ge = this->geometry(e0,e1);
  this->set_entity(ge);

  if (this->topology_.element().parameter()) {
    avro_assert( ge!=nullptr );
    avro_assert( ge->number()==2 ); // otherwise this should have been caught by 'valid'
    avro_assert( this->C_.size()==2 );

    bool accept;

    this->Cavity<type>::clear();
    for (index_t j = 0; j < this->C_.size(); j++) {
      this->add_cavity( this->C_[j] );
    }

    if (!visible_geometry(p,e0,e1,ge)) {
      return false;
    }

    this->check_visibility(false);
    accept = this->compute( p ,this->topology_.points()[p],this->C_);
    avro_assert(accept);
    this->check_visibility(true);

    return true;
  }

  // apply the operator, checking visibility
  this->enlarge_ = false;
  bool accept = this->compute( p , this->topology_.points()[p] , this->C_ );
  if (!accept) return false;

  // check if all produce elements have a positive determinant of implied metric
  if (!this->positive_implied_metrics())
    return false;

  // check visibility in the parametric space
  if (ge != NULL && ge->number() == 2 && this->topology_.number() > 2) {

    if (!visible_geometry(p,e0,e1,ge)) {
      //printf("swap not visible in parameter space!\n");
      return false;
    }
  }

  // count the resulting number of ghost elements
  index_t nb_ghost = 0;
  for (index_t j = 0; j < this->nb(); j++)
    if (this->ghost(j))
      nb_ghost++;

  // check if only ghosts are created
  index_t only_ghost = true;
  for (index_t k = 0; k < this->nb(); k++) {
    if (!this->ghost(k)) {
      only_ghost = false;
      break;
    }
  }
  if (only_ghost) return false;

  // do not allow swaps which change the number of ghosts for 3-simplices
  if (this->topology_.number() == 3 && nb_ghost0 != nb_ghost)
    return false;

  if (!accept) return false;

  return accept;
}

template<typename type>
FacetSwap<type>::FacetSwap( Topology<type>& _topology ) :
  Primitive<type>(_topology)
{
  this->setName("facet swapper");
}

template<typename type>
bool
FacetSwap<type>::valid( const index_t p , const index_t k0 , const index_t k1 )
{
  if (p<=this->topology_.points().nb_ghost()) return false;
  if (this->topology_.ghost(k0) || this->topology_.ghost(k1))
    return false;
  return true;
}

template<typename type>
bool
FacetSwap<type>::apply( const index_t p , const index_t k1 , const index_t k2 )
{
  // check if the facet swap is valid
  if (!valid(p,k1,k2)) return false;
  this->enlarge_ = false;

  // check if the swap is visible
  std::vector<index_t> elems = {k1,k2};
  bool accept = this->compute( p , this->topology_.points()[p] , elems );

  // check if all produced elements have a positive determinant of implied metric
  if (!this->positive_implied_metrics())
    return false;
  return accept;
}

template class FacetSwap<Simplex>;
template class EdgeSwap<Simplex>;
template class AdaptThread<Simplex>;

} // avro
