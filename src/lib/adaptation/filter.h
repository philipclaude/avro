//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_ADAPTATION_FILTER_H_
#define avro_LIB_ADAPTATION_FILTER_H_

#include "mesh/points.h"

#include <memory>

namespace avro
{

template<typename type> class Topology;
class Simplex;
template<typename type> class MetricField;

class PointCloud;
class KdTreeNd;

template<typename type> class Insert;

class Filter : public Points
{

public:
  Filter( const coord_t dim );

  void createPermanent( const real_t* x );
  void createCandidate( const index_t n0 , const index_t n1 ,
                        const real_t s , const real_t* x ,
                        Entity * e , const real_t* params );
  void accept( const index_t k , const index_t idx );

  index_t edge( const index_t k , const coord_t p )
  {
    avro_assert( p<=1 );
    avro_assert( !permanent(k) );
    int n;
    if (p==0)
      n = node0_[k];
    else
      n = node1_[k];
    avro_assert(n>=0);
    return index_t(n);
  }

  void generateCandidates( Topology<Simplex>& topology ,
                           MetricField<Simplex>& , real_t lmax , Insert<Simplex>& inserter );

  bool tooclose( Points& v , const index_t k , MetricField<Simplex>& metric ,
                 const real_t lmin , index_t nn );

  bool permanent( const index_t k ) const
    { return node0_[k]<0; }

  index_t candidate( const index_t k ) const
    { return candidates_[k]; }

  index_t nb_candidates() const
    { return candidates_.size(); }

  void clearCandidates();
  void print() const;

  real_t minlength() const { return lmin_; }
  real_t maxlength() const { return lmax_; }
  index_t nb_long() const { return nb_long_; }

private:
  // edges for the insertion
  std::vector<int> node0_;
  std::vector<int> node1_;
  std::vector<real_t> s_;
  std::vector<index_t> candidates_;
  std::vector<int> idx_;

  real_t lmin_,lmax_;
  index_t nb_long_;
};

} // avro

#endif
