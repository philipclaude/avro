//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MESH_INVERSE_H_
#define avro_LIB_MESH_INVERSE_H_

#include "avro_types.h"

#include <unordered_set>
#include <vector>

namespace avro
{

template<typename type> class Topology;

template<typename type> class Cavity;

template<typename type>
class InverseTopology
{

public:
  InverseTopology( const Topology<type>& topology );

  void build();

  index_t elem( const index_t k ) const { return elem_[k]; }
  std::unordered_set<index_t> ball( index_t k ) const { return ball_[k]; }

  void ball( index_t k , std::vector<index_t>& B ) const;
  void getball( index_t p , index_t k , std::unordered_set<index_t>& B ) const;

  void shell( index_t p , index_t k , std::vector<index_t>& S ) const;
  void shell( index_t t0, index_t t1, index_t t2 , std::vector<index_t>& S ) const;

  index_t nb() const;
  bool created() const;

  void update( Cavity<type>& cavity , bool delay=false );
  void create( index_t nb_new );
  void remove( index_t k );
  void decrement( index_t bar );
  void add( index_t k );
  index_t find( index_t k ) const;

  void clear() { elem_.clear(); }

  bool check() const;

  void copy( const InverseTopology<type>& inverse );

  void print( bool balls=true ) const;

  void use_ball( bool x ) { use_ball_ = x; }


private:

  void getshell( index_t p , index_t q , index_t k , std::unordered_set<index_t>& B ) const;
  void getshell( index_t t0 , index_t t1 , index_t t2 , index_t k , std::unordered_set<index_t>& B ) const;

  const Topology<type>& topology_;
  std::vector<index_t> elem_;
  std::vector<std::unordered_set<index_t>> ball_;

  bool use_ball_;

};

} // avro

#endif
