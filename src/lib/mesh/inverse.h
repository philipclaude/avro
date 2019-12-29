#ifndef luma_LIB_MESH_INVERSE_H_
#define luma_LIB_MESH_INVERSE_H_

#include "common/types.h"

#include <unordered_set>
#include <vector>

namespace luma
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

  void ball( index_t k , std::vector<index_t>& B ) const;
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

private:

  void getball( index_t p , index_t k , std::unordered_set<index_t>& B ) const;
  void getshell( index_t p , index_t q , index_t k , std::unordered_set<index_t>& B ) const;
  void getshell( index_t t0 , index_t t1 , index_t t2 , index_t k , std::unordered_set<index_t>& B ) const;

  const Topology<type>& topology_;
  std::vector<index_t> elem_;

};

} // luma

#endif
