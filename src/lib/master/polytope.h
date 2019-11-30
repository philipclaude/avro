#ifndef LUNA_LIB_MASTER_POLYTOPE_H_
#define LUNA_LIB_MASTER_POLYTOPE_H_

#include "common/table.h"

#include "master/master.h"
#include "master/simplex.h"

namespace luna
{

template<typename type> class Topology;

class Polytope : public Master<Polytope>
{

public:
  Polytope( coord_t number , coord_t order , const Table<int>& incidence );
  Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence );

  const Table<int>& incidence() const { return incidence_; }

  void get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const;

  void facet( const index_t* v , index_t j , std::vector<index_t>& f ) const
    { luna_implement; }

private:
  Simplex simplex_;
  const Table<int>& incidence_;
};

} // luna

#endif
