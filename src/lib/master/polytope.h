#ifndef LUNA_LIB_MASTER_POLYTOPE_H_
#define LUNA_LIB_MASTER_POLYTOPE_H_

#include "common/array.h"

#include "master/master.h"
#include "master/simplex.h"

namespace luna
{

template<typename type> class Topology;

class Polytope : public Master<Polytope>
{

public:
  Polytope( coord_t number , coord_t order , const Array<int>& incidence );
  Polytope( Topology<Polytope>& topology , const coord_t order , const Array<int>& incidence );

  const Array<int>& incidence() const { return incidence_; }

  void get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const;

private:
  Simplex simplex_;
  const Array<int>& incidence_;
};

} // luna

#endif
