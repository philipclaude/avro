#ifndef LUNA_LIB_MASTER_POLYTOPE_H_
#define LUNA_LIB_MASTER_POLYTOPE_H_

#include "common/data.h"

#include "master/master.h"
#include "master/simplex.h"

namespace luna
{

template<typename type> class Topology;

class Polytope : public Master<Polytope>
{

public:
  Polytope( coord_t number , coord_t order );
  Polytope( Topology<Polytope>& topology , const coord_t order );

  Data<int>& incidence() { return *incidence_.get(); }
  const Data<int>& incidence() const { return *incidence_.get(); }

  void get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const;

private:
  Simplex simplex_;
  std::shared_ptr<Data<int>> incidence_;
};

} // luna

#endif
