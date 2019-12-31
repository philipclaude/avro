#ifndef avro_LIB_MASTER_POLYTOPE_H_
#define avro_LIB_MASTER_POLYTOPE_H_

#include "common/table.h"

#include "master/master.h"
#include "master/simplex.h"

namespace avro
{

template<typename type> class Topology;

class Polytope : public Master<Polytope>
{

public:
  Polytope( coord_t number , coord_t order , const Table<int>& incidence );
  Polytope( Topology<Polytope>& topology , const coord_t order , const Table<int>& incidence );

  static std::string type_name() { return "polytope"; }

  const Table<int>& incidence() const { return incidence_; }

  void get_edges( const index_t* v0 , index_t nv , std::vector<index_t>& ek ) const;

  real_t volume( const Points& points , const index_t* v , index_t nv ) const;

  void facet( const index_t* v , index_t j , std::vector<index_t>& f ) const
    { avro_implement; }

private:
  Simplex simplex_;
  const Table<int>& incidence_;
};

} // avro

#endif
