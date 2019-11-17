#ifndef LUNA_LIB_MESH_BUILDER_H_
#define LUNA_LIB_MESH_BUILDER_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <vector>

namespace luna
{

template<typename Shape_t,typename Master_t>
class Builder : public Data<index_t>
{

public:
  Builder( const Topology<Shape_t>& topology , const Master_t& master );

  void build();

  void transfer( Topology<Master_t>& F ) const;
  template<typename T> void transfer( Field<Shape_t,Master_t,T>& F ) const;
  template<typename MasterFrom_t,typename T> void transfer( const Field<Shape_t,MasterFrom_t,T>& Fx , Field<Shape_t,Master_t,T>& Fy ) const;

private:
  const Topology<Shape_t>& topology_;
  const Master_t&          master_;
};

}

#endif
