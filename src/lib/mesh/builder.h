#ifndef URSA_LIB_MESH_BUILDER_H_
#define URSA_LIB_MESH_BUILDER_H_

#include "mesh/topology.h"
#include "mesh/vertices.h"

#include <vector>

namespace ursa
{

template<typename Shape_t,typename Master_t>
class Builder : public Data<index_t>
{

public:
  Builder( const Topology<Shape_t>& topology , const Master_t& master );

  void build();

  void transfer( Topology<Master_t>& F ) const;
  //template<typename Basis_t,typename T> void transfer( const Field<Shape_t,Basis_t,T>& Fx , Field<Shape_t,Master_t,T>& Fy ) const;

private:
  const Topology<Shape_t>& topology_;
  const Master_t&          master_;
};

}

#endif
