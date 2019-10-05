#ifndef URSA_LIB_MESH_BUILDER_H_
#define URSA_LIB_MESH_BUILDER_H_

#include "mesh/topology.h"
#include "mesh/vertices.h"

#include <vector>

namespace ursa
{

template<typename MasterFrom_t,typename MasterTo_t,typename dof_t>
class Builder : public Data<index_t>
{

public:
  Builder( const Topology<MasterFrom_t>& topology , const MasterTo_t& masterTo );
  Builder( const Field<MasterFrom_t,dof_t>& field , const MasterTo_t& masterTo );

  void build();

private:
  std::vector<dof_t> dofFrom_;
  std::vector<dof_t> dofTo_;

  const Data<index_t>& topology_;
  const MasterFrom_t&  masterFrom_;
  const MasterTo_t&    masterTo_;
};

}

#endif
