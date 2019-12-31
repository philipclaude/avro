#ifndef avro_LIB_MESH_BUILDER_H_
#define avro_LIB_MESH_BUILDER_H_

#include "mesh/topology.h"
#include "mesh/points.h"

#include <vector>

namespace avro
{

template<typename type>
class Builder : public Table<index_t>
{

public:
  Builder( const Topology<type>& topology , coord_t order , BasisFunctionCategory category );

  void build();

  void transfer( Topology<type>& F ) const;
  template<typename T> void transfer( Field<type,T>& F ) const;
  template<typename T> void transfer( const Field<type,T>& Fx , Field<type,T>& Fy ) const;

private:
  const Topology<type>& topology_;
   type master_;
};

}

#endif
