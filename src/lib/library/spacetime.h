#ifndef avro_LIB_LIBRARY_SPACETIME_H_
#define avro_LIB_LIBRARY_SPACETIME_H_

#include "mesh/mesh.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include <map>

namespace avro
{

class Entity;

template<typename type>
class Topology_Spacetime : public Topology<type>
{
public:
  Topology_Spacetime( const Topology<type>& topology );

  void extract();

  void write( const std::string& filename );
  void read( const std::string& filename );

  void clear();

private:
  Points points_;
  const Topology<type>& topology_;
  //Topology<type> slice_;
  std::map<Entity*,Topology<type>*> entity2topology_;
};

} // avro

#endif
