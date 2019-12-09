#ifndef LUNA_LIB_LIBRARY_MESHB_H_
#define LUNA_LIB_LIBRARY_MESHB_H_

#include "mesh/topology.h"

#include <string>

namespace luna
{

namespace meshb
{

template<typename type>
void read_meshb( const std::string& meshfile , Topology<type>& topology );

template<typename type>
void read_solb( const std::string& solfile );

} // meshb

} // luna

#endif
