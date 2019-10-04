#ifndef URSA_LIB_LIBRARY_FIELD_H_
#define URSA_LIB_LIBRARY_FIELD_H_

#include "mesh/field.h"
#include "mesh/topology.h"

#include <vector>

namespace ursa
{

namespace library
{

class TriangleTextureField : public Field<Simplex<Lagrange>,std::vector<real_t>>
{
public:
  TriangleTextureField( const Topology<Simplex<Lagrange>>& topology );

private:
};



} // library

} // ursa

#endif
