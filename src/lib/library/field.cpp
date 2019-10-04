#include "library/field.h"

namespace ursa
{

namespace library
{

TriangleTextureField::TriangleTextureField( const Topology<Simplex<Lagrange>>& topology ) :
  Field<Simplex<Lagrange>,std::vector<real_t>>(topology,1)
{}

} // library

} // ursa
