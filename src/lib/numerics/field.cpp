#include "common/data.h"

#include "geometrics/egads.h"
#include "geometrics/primitive.h"
#include "geometrics/plc.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/topology.h"

#include "numerics/field.h"
#include "numerics/matrix.h"

namespace ursa
{

template<typename Basis,typename T>
Field<Simplex<Basis>,T>::Field( const Topology<Shape_t>& topology , coord_t order ) :
  topology_(topology),
  master_(topology.number(),order)
{
  printf("constructing field with order %u\n",master_.order());
}

template<typename T>
Field<Polytope,T>::Field( Topology<Polytope>& topology , coord_t order ) :
  topology_(topology),
  master_(topology.number(),topology.order(),topology.master().incidence())
{
  printf("constructing field with order %u\n",master_.order());
}

template class Field< Simplex<Lagrange> , real_t >;
template class Field< Simplex<Bezier> , real_t >;
template class Field< Polytope , real_t >;

template class Field< Simplex<Lagrange> , geometrics::Primitive* >;
template class Field< Simplex<Lagrange> , numerics::SymMatrixD<real_t> >;

} // ursa
