#include "common/tree.hpp"

#include "mesh/builder.h"
#include "mesh/topology.h"
#include "mesh/topology.hpp"

#include "numerics/geometry.h"

namespace luna
{

template<>
Topology<Simplex>::Topology( Points& vertices , coord_t number , coord_t order ) :
  TopologyBase(vertices,number,TableLayout_Rectangular,"simplex"),
  master_( number , order ),
  neighbours_(*this),
  inverse_(*this)
{}

template<>
Topology<Simplex>::Topology( Points& points , const Topology<Simplex>& linear , coord_t order ) :
 Topology(points,linear.number(),order)
{
  printf("converting to order %u...\n",master_.order());
  Builder<Simplex> builder(linear,master_.order(),BasisFunctionCategory_Lagrange);
  builder.transfer(*this);
}

template<>
void
Topology<Simplex>::orient( index_t* v , const index_t nv , real_t* q )
{
  std::vector<const real_t*> x(nv);
  std::vector<index_t> p = linspace(nv);

  // options to orient with a given vertex,
  // e.g. when a 2d topology is oriented in 3d
  if (q!=NULL) x.resize(nv+1);

  // loop through the permutations
  do
  {
    for (coord_t j=0;j<nv;j++)
      x[j] = points_[ v[ p[j] ] ];

    if (q!=NULL) x[nv] = q;

    if (numerics::simplex_volume(x,points_.dim())>0.)
      break;

  } while (std::next_permutation(p.begin(),p.end()));

  // save the order
  std::vector<index_t> u(nv);
  for (index_t j=0;j<nv;j++)
    u[j] = v[p[j]];
  for (index_t j=0;j<nv;j++)
    v[j] = u[j];
}

template class Topology<Simplex>;
template class Tree<Topology<Simplex>>;

} // luna
