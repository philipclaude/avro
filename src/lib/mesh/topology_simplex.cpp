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

template<>
void
Topology<Simplex>::get_triangles( std::vector<index_t>& t ) const
{
  std::vector<index_t> tk;
  std::set<std::string> MAP;
  std::vector<index_t> triangle(3);
  std::string s;

  // loop through all the cells
  for (index_t k=0;k<nb();k++)
  {
    if (ghost(k)) continue;

    const index_t* v0 = (*this)(k);

    // get the edges of this cell
    master_.get_triangles( v0 , nv(k) , tk );

    // add the edges
    for (index_t j=0;j<tk.size()/3;j++)
    {
      index_t p0 = tk[3*j];
      index_t p1 = tk[3*j+1];
      index_t p2 = tk[3*j+2];

      if (p0<points_.nb_ghost() || p1<points_.nb_ghost() || p2<points_.nb_ghost())
        continue;

      triangle[0] = p0;
      triangle[1] = p1;
      triangle[2] = p2;
      s = unique_label(triangle);

      if (MAP.find(s)==MAP.end())
      {
        MAP.insert(s);
        t.push_back(p0);
        t.push_back(p1);
        t.push_back(p2);
        //printf("adding triangle %s\n",s.c_str());
      }
      //else
      //printf("triangle %s exists!\n",s.c_str());
    }
  }
}

template class Topology<Simplex>;
template class Tree<Topology<Simplex>>;

} // luna
