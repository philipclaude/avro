#include "adaptation/cavity.h"

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
{
  set_rank( master_.nb_basis() );
}

template<>
Topology<Simplex>::Topology( Points& points , const Topology<Simplex>& linear , coord_t order ) :
 Topology(points,linear.number(),order)
{
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
      }
    }
  }
}

template<>
void
Topology<Simplex>::apply( Cavity<Simplex>& cavity )
{
  // this simply applies the cavity operator in terms of element
  // removals and insertions, including the neighbour relationships
  // it is the responsibility of the caller to do the vertex removals

  // the cache needs to be empty if the topology is closed (null boundary)
  if (this->closed())
    luna_assert( this->neighbours_.cache().nb()==0 );

  index_t elem;

  if (cavity.nb_insert()<=cavity.nb_cavity())
  {
    // fewer inserted elements than we are deleting
    for (index_t k=0;k<cavity.nb_insert();k++)
    {
      elem = cavity.cavity(k);

      this->neighbours_.remove( elem , false  ); // no erase
      for (index_t j=0;j<cavity.nv(k);j++)
        this->operator()( elem , j ) = cavity(k,j);
    }

    // delete the remaining cells
    index_t count = 0;
    for (index_t k=cavity.nb_insert();k<cavity.nb_cavity();k++)
    {
      elem = cavity.cavity(k);
      this->neighbours_.remove( elem-count ); // erase
      Topology<Simplex>::remove( elem-count );
      count++;
    }

    // compute the neighbours of the replaced elements
    for (index_t k=0;k<cavity.nb_insert();k++)
    {
      this->neighbours_.addElement( cavity.cavity(k) );
      cavity.inserted()[k] = cavity.cavity(k);
    }
  }
  else
  {

    // more inserted elements than cavity elements
    for (index_t k=0;k<cavity.nb_cavity();k++)
    {
      elem = cavity.cavity(k);

      this->neighbours_.remove( elem , false  );
      for (index_t j=0;j<cavity.nv(elem);j++)
        this->operator()( elem , j ) = cavity(k,j);
    }

    // add the remaining cells
    this->neighbours_.enlarge( cavity.nb_insert() -cavity.nb_cavity() );
    elem = this->nb();
    for (index_t k=cavity.nb_cavity();k<cavity.nb_insert();k++)
    {
      Topology<Simplex>::add( cavity(k) , cavity.nv(k) );
    }

    for (index_t k=0;k<cavity.nb_cavity();k++)
    {
      this->neighbours_.addElement( cavity.cavity(k) );
      cavity.inserted()[k] = cavity.cavity(k);
    }
    index_t count = 0;
    for (index_t k=cavity.nb_cavity();k<cavity.nb_insert();k++)
    {
      this->neighbours_.addElement( elem+count );
      cavity.inserted()[k] = elem+count;
      count++;
    }
  }

  // update the inverse topology
  this->inverse_.update( cavity , true ); // delay the removal of vertices

  // ensure the neighbours cache is empty (i.e. null boundary)
  if (this->closed())
  {
    if (this->neighbours_.cache().nb()!=0)
    {
      //this->neighbours_.cache().print();
      cavity.print();
    }
    luna_assert( this->neighbours_.cache().nb()==0 );

    // debug check (slow) REMOVE ME!
    //for (index_t k=0;k<this->neighbours_.nb();k++)
    //for (index_t j=0;j<this->number_+1;j++)
    //{
    //  if (this->neighbours_(k,j)<0)
    //    throw "oops";
    //}
  }
}

template class Topology<Simplex>;
template class Tree<Topology<Simplex>>;

} // luna
