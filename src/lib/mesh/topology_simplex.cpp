//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/cavity.h"

#include "common/tree.hpp"

#include "mesh/builder.h"
#include "mesh/facets.h"
#include "mesh/topology.h"
#include "mesh/topology.hpp"

#include "numerics/geometry.h"

namespace avro
{

template<>
Topology<Simplex>::Topology( Points& vertices , coord_t number , coord_t order ) :
  TopologyBase(vertices,number,TableLayout_Rectangular,"simplex"),
  element_( number , order ),
  neighbours_(*this),
  inverse_(*this)
{
  set_rank( element_.nb_basis() );
}

template<>
Topology<Simplex>::Topology( Points& points , const Topology<Simplex>& linear , coord_t order ) :
 Topology(points,linear.number(),order)
{
  Builder<Simplex> builder(linear,element_.order(),BasisFunctionCategory_Lagrange);
  builder.transfer(*this);
}

template<>
void
Topology<Simplex>::orient( index_t* v , const index_t nv , real_t* q )
{
  std::vector<const real_t*> x(nv);
  std::vector<index_t> p = linspace(nv);
  std::vector<index_t> idx(nv);

  if (q!=NULL) avro_implement;

  // options to orient with a given vertex,
  // e.g. when a 2d topology is oriented in 3d
  if (q!=NULL) x.resize(nv+1);

  // loop through the permutations
  do
  {
    for (coord_t j=0;j<nv;j++)
      x[j] = points_[ v[ p[j] ] ];

    if (q!=NULL) x[nv] = q;

    for (coord_t j=0;j<nv;j++)
      idx[j] = v[p[j]];

    if (element_.volume(points_,idx.data(),nv)>0.) break;
    //if (numerics::simplex_volume(x,points_.dim())>0.)
    //  break;

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
Topology<Simplex>::apply( Cavity<Simplex>& cavity )
{
  // this simply applies the cavity operator in terms of element
  // removals and insertions, including the neighbour relationships
  // it is the responsibility of the caller to do the vertex removals

  // the cache needs to be empty if the topology is closed (null boundary)
  if (this->closed())
    avro_assert( this->neighbours_.cache().nb()==0 );

  index_t elem;
  if (cavity.nb_insert() <= cavity.nb_cavity()) {

    // fewer inserted elements than we are deleting
    for (index_t k = 0; k < cavity.nb_insert(); k++) {

      elem = cavity.cavity(k);
      this->neighbours_.remove( elem , false  ); // no erase
      for (index_t j = 0; j < cavity.nv(k); j++)
        this->operator()( elem , j ) = cavity(k,j);
    }

    // delete the remaining cells
    index_t count = 0;
    for (index_t k = cavity.nb_insert(); k < cavity.nb_cavity(); k++) {

      elem = cavity.cavity(k);
      this->neighbours_.remove( elem-count ); // erase
      Topology<Simplex>::remove( elem-count );
      count++;
    }

    // compute the neighbours of the replaced elements
    for (index_t k = 0; k < cavity.nb_insert(); k++)
    {
      this->neighbours_.addElement( cavity.cavity(k) );
      cavity.inserted()[k] = cavity.cavity(k);
    }
  }
  else {

    // more inserted elements than cavity elements
    for (index_t k = 0; k < cavity.nb_cavity(); k++) {
      elem = cavity.cavity(k);
      this->neighbours_.remove( elem , false  );
      for (index_t j = 0; j < cavity.nv(elem); j++)
        this->operator()( elem , j ) = cavity(k,j);
    }

    // add the remaining cells
    this->neighbours_.enlarge( cavity.nb_insert() -cavity.nb_cavity() );
    elem = this->nb();
    for (index_t k = cavity.nb_cavity(); k < cavity.nb_insert(); k++) {
      Topology<Simplex>::add( cavity(k) , cavity.nv(k) );
    }

    for (index_t k = 0; k < cavity.nb_cavity(); k++) {
      this->neighbours_.addElement( cavity.cavity(k) );
      cavity.inserted()[k] = cavity.cavity(k);
    }
    index_t count = 0;
    for (index_t k = cavity.nb_cavity(); k < cavity.nb_insert(); k++) {
      this->neighbours_.addElement( elem+count );
      cavity.inserted()[k] = elem+count;
      count++;
    }
  }

  // update the inverse topology
  this->inverse_.update( cavity , true ); // delay the removal of vertices

  // ensure the neighbours cache is empty (i.e. null boundary)
  if (this->closed()) {
    if (this->neighbours_.cache().nb() != 0) {
      cavity.print();
    }
    avro_assert( this->neighbours_.cache().nb() == 0 );
  }
}

template<>
void
Topology<Simplex>::get_boundary( Topology<Simplex>& bnd ) const
{
  avro_assert( bnd.number()==number_-1 );

  // compute the facets
  Facets facets(*this);
  facets.compute();

  // retrieve boundary facets
  std::vector<index_t> f(number_);
  for (index_t k=0;k<facets.nb();k++)
  {
    if (facets.boundary(k))
    {
      facets.retrieve(k,f);
      bnd.add( f.data() , f.size() );
    }
  }
}

template<>
void
Topology<Simplex>::reserve( index_t n ) {
  Table<index_t>::reserve( n * element_.nb_basis() );
}

template class Topology<Simplex>;
template class Tree<Topology<Simplex>>;
template void Topology<Simplex>::construct( std::shared_ptr<Topology<Polytope>>& node , Topology<Polytope>& ) const;
template void Topology<Simplex>::construct( std::shared_ptr<Topology<Simplex>>& node , Topology<Simplex>& ) const;

} // avro
