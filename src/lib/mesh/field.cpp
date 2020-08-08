//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/metric.h"

#include "geometry/entity.h"

#include "element/polytope.h"
#include "element/simplex.h"

#include "mesh/builder.h"
#include "mesh/builder.hpp"
#include "mesh/field.h"
#include "mesh/field.hpp"
#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/matrix.h"

namespace avro
{

void
Fields::get_names( std::vector<std::string>& names , std::vector<std::string>& ids ) const
{
  // get a json describing all the field names with ranks
  std::map<std::string,std::shared_ptr<FieldHolder>>::const_iterator it;
  for (it=fields_.begin();it!=fields_.end();it++)
  {
    index_t rank = it->second->nb_rank();
    for (index_t j=0;j<rank;j++)
    {
      names.push_back( it->first + "-" + it->second->get_name(j) );
      ids.push_back( identifiers_.at(it->first) + "-" + std::to_string(j) );
    }
  }
}

void
Fields::print() const
{
  printf("nb_fields = %lu\n",fields_.size());
  std::map<std::string,std::shared_ptr<FieldHolder>>::const_iterator it;
  for (it=fields_.begin();it!=fields_.end();it++)
  {
    index_t rank = it->second->nb_rank();
    for (index_t j=0;j<rank;j++)
    {
      printf("field[%s] at rank %lu (%s) with id %s\n",it->first.c_str(),j,it->second->get_name(j).c_str(),identifiers_.at(it->first).c_str() );
    }
  }
}

std::string
Fields::id2name( const std::string& id ) const
{
  std::map<std::string,std::string>::const_iterator it;
  for (it=identifiers_.begin();it!=identifiers_.end();it++)
  {
    if (it->second==id)
      return it->first;
  }
  printf("could not find field name associated with id %s\n",id.c_str());
  avro_assert_not_reached;
  return "no-name";
}

template<typename T>
FieldBase<T>::FieldBase( FieldType type , ElementBase& element , TableLayoutCategory category ) :
  Table<index_t>(category),
  data_(1),
  type_(type),
  element_(element)
{}

template<typename T>
Field<Simplex,T>::Field( const Topology<Simplex>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type,element_,TableLayout_Rectangular),
  topology_(topology),
  element_(topology.number(),order)
{
  Table<index_t>::set_rank( element_.nb_basis() );
}

template<typename T>
Field<Polytope,T>::Field( Topology<Polytope>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type,element_),
  topology_(topology),
  element_(topology,order,topology.points().incidence())
{}

template<typename T>
void
Field<Simplex,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    if (topology_.element().order() == element_.order() )
    {
      // a simpler way to build the field
      // copy the connectivity over from the topology
      for (index_t k=0;k<topology_.nb();k++)
      {
        const index_t* v = topology_(k);
        const index_t nv = topology_.nv(k);
        Table<index_t>::add( v , nv );
      }

      T x0(0);
      for (index_t k=0;k<topology_.points().nb();k++)
      {
        this->data_.add( &x0 , 1 );
      }
    }
    else
    {
      // a more complicated way to build the field
      // get the number of unique entries in the field
      // using the number of elements and nb_poly, for now assume p = 1
      Builder<Simplex> builder(topology_,element_.order(),BasisFunctionCategory_None);
      builder.template transfer<T>(*this);
    }
  }
  else if (this->type()==DISCONTINUOUS)
  {
    index_t n = 0;
    std::vector<index_t> dof( element_.nb_basis() );
    for (index_t k=0;k<topology_.nb();k++)
    {
      for (index_t j=0;j<element_.nb_basis();j++)
        dof[j] = n++;
      Table<index_t>::add( dof.data() , dof.size() );
      const std::vector<index_t>& idx = Table<index_t>::data();
      index_t nb_dof = * std::max_element( idx.begin() , idx.end() ) +1;
      this->allocate( nb_dof );
    }
  }
  else
    avro_assert_not_reached;
}

template<typename T>
void
Field<Polytope,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    avro_assert( element_.order()==1 );

    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_(k);
      const index_t nv = topology_.nv(k);
      Table<index_t>::add( v , nv );
    }

    T x0(0);
    for (index_t k=0;k<topology_.points().nb();k++)
    {
      this->data_.add( &x0 , 1 );
    }
  }
  else if (this->type()==DISCONTINUOUS)
  {
    avro_assert( element_.order()==0 );
    T x0(0);
    for (index_t k=0;k<topology_.nb();k++)
    {
      this->data_.add( &x0 , 1 );
    }
  }
  else
    avro_assert_not_reached;
}

template<>
real_t
__at_rank__( const real_t& x , index_t r )
{
  return x;
}

template<>
real_t
__at_rank__( const numerics::SymMatrixD<real_t>& x , index_t r )
{
  avro_implement;
  return x(0,0);
}

template<>
real_t
__at_rank__( const Metric& x , index_t r )
{
  return __at_rank__<numerics::SymMatrixD<real_t>>(x,r);
}

template<>
real_t
__at_rank__( const std::vector<real_t>& x , index_t r )
{
  return x[r];
}

template<>
real_t
__at_rank__( const std::vector<index_t>& x , index_t r )
{
  return x[r];
}

template<>
real_t
__at_rank__( Entity* const& entity , index_t r )
{
  return real_t(entity->identifier());
}

template<>
real_t
__at_rank__( const SurrealS<1,real_t>& x , index_t r )
{
  return x.value();
}

template<>
real_t
__at_rank__( const SurrealD& x , index_t r )
{
  return x.value();
}

template class FieldBase<Metric>;

template class Field< Simplex , real_t >;
template class Field< Polytope , real_t >;

template class Field< Simplex , std::vector<real_t> >;
template class Field< Simplex , std::vector<index_t> >;

template class Field< Simplex , Entity* >;
template class Field< Simplex , numerics::SymMatrixD<real_t> >;

template class Field< Simplex , Metric >;

template class Field< Simplex , SurrealS<1,real_t> >;
template class Field< Simplex , SurrealD >;

} // avro
