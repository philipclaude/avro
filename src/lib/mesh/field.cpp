#include "common/data.h"

#include "geometrics/egads.h"
#include "geometrics/primitive.h"
#include "geometrics/plc.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/field.h"
#include "mesh/topology.h"
#include "mesh/vertices.h"

#include "numerics/matrix.h"

namespace ursa
{

template<typename T>
FieldBase<T>::FieldBase( FieldType type ) :
  type_(type)
{}

template<typename Basis,typename T>
Field<Simplex<Basis>,T>::Field( const Topology<Shape_t>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type),
  topology_(topology),
  master_(topology.number(),order)
{
  printf("constructing field with order %u\n",master_.order());
}

template<typename Basis,typename T>
void
Field<Simplex<Basis>,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    // get the number of unique entries in the field
    // using the number of elements and nb_poly, for now assume p = 1
    ursa_assert( master_.order()==1 );

    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_.template Data<index_t>::operator()(k);
      const index_t nv = topology_.nv(k);
      Data<index_t>::add( v , nv );
    }

    for (index_t k=0;k<topology_.vertices().nb();k++)
    {
      this->data_.push_back( T(0) );
    }
  }
  else if (this->type()==DISCONTINUOUS)
  {
    ursa_assert( master_.order()==0 );
    const index_t nb_poly = 1; // assume order zero
    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_.template Data<index_t>::operator()(k);
      const index_t nv = topology_.nv(k);
      Data<index_t>::add( v , nv );
      for (index_t j=0;j<nb_poly;j++)
        this->data_.push_back( T(0) );
    }
  }
  else
    ursa_assert_not_reached;
}

template<typename T>
Field<Polytope,T>::Field( Topology<Polytope>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type),
  topology_(topology),
  master_(topology.number(),topology.order(),topology.master().incidence())
{
  printf("constructing field with order %u\n",master_.order());
}

template<typename T>
void
Field<Polytope,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    ursa_assert( master_.order()==1 );

    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_.template Data<index_t>::operator()(k);
      const index_t nv = topology_.nv(k);
      Data<index_t>::add( v , nv );
    }

    for (index_t k=0;k<topology_.vertices().nb();k++)
    {
      this->data_.push_back( T(0) );
    }
  }
  else if (this->type()==DISCONTINUOUS)
  {
    ursa_assert( master_.order()==0 );
    for (index_t k=0;k<topology_.nb();k++)
    {
      this->data_.push_back( T(0) );
    }
  }
  else
    ursa_assert_not_reached;
}

template class Field< Simplex<Lagrange> , real_t >;
template class Field< Simplex<Bezier> , real_t >;
template class Field< Polytope , real_t >;

template class Field< Simplex<Lagrange> , std::vector<real_t> >;
template class Field< Simplex<Lagrange> , std::vector<index_t> >;

template class Field< Simplex<Lagrange> , geometrics::Primitive* >;
template class Field< Simplex<Lagrange> , numerics::SymMatrixD<real_t> >;

} // ursa
