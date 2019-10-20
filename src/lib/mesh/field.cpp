#include "common/data.h"

#include "geometrics/egads.h"
#include "geometrics/primitive.h"
#include "geometrics/plc.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/builder.h"
#include "mesh/builder.hpp"
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

template<typename ShapeBasis_t,typename FieldBasis_t,typename T>
Field<Simplex<ShapeBasis_t>,Simplex<FieldBasis_t>,T>::Field( const Topology<Simplex<ShapeBasis_t>>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type),
  topology_(topology),
  master_(topology.number(),order)
{
  printf("constructing field with order %u\n",master_.order());
}

template<typename T>
Field<Polytope,Polytope,T>::Field( Topology<Polytope>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type),
  topology_(topology),
  master_(topology,order)
{
  printf("constructing field with order %u\n",master_.order());
}

template<typename ShapeBasis_t,typename FieldBasis_t,typename T>
void
Field<Simplex<ShapeBasis_t>,Simplex<FieldBasis_t>,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    // get the number of unique entries in the field
    // using the number of elements and nb_poly, for now assume p = 1
    ursa_assert( master_.order()==1 );

    Builder<Simplex<ShapeBasis_t>,Simplex<FieldBasis_t>> builder(topology_,master_);
    builder.template transfer<T>(*this);
    return;

    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_(k);
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
    printf("heerrr\n");
    ursa_assert( master_.order()==0 );
    const index_t nb_poly = 1; // assume order zero
    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_(k);
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
void
Field<Polytope,Polytope,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    ursa_assert( master_.order()==1 );

    for (index_t k=0;k<topology_.nb();k++)
    {
      const index_t* v = topology_(k);
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

template class Field< Simplex<Lagrange> , Simplex<Lagrange> , real_t >;
template class Field< Simplex<Lagrange> , Simplex<Bezier> , real_t >;
template class Field< Polytope , Polytope , real_t >;

template class Field< Simplex<Lagrange> , Simplex<Lagrange> , std::vector<real_t> >;
template class Field< Simplex<Lagrange> , Simplex<Lagrange> , std::vector<index_t> >;

template class Field< Simplex<Lagrange> , Simplex<Lagrange> , geometrics::Primitive* >;
template class Field< Simplex<Lagrange> , Simplex<Lagrange> , numerics::SymMatrixD<real_t> >;

} // ursa
