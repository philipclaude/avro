#include "adaptation/metric.h"

#include "geometry/entity.h"

#include "master/polytope.h"
#include "master/simplex.h"

#include "mesh/builder.h"
#include "mesh/builder.hpp"
#include "mesh/field.h"
#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/matrix.h"

namespace luma
{

template<typename T>
FieldBase<T>::FieldBase( FieldType type , TableLayoutCategory category ) :
  Table<index_t>(category),
  data_(1),
  type_(type)
{}

template<typename T>
Field<Simplex,T>::Field( const Topology<Simplex>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type,TableLayout_Rectangular),
  topology_(topology),
  master_(topology.number(),order)
{
  Table<index_t>::set_rank( master_.nb_basis() );
  //build();
}

template<typename T>
Field<Polytope,T>::Field( Topology<Polytope>& topology , coord_t order , FieldType type ) :
  FieldBase<T>(type),
  topology_(topology),
  master_(topology,order,topology.points().incidence())
{}

template<typename T>
void
Field<Simplex,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    // get the number of unique entries in the field
    // using the number of elements and nb_poly, for now assume p = 1
    Builder<Simplex> builder(topology_,master_.order(),BasisFunctionCategory_None);
    builder.template transfer<T>(*this);
  }
  else if (this->type()==DISCONTINUOUS)
  {
    index_t n = 0;
    std::vector<index_t> dof( master_.nb_basis() );
    for (index_t k=0;k<topology_.nb();k++)
    {
      for (index_t j=0;j<master_.nb_basis();j++)
        dof[j] = n++;
      Table<index_t>::add( dof.data() , dof.size() );
      const std::vector<index_t>& idx = Table<index_t>::data();
      index_t nb_dof = * std::max_element( idx.begin() , idx.end() ) +1;
      this->allocate( nb_dof );
    }
  }
  else
    luma_assert_not_reached;
}

template<typename T>
void
Field<Polytope,T>::build()
{
  if (this->type()==CONTINUOUS)
  {
    luma_assert( master_.order()==1 );

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
    luma_assert( master_.order()==0 );
    T x0(0);
    for (index_t k=0;k<topology_.nb();k++)
    {
      this->data_.add( &x0 , 1 );
    }
  }
  else
    luma_assert_not_reached;
}

template class Field< Simplex , real_t >;
template class Field< Polytope , real_t >;

template class Field< Simplex , std::vector<real_t> >;
template class Field< Simplex , std::vector<index_t> >;

template class Field< Simplex , Entity* >;
template class Field< Simplex , numerics::SymMatrixD<real_t> >;

template class Field< Simplex , Metric >;

} // luma
