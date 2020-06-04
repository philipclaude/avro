#include "adaptation/metric.h"

#include "mesh/interpolation.h"

namespace avro
{

template<typename type,typename T>
FieldInterpolation<type,T>::FieldInterpolation( const Field<type,T>& fld ) :
  field_(fld),
  searcher_(field_.topology())
{}

template<typename type,typename T>
int
FieldInterpolation<type,T>::eval( const Points& points ,  index_t p , const std::vector<index_t>& guesses , T& tp )
{
  avro_assert( p < points.nb() );
  avro_assert( p >= points.nb_ghost() );
  const Topology<type>& topology = field_.topology();
  std::vector<real_t> phi( field_.master().nb_basis() , 0. );
  std::vector<real_t> xref( topology.master().number() + 1 );
  bool success;

  int ielem = -1;
  for (index_t iguess=0;iguess<guesses.size();iguess++)
  {
    ielem = searcher_.find( points[p] , guesses[iguess] );
    if (ielem<0)
    {
      // point is probably outside domain
  		// let's make sure by first brute forcing the check
      ielem = searcher_.brute( points[p] );
      if (ielem<0)
      {
        // get the reference coordinates of the closest element
  			ielem = searcher_.closest( points[p] , xref );
        topology.master().basis().evaluate( xref.data() , phi.data() );

        // perform the interpolation and return the element containing the point
  			index_t elem = index_t(ielem);
  			success = field_.dof().interpolate( field_[elem] , field_.nv(elem) , phi , &tp );
  			if (success) return ielem;
      }
    }

    // get the reference coordinates of the element and evaluate the basis functions for interpolation
    index_t elem = index_t(ielem);
    topology.master().physical_to_reference( topology.points() , topology(elem) , topology.nv(elem) , points[p] , xref.data() );
    topology.master().basis().evaluate( xref.data() , phi.data() );

    // perform the interpolation and return the element containing the point
    success = field_.dof().interpolate( field_[elem] , field_.nv(elem) , phi , &tp );
    if (success) return ielem;
  }
  return -1; // indicate there was a problem interpolating
}

template<typename type>
GeometryMetric<type>::GeometryMetric( const Field<type,real_t>& fld ) :
  FieldInterpolation<type,real_t>(fld)
{}

template<typename type>
int
GeometryMetric<type>::eval( const Points& points ,  index_t p , const std::vector<index_t>& guesses , real_t& mp )
{
  printf("evaluating analytically!\n");
  return -1;
}

template class FieldInterpolation<Simplex,real_t>;
template class FieldInterpolation<Simplex,Metric>;
template class GeometryMetric<Simplex>;

} // avro
