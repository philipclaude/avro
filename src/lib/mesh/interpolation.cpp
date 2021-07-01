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

#include "mesh/interpolation.h"
#include "mesh/search.h"

namespace avro
{

template<typename type,typename T>
FieldInterpolation<type,T>::FieldInterpolation( const Field<type,T>* fld ) :
  analytic_(false),
  pfield_(fld)
{
  if (fld!=nullptr)
    searcher_ = std::make_shared<ElementSearch<type>>(pfield_->topology());
}

template<typename type,typename T>
int
FieldInterpolation<type,T>::eval( const Points& points ,  index_t p , const std::vector<index_t>& guesses , T& tp )
{
  avro_assert( p < points.nb() );
  avro_assert( p >= points.nb_ghost() );
  avro_assert( pfield_!=nullptr );
  const Field<type,T>& field_ = *pfield_;
  const Topology<type>& topology = field_.topology();
  std::vector<real_t> phi( field_.element().nb_basis() , 0. );
  std::vector<real_t> xref( topology.element().number() + 1 );
  bool success;

  int ielem = -1;
  for (index_t iguess=0;iguess<guesses.size();iguess++)
  {
    ielem = searcher_->find( points[p] , guesses[iguess] );
    if (ielem<0)
    {
      // point is probably outside domain
  		// let's make sure by first brute forcing the check
      ielem = searcher_->brute( points[p] );
      if (ielem<0)
      {
        // get the reference coordinates of the closest element
  			ielem = searcher_->closest( points[p] , xref );
        topology.element().basis().evaluate( xref.data() , phi.data() );

        // perform the interpolation and return the element containing the point
  			index_t elem = index_t(ielem);
  			success = field_.dof().interpolate( field_[elem] , field_.nv(elem) , phi , &tp );
  			if (success) return ielem;
      }
    }

    // get the reference coordinates of the element and evaluate the basis functions for interpolation
    index_t elem = index_t(ielem);
    topology.element().physical_to_reference( topology.points() , topology(elem) , topology.nv(elem) , points[p] , xref.data() );
    topology.element().basis().evaluate( xref.data() , phi.data() );

    // perform the interpolation and return the element containing the point
    success = field_.dof().interpolate( field_[elem] , field_.nv(elem) , phi , &tp );

    if (success) return ielem;
  }
  return -1; // indicate there was a problem interpolating
}

template class FieldInterpolation<Simplex,real_t>;
template class FieldInterpolation<Simplex,Metric>;

} // avro
