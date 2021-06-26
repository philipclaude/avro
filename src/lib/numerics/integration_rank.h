//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_INTEGRATION_RANK_H_
#define avro_LIB_NUMERICS_INTEGRATION_RANK_H_

#include "common/error.h"
#include "types.h"

#include "mesh/dof.h"
#include "mesh/field.h"

#include "numerics/integration.h"

namespace avro
{

template<typename> class Topology;

template<typename type,typename T>
class ElementIntegral_Ranked : public ElementIntegralBase
{
public:
  ElementIntegral_Ranked( const Topology<type>& topology , const DOF<T>& dof,
           const Table<index_t>& field , const type& element , index_t nrank ) :
    topology_(topology),
    field_(field),
    element_(element),
    dof_(dof),
    nrank_(nrank)
  {
    if (topology_.layout()==TableLayout_Rectangular )
    {
      avro_assert( field.layout()==TableLayout_Rectangular );
      avro_assert( element_.nb_basis()==field_.rank() );
      x_.resize( topology_.rank() );
      f_.resize( field_.rank() );
    }
    else
      avro_implement
  }

  void get( index_t k )
  {
    if (topology_.layout()==TableLayout_Rectangular )
    {
      for (index_t j=0;j<field_.nv(k);j++)
        f_[j] = dof_[field_(k,j)];
      for (index_t j=0;j<topology_.nv(k);j++)
        x_[j] = topology_.points()[topology_(k,j)];
    }
    else
    {
      // need to query how many dof/points we need to store
      avro_implement;
    }
  }

  template<typename Integrand>
  void
  integrate( index_t k , const Integrand& integrand , std::vector<T>& f )
  {
    get(k);

    std::vector<T> I(nrank_,0);
    avro_assert( f.size() == nrank_ );

    std::fill( f.begin() , f.end() , 0 );
    std::vector<real_t> x(topology_.points().dim());
    std::vector<real_t> phi( topology_.element().nb_basis() );
    for (index_t i=0;i<element_.nb_quad();i++)
    {
      // retrieve the quadrature point and weight
      real_t w = element_.quad_weight(i);
      const real_t* xref = element_.quad_point(i);

      // evaluate the basis functions at the quadrature point
      topology_.element().basis().evaluate( xref , phi.data() );

      // evaluate the physical coordinates
      topology_.points().interpolate( x_ , phi , x.data() );

      // evaluate the jacobian at the reference point (for now assume constant)
      real_t dj = topology_.element().jacobian(x_,topology_.points().dim());

      // evaluate the integrand at the quadrature point
      integrand( k , xref , x.data() , I );

      for (index_t r = 0; r < nrank_; r++)
        f[r] += w*I[r]*dj;

      //avro_assert( dj > 0.0 );
    }
  }

private:
  const Topology<type>& topology_;
  const Table<index_t>& field_;
  const type& element_;
  const DOF<T>& dof_;

  std::vector<const real_t*> x_;
  std::vector<const T*>      f_;
  index_t nrank_;

};

template<typename Integrand>
class Functional_Ranked
{
private:
  typedef typename Integrand::T T;
  typedef Functional_Ranked thisclass;

public:
  Functional_Ranked( const Integrand& integrand , index_t nrank ) :
    integrand_(integrand),
    functional_(nrank,0),
    nrank_(nrank)
  {}

  template<typename type>
  Functional_Ranked( const Integrand& integrand , index_t nrank , const Topology<type>& topology ) :
    integrand_(integrand),
    functional_(nrank,0),
    nrank_(nrank)
  {
    integrators_.resize( topology.nb() );
    for (index_t k = 0; k < topology.nb(); k++)
      integrators_[k] = std::make_shared<ElementIntegral_Ranked<type,T>>(topology,topology.points(),topology,topology.element(),nrank_);
  }

  template<typename type>
  void
  integrate_element( index_t k )
  {
    ElementIntegral_Ranked<type,T>& elem = static_cast<ElementIntegral_Ranked<type,T>&>(*integrators_[k].get());
    std::vector<T> df(nrank_,0);
    elem.integrate( k , integrand_ , df );
    for (index_t r = 0; r < nrank_; r++)
      values_[k*nrank_+r] = df[r];
  }

  template<typename type>
  void
  integrate( const Topology<type>& topology , T* values=nullptr )
  {
    #if 0
    ElementIntegral_Ranked<type,T> elem( topology , topology.points() , topology , topology.element() , nrank_ );
    for (index_t k=0;k<topology.nb();k++)
    {
      std::vector<T> df(nrank_,0);
      elem.integrate( k , integrand_ , df );
      for (index_t r = 0; r < nrank_; r++)
        functional_[r] += df[r];
      if (values != nullptr)
      {
        for (index_t r = 0; r < nrank_; r++)
          values[k*nrank_+r] = df[r];
      }
    }
    #else
    avro_assert( topology.nb() == integrators_.size() );
    values_.clear();
    values_.resize( topology.nb()*nrank_ );
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::integrate_element<type> ),
      0,topology.nb() );
    for (index_t k = 0; k < topology.nb(); k++)
    {
      for (index_t r = 0; r < nrank_; r++)
        functional_[r] += values_[k*nrank_+r];
    }
    if (values != nullptr)
    {
      for (index_t i = 0; i < values_.size(); i++)
        values[i] = values_[i];
    }

    #endif
  }

  std::vector<T> value() const { return functional_; }

private:
  const Integrand& integrand_;
  std::vector<real_t> functional_;
  index_t nrank_;
  std::vector<std::shared_ptr<ElementIntegralBase>> integrators_;
  std::vector<T> values_;

};

} // avro

#endif
