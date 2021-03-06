//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_NUMERICS_INTEGRATION_H_
#define avro_LIB_NUMERICS_INTEGRATION_H_

#include "common/error.h"
#include "common/parallel_for.h"
#include "avro_types.h"

#include "mesh/dof.h"
#include "mesh/field.h"

#define QUAD_ORDER 4

namespace avro
{

template<typename> class Topology;

template<typename Derived>
struct Integrand
{
  Integrand() {}
  Integrand(const Integrand&) = delete;
  Integrand operator=(const Integrand&) = delete;
  const Derived& cast() const { return static_cast<const Derived&>(*this); }
};

template<typename type>
class Integrand_Volume : public Integrand<Integrand_Volume<type>>
{
public:
  typedef real_t T;

public:
  Integrand_Volume( const Topology<type>& topology ) :
    topology_(topology)
  {}

  T operator()( index_t k , const real_t* xref , const real_t* x ) const
  {
    return 1.0;
  }

private:
  const Topology<type>& topology_;
};

template<typename type, typename _T,typename Functor>
class Integrand_Field : public Integrand<Integrand_Field<type,_T,Functor>> {

public:

  typedef _T T;

  Integrand_Field( const Field<type,T>& field ) :
    field_(field)
  {}

  T operator()( index_t k , const QuadraturePoint& point , const real_t* x ) const {

    const matd<real_t>& Bu = field_.element().reference().quadrature().get_basis( field_.element().number() , field_.element().order() , QUAD_ORDER );

    std::vector<real_t> phi( field_.element().nb_basis(), 0 );
    std::vector<real_t> phix,phixx;

    std::vector<T> u(field_.dof().rank());
    std::vector<T> ux,uxx;

    for (index_t j = 0; j < phi.size(); j++)
      phi[j] = Bu(j,point.index());

    const real_t* xref = point.coordinate();

    if (functor_.needs_solution()) {
      field_.dof().interpolate( field_[k] , field_.nv(k) , phi , u.data() );
    }
    if (functor_.needs_gradient()) {
      avro_implement;
      field_.element().reference().basis().evaluate( xref , phix.data() );
      //field_.dof().interpolate( field_(k) , field_.nv(k) , phix , ux.data() );
    }
    if (functor_.needs_hessian()) {
      avro_implement;
      field_.element().reference().basis().evaluate( xref , phixx.data() );
      //field_.dof().interpolate( field_(k) , field_.nv(k) , phixx , uxx.data() );
    }

    // evaluate the function
    return functor_(x,u,ux,uxx);
  }

private:
  const Field<type,T>& field_;
  Functor functor_;
};

class ElementIntegralBase {};

template<typename type,typename T>
class ElementIntegral : public ElementIntegralBase
{
public:
  ElementIntegral( const Topology<type>& topology , const DOF<T>& dof,
           const Table<index_t>& field , const type& element ) :
    topology_(topology),
    field_(field),
    element_(element),
    dof_(dof)
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
  integrate( index_t k , const Integrand& integrand , T& f )
  {
    get(k);

    f = 0;
    std::vector<real_t> x(topology_.points().dim());
    std::vector<real_t> phix( topology_.element().nb_basis() );

    const QuadratureStore<type>& quadrature_x = topology_.element().reference().quadrature();
    const matd<real_t>& Bx = quadrature_x.get_basis( topology_.element().number() , topology_.element().order() , QUAD_ORDER );
    const Quadrature& quad = quadrature_x.quadrature( element_.number() , QUAD_ORDER );
    QuadraturePoint point(quad);

    for (index_t i = 0; i < quad.nb(); i++) {

      point.set_index(i);
      real_t w = point.weight();

      // retrieve the basis functions evaluated at this quadrature point
      for (index_t j = 0; j < phix.size(); j++)
        phix[j] = Bx(j,i);

      // evaluate the physical coordinates
      topology_.points().interpolate( x_ , phix , x.data() );

      // evaluate the jacobian at the reference point (for now assume constant)
      real_t dj = topology_.element().jacobian(x_,topology_.points().dim());

      // evaluate the integrand at the quadrature point
      f += w*integrand( k , point , x.data() )*dj;
    }
  }

private:
  const Topology<type>& topology_;
  const Table<index_t>& field_;
  const type& element_;
  const DOF<T>& dof_;

  std::vector<const real_t*> x_;
  std::vector<const T*>      f_;

};

template<typename Integrand>
class Functional
{
private:
  typedef typename Integrand::T T;
  typedef Functional thisclass;

public:
  Functional( const Integrand& integrand ) :
    integrand_(integrand),
    functional_(0)
  {}

  template<typename type>
  Functional( const Integrand& integrand , const Topology<type>& topology ) :
    integrand_(integrand),
    functional_(0)
  {
    initialize(topology);
  }

  template<typename type>
  void
  initialize( const Topology<type>& topology )
  {
    integrators_.resize( topology.nb() );
    for (index_t k = 0; k < topology.nb(); k++)
      integrators_[k] = std::make_shared<ElementIntegral<type,T>>(topology,topology.points(),topology,topology.element());
  }

  template<typename type>
  void
  integrate_element( index_t k )
  {
    ElementIntegral<type,T>& elem = static_cast<ElementIntegral<type,T>&>(*integrators_[k].get());
    T df = T(0);
    elem.integrate( k , integrand_ , df );
    values_[k] = df;
  }

  template<typename type>
  void
  integrate( const Topology<type>& topology , T* values=nullptr )
  {
    if (integrators_.size() != topology.nb() ) initialize(topology);
    values_.clear();
    values_.resize( topology.nb() );
    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::integrate_element<type> ),
      0,values_.size() );
    for (index_t k = 0; k < values_.size(); k++)
    {
      functional_ += values_[k];
      if (values != nullptr) values[k] = values_[k];
    }
  }

  template<typename type>
  void
  integrate( const Field<type,T>& field )
  {
    const Topology<type>& topology = field.topology();
    avro_assert( topology.nb() == field.nb() );

    ElementIntegral<type,T> elem( topology , field.dof() , field , field.element() );
    for (index_t k=0;k<field.nb_elem();k++)
    {
      T df = 0;
      elem.integrate( k , integrand_ , df );
      functional_ += df;
    }
  }

  T value() const { return functional_; }

private:
  const Integrand& integrand_;
  std::vector<T> values_;
  T functional_;
  std::vector<std::shared_ptr<ElementIntegralBase>> integrators_;
};

#if 0
template<typename Integrand,typename type>
class Functional
{
private:
  typedef typename Integrand::T T;
  typedef Functional_typed thisclass;

public:
  Functional_typed( const Integrand& integrand ) :
    integrand_(integrand),
    functional_(0),
    topology_(nullptr)
  {}

  void
  integrate( const Topology<type>& topology , T* values=nullptr )
  {
    #if 0
    ElementIntegral<type,T> elem( topology , topology.points() , topology , topology.element() );
    //#pragma omp parallel for
    for (index_t k=0;k<topology.nb();k++)
    {
      T df = T(0);
      elem.integrate( k , integrand_ , df );
      functional_ += df;
      if (values != nullptr) values[k] = df;
    }
    #else
    values_.clear();
    values_.resize( topology.nb() );

    // initialize the element integrals

    ProcessCPU::parallel_for(
      parallel_for_member_callback( this , &thisclass::integrate_element ),
      0,values_.size() );
    #endif
  }

  T value() const { return functional_; }

private:
  const Integrand& integrand_;
  std::vector<T> values_;
  T functional_;
  const Topology<type>& topology_;
};
#endif


} // avro

#endif
