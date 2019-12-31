#ifndef avro_LIB_NUMERICS_INTEGRATION_H_
#define avro_LIB_NUMERICS_INTEGRATION_H_

#include "common/error.h"
#include "common/types.h"

#include "mesh/dof.h"
#include "mesh/field.h"

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
class Integrand_Field : public Integrand<Integrand_Field<type,_T,Functor>>
{
public:

  typedef _T T;

  Integrand_Field( const Field<type,T>& field ) :
    field_(field)
  {}

  T operator()( index_t k , const real_t* xref , const real_t* x ) const
  {
    std::vector<real_t> phi( field_.master().nb_basis(), 0 );
    std::vector<real_t> phix,phixx;

    std::vector<T> u(field_.dof().rank());
    std::vector<T> ux,uxx;

    if (functor_.needs_solution())
    {
      field_.master().basis().evaluate( xref , phi.data() );
      field_.dof().interpolate( field_[k] , field_.nv(k) , phi , u.data() );
    }
    if (functor_.needs_gradient())
    {
      field_.master().basis().evaluate( xref , phix.data() );
      //field_.dof().interpolate( field_(k) , field_.nv(k) , phix , ux.data() );
    }
    if (functor_.needs_hessian())
    {
      field_.master().basis().evaluate( xref , phixx.data() );
      //field_.dof().interpolate( field_(k) , field_.nv(k) , phixx , uxx.data() );
    }

    // evaluate the function
    return functor_(x,u,ux,uxx);
  }

private:
  const Field<type,T>& field_;
  Functor functor_;
};

template<typename type,typename T>
class ElementIntegral
{
public:
  ElementIntegral( const Topology<type>& topology , const DOF<T>& dof,
           const Table<index_t>& field , const type& master ) :
    topology_(topology),
    field_(field),
    master_(master),
    dof_(dof)
  {
    if (topology_.layout()==TableLayout_Rectangular )
    {
      avro_assert( field.layout()==TableLayout_Rectangular );
      avro_assert( master_.nb_basis()==field_.rank() );
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
    std::vector<real_t> phi( topology_.master().nb_basis() );
    for (index_t i=0;i<master_.nb_quad();i++)
    {
      // retrieve the quadrature point and weight
      real_t w = master_.quad_weight(i);
      const real_t* xref = master_.quad_point(i);

      // evaluate the basis functions at the quadrature point
      topology_.master().basis().evaluate( xref , phi.data() );

      // evaluate the physical coordinates
      topology_.points().interpolate( x_ , phi , x.data() );

      // evaluate the jacobian at the reference point (for now assume constant)
      real_t dj = topology_.master().jacobian(x_,topology_.points().dim());

      // evaluate the integrand at the quadrature point
      f += w*integrand( k , xref , x.data() )*dj;
    }
  }

private:
  const Topology<type>& topology_;
  const Table<index_t>& field_;
  const type& master_;
  const DOF<T>& dof_;

  std::vector<const real_t*> x_;
  std::vector<const T*>      f_;

};

template<typename Integrand>
class Functional
{
private:
  typedef typename Integrand::T T;

public:
  Functional( const Integrand& integrand ) :
    integrand_(integrand),
    functional_(0)
  {}

  template<typename type>
  void
  integrate( const Topology<type>& topology )
  {
    ElementIntegral<type,T> elem( topology , topology.points() , topology , topology.master() );
    for (index_t k=0;k<topology.nb();k++)
    {
      T df = 0;
      elem.integrate( k , integrand_ , df );
      functional_ += df;
    }
  }

  template<typename type>
  void
  integrate( const Field<type,T>& field )
  {
    const Topology<type>& topology = field.topology();
    avro_assert( topology.nb() == field.nb() );

    ElementIntegral<type,T> elem( topology , field.dof() , field , field.master() );
    for (index_t k=0;k<field.nb();k++)
    {
      T df = 0;
      elem.integrate( k , integrand_ , df );
      functional_ += df;
    }
  }

  T value() const { return functional_; }

private:
  const Integrand& integrand_;
  T functional_;

};

} // avro

#endif
