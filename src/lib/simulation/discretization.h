#ifndef AVRO_SRC_LIB_DISCRETIZATION_H_
#define AVRO_SRC_LIB_DISCRETIZATION_H_

#include "common/error.h"
#include "common/types.h"

#include <numpack/dense/tools/index.h>

namespace avro
{

template<typename> class FieldBase;
template<typename,typename> class Field;
class Equations;

template<typename derived_t>
class Discretization
{
public:
  Discretization();

  index_t nb() const { return derived().nb(); }

  derived_t& derived() { return static_cast<derived_t&>(*this); }
  const derived_t& derived() const { return static_cast<const derived_t&>(*this); }

  template<typename shape_t,typename field_t>
  void assemble( const Field<shape_t,field_t>& field )
  {
    for (index_t k=0;k<nb();k++)
      derived().assemble(field,k);
  }

};

template<typename PDE_t>
class IntegrandCell
{
public:
  static const int N_ = PDE_t::N;

  // solution, gradient and jacobian arrays
  template<typename T> using ArrayQ       = typename PDE_t::template ArrayQ<T>;
  template<typename T> using VectorArrayQ = typename PDE_t::template VectorArrayQ<T>;
  template<typename T> using MatrixQ      = typename PDE_t::template MatrixQ<T>;

  IntegrandCell( const PDE_t& pde ) :
    pde_(pde)
  {}

  template<typename T,typename type>
  class BasisWeighted
  {
  public:
    BasisWeighted( const PDE_t& pde , const Field<T,type>& field ) :
      pde_(pde),
      field_(field),
      phi_(field.element().nb_basis()),
      dphi_(field.element().nb_basis())
    {}

    bool needs_evaluation() const;
    coord_t nb_pde() const { return N_; }
    index_t nb_dof() const { return phi_.size(); }

    // residual integrand
    template<class Ti> void operator()( const real_t dJ , const real_t* xref , numpack::DLA::VectorD<ArrayQ<Ti> >& elemR ) const;

    // jacobian integrand
    void operator()( const real_t dJ , const real_t* xref , numpack::DLA::MatrixD<MatrixQ<real_t>>& elemJ ) const;

  protected:
    template<class Tq, class Tg, class Ti>
    void weighted_integrand( const ArrayQ<Tq>& q , const VectorArrayQ<Tg>& gradq , ArrayQ<Ti> integrand[] , int neqn ) const;

    const PDE_t& pde_;
    const Field<T,type>& field_;
    std::vector<real_t> phi_;
    std::vector<real_t> dphi_;
  };

protected:
  const PDE_t& pde_;
};


template<typename PDE_t>
template<typename T,typename type>
template<typename Ti>
void
IntegrandCell<PDE_t>::BasisWeighted<T,type>::operator()( const real_t dJ , const real_t* xref , numpack::DLA::VectorD<ArrayQ<Ti>>& elemR ) const
{
  avro_assert( elemR.m() == nb_dof() );

  ArrayQ<T> q;        // solution
  VectorArrayQ<T> dq; // gradient

  const bool needs_gradient = pde_.has_flux_viscous() || (pde_.has_source() && pde_.needs_source_gradient() );

  // evaluate elemental parameters
  // TODO

  // evaluate the basis functions and gradients
  field_.element().eval( xref , phi_.data() );
  field_.element().eval_gradient( xref, dphi_.data() );

  // evaluate the solution and solution gradient
  // TODO
  if (needs_gradient)
  {
    // TODO
  }

  // compute the residual
  numpack::DLA::VectorD<ArrayQ<T>> integrand( nb_dof() );
  weighted_integrand( q , dq , integrand.data() , integrand.size() );

  // accumulate the weighted residual
  for (index_t k=0;k<nb_dof();k++)
    elemR[k] += dJ*integrand(k);
}

template<typename PDE_t>
template<typename T,typename type>
void
IntegrandCell<PDE_t>::BasisWeighted<T,type>::operator()( const real_t dJ , const real_t* xref , numpack::DLA::MatrixD<MatrixQ<real_t>>& elemJ ) const
{
  typedef SurrealS<PDE_t::N> SurrealClass;

  avro_assert( elemJ.m() == nb_dof() );
  avro_assert( elemJ.n() == nb_dof() );

  static const int nb_deriv = SurrealClass::N;
  static_assert( nb_deriv & PDE_t::N == 0 , "derivatives must be a multiple of PDE::N" );

  ArrayQ<real_t> q;        // solution
  VectorArrayQ<real_t> dq; // gradient

  ArrayQ<SurrealClass> qs;        // solution
  VectorArrayQ<SurrealClass> dqs; // gradient

  MatrixQ<real_t> PDE_q = 0;
  MatrixQ<real_t> PDE_dq = 0;

  const bool needs_gradient = pde_.has_flux_viscous() || (pde_.has_source() && pde_.needs_source_gradient() );

  // evaluate elemental parameters
  // TODO

  // evaluate the basis functions and gradients
  field_.element().eval( xref , phi_.data() );
  field_.element().eval_gradient( xref, dphi_.data() );

  // evaluate the solution and solution gradient
  // TODO
  qs  = q;
  if (needs_gradient)
  {
    // TODO
    dqs = dq;
  }

  numpack::DLA::VectorD<ArrayQ<SurrealClass>> integrand_surreal( nb_dof() );

  for (index_t nchunk=0;nchunk<PDE_t::N;nchunk+=nb_deriv)
  {
    // associate derivative slots with solution variables
    index_t slot = 0;
    if ((slot >= nchunk) && (slot < nchunk + nb_deriv))
      for (index_t n=0;n<PDE_t::N;n++)
        numpack::DLA::index(qs,n).deriv(slot+n-nchunk) = 1;
    slot += PDE_t::N;

    // evaluate the integrand
    integrand_surreal = 0;
    weighted_integrand( qs , dq , integrand_surreal.data() , integrand_surreal.size() );

    // accumulate the derivatives into the element jacobian
    slot = 0;
    if ((slot >= nchunk) && (slot < nchunk + nb_deriv))
    {
      for (index_t n=0;n<PDE_t::N;n++)
        numpack::DLA::index(qs,n).deriv(slot+n-nchunk) = 0; // reset derivative

      for (index_t i=0;i<nb_dof();i++)
      {
        for (index_t m=0;m<PDE_t::N;m++)
        for (index_t n=0;n<PDE_t::N;n++)
          numpack::DLA::index(PDE_q,m,n) = numpack::DLA::index(integrand_surreal[i],m).deriv(slot+n-nchunk);

        for (index_t j=0;j<nb_dof();j++)
          elemJ(i,j) += dJ*phi_[j]*PDE_q;
      }
    }
  }
}

template<typename PDE_t>
template<typename T,typename type>
template<class Tq,class Tg,class Ti>
void
IntegrandCell<PDE_t>::BasisWeighted<T,type>::weighted_integrand( const ArrayQ<Tq>& q , const VectorArrayQ<Tg>& dq , ArrayQ<Ti> integrand[] , int neqn ) const
{
  VectorArrayQ<Ti> F;
  ArrayQ<Ti> source;

  for (index_t k=0;k<neqn;k++)
    integrand[k] = 0;

  if (pde_.has_flux_advective() || pde_.has_flux_viscous())
  {
    F = 0;

    // advective flux
    if (pde_.has_flux_advective())
      pde_.get_flux_advective( q , F );

    if (pde_.has_viscous_flux())
      pde_.get_flux_viscous( q, dq , F );

    for (index_t k=0;k<neqn;k++)
      integrand[k] -= dot(dphi_[k],F);
  }

  if (pde_.has_source())
  {
    source = 0;
    pde_.eval_source( q , dq , source );

    for (index_t k=0;k<neqn;k++)
      integrand[k] += phi_[k]*source;
  }

  if (pde_.has_forcing_function())
  {
    ArrayQ<real_t> forcing = 0;
    pde_.eval_forcing( forcing );

    for (index_t k=0;k<neqn;k++)
      integrand[k] -= phi_[k]*forcing;
  }
}

template<typename discretization_t,typename pde_t>
class Assembler
{
public:
  typedef typename pde_t::field_t field_t;

  Assembler( const FieldBase<field_t>& field , Equations& equations );

  void assemble();

private:

  template<typename shape_t> const Field<shape_t,field_t>& field_cast()
    { return static_cast<const Field<shape_t,field_t>&>(field_); }

  Discretization<discretization_t> discretization_;
  const FieldBase<field_t>& field_;
  Equations& equations_;
};

class CG : public Discretization<CG>
{

public:
  void assemble();

  template<typename shape_t,typename field_t>
  void assemble( const Field<shape_t,field_t>& field , index_t k )
  {
    avro_implement;
  }

  index_t nb() const { avro_implement; return 0; }
};

class DG : public Discretization<DG>
{

public:
  void assemble();
  void assemble( index_t k ); // k is an element in the mesh

  index_t nb() const { avro_implement; return 0; }
};

class FV : public Discretization<FV>
{

public:
  void assemble();
  void assemble( index_t k ); // k is an element in the mesh

  index_t nb() const { avro_implement; return 0; }
};

class FD : public Discretization<FD>
{
public:
  void assemble();
  void assemble( index_t k ); // k is a vertex in the mesh

  index_t nb() const { avro_implement; return 0; }
};

} // avro

#endif
