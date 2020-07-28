#ifndef AVRO_LIB_SIMULATION_EQUATIONS_H_
#define AVRO_LIB_SIMULATION_EQUATIONS_H_

#include "simulation/equation_traits.h"
#include "simulation/jacobian.h"

#include <numpack/AlgebraicEquationSetBase.h>

class Residual;
class Jacobian;

namespace avro
{

class IntegrateCellGroups
{
    template<class IntegralCell,class Shape_t,class Field_t>
    static void
    integrate( CellIntegralType<IntegralCell>&& integral_type )
    {
      IntegralCell& integral = integral_type.cast();
      avro_implement;

      integral.template integrate<Shape_t>();
    }
};

template<typename PDE>
class EquationSet_LineSearch
{
public:

  EquationSet_LineSearch( const PDE& pde , const std::vector<real_t>& tol ) :
    pde_(pde),
    tol_(tol)
  {}

  virtual ~EquationSet_LineSearch() {}

  virtual void dump_solution( const std::string& filename ) { avro_assert_not_reached; }

  virtual bool residual_converged( const std::vector<std::vector<real_t>>& rsd ) const
  {
    avro_implement;
    return false;
  }

  virtual bool residual_converged( const std::vector<std::vector<real_t>>& rsd , index_t ieq , index_t imon ) const
  {
    avro_implement;
    return false;
  }

  virtual bool residual_decreased( const std::vector<std::vector<real_t>>& rsd0 ,
                                  const std::vector<std::vector<real_t>>& rsd1 ) const
  {
    avro_implement;
    return false;
  }

private:
  const PDE& pde_;
  const std::vector<real_t> tol_;
};

template<typename PDE,class Traits>
class EquationSet_Galerkin : public EquationSet_LineSearch<PDE>
{
public:
  typedef typename PDE::template ArrayQ<real_t> ArrayQ;
  typedef typename PDE::template MatrixQ<real_t> MatrixQ;

  typedef EquationSetTraits<MatrixQ,ArrayQ,Traits> Traits_t;
  typedef typename Traits_t::AlgebraicEquationSetBaseClass BaseType;

  typedef typename Traits_t::VectorSizeClass VectorSizeClass;
  typedef typename Traits_t::MatrixSizeClass MatrixSizeClass;

  typedef typename Traits_t::SystemMatrix SystemMatrix;
  typedef typename Traits_t::SystemVector SystemVector;

  typedef typename Traits_t::SystemMatrixView SystemMatrixView;
  typedef typename Traits_t::SystemVectorView SystemVectorView;

  typedef IntegrandCell<PDE> IntegrandCellClass;

  virtual void residual( SystemVectorView& rsd ) override;
  virtual void jacobian( SystemMatrixView& mtx ) override { this->template jacobian< SystemMatrixView&>(mtx); }

  static const int iPDE = 0;
  static const int iBC  = 1;
  static const int iq   = 0;
  static const int ilg  = 1;

protected:
  template<class SparseMatrix_t> void jacobian( SparseMatrix_t mtx );

  IntegrandCellClass fcn_cell_;
  Field<Simplex,ArrayQ>& qfld_;
  Field<Simplex,ArrayQ>& lgfld_;
  const PDE& pde_;
  const std::vector<real_t> tol_;
};

template<typename PDE,class Traits_t>
template<typename SparseMatrix_t>
void
EquationSet_Galerkin<PDE,Traits_t>::jacobian( SparseMatrix_t jac )
{
  avro_assert( jac.m() == 2 );
  avro_assert( jac.n() == 2 );

  typedef typename std::result_of<SparseMatrix_t(const int,const int)>::type Matrix;

  Matrix jacPDE_q  = jac(iPDE,iq);
  Matrix jacPDE_lg = jac(iPDE,ilg);
  Matrix jacBC_q   = jac(iBC,iq);
  Matrix jacBC_lg  = jac(iBC,ilg);

  typedef SurrealS<PDE::N> SurrealClass;

  const int qorder = -1;

  // dispatch the integration over the interior cells
  IntegrateCellGroups::integrate( JacobianCell_impl<IntegrandCellClass>(fcn_cell_,jacPDE_q) );//, qfld_.topology() , qfld_ , qorder ) );

  // dispatch the interation over the boundary traces
  // TODO
}


} // avro

#endif
