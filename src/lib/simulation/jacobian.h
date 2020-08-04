#ifndef AVRO_LIB_JACOBIAN_H_
#define AVRO_LIB_JACOBIAN_H_

#include "common/types.h"

namespace avro
{

template<class Derived>
struct CellIntegralType
{
  CellIntegralType() {}
  CellIntegralType(const CellIntegralType&&) {}
  CellIntegralType(const CellIntegralType&) = delete;
  const CellIntegralType& operator=(const CellIntegralType&) = delete;

  //A convenient method for casting to the derived type
  inline const Derived& cast() const { return static_cast<const Derived&>(*this); }
  inline       Derived& cast()       { return static_cast<      Derived&>(*this); }
};

template<typename type,typename T>
class FieldElement
{
public:
  FieldElement( const Field<type,T>& field ) :
    field_(field)
  {}

  void retrieve( index_t elem );

private:
  const Field<type,T>& field_;
};



template<class IntegrandCell>
class JacobianCell_impl : public CellIntegralType< JacobianCell_impl<IntegrandCell> >
{
public:
  typedef typename IntegrandCell::template ArrayQ<real_t> ArrayQ;
  typedef typename IntegrandCell::template MatrixQ<real_t> MatrixQ;

  typedef numpack::DLA::MatrixD<MatrixQ> MatrixElemClass;

  JacobianCell_impl( const IntegrandCell& fcn , numpack::MatrixScatterAdd<MatrixQ>& matrix ) :
    fcn_(fcn),
    matrix_(matrix)
  {}

  template<typename type>
  void integrate( const Field<type,ArrayQ>& field , index_t k )
  {

    FieldElement<type,ArrayQ> qelem( field );
    FieldElement<type,real_t> xelem( field.topology() );

    const index_t nb_dof = field.element().nb_dof();
    std::vector<index_t> global( nb_dof );

    MatrixElemClass matrix_pde_elem( nb_dof , nb_dof );

    auto integrand = fcn_.integrand( field );

    for (index_t elem = 0; elem < field.topology().nb(); elem++)
    {

      qelem.retireve( elem );
      xelem.retrieve( elem );

      integral( integrand , xelem , matrix_pde_elem );

      scatter_add( elem , global , nb_dof , matrix_pde_elem , matrix_ );
    }
  }

protected:
  template<class MatrixQ, template<class> class SparseMatrixType>
  void
  scatter_add( const index_t elem , std::vector<index_t>& global , const index_t nb_dof , numpack::DLA::MatrixD<MatrixQ>& matrix_pde_elem , SparseMatrixType<MatrixQ>& matrix )
  {
    matrix.scatter_add( matrix_pde_elem , nb_dof , global.data() , global.size() );
  }

private:
  const IntegrandCell& fcn_;
  numpack::MatrixScatterAdd<MatrixQ>& matrix_;

};

}

#endif
