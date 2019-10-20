#ifndef URSA_LIB_MASTER_SIMPLEX_H_
#define URSA_LIB_MASTER_SIMPLEX_H_

#include "common/error.h"

#include "master/basis.h"
#include "master/master.h"

#include "numerics/matrix.h"
#include "numerics/types.h"

#include <vector>

namespace ursa
{

template<typename Shape_t> class Topology;

template<typename Basis> class Simplex;
class Quadrature;

inline index_t
factorial( index_t n )
{
  if (n==0) return 1;
  return n*factorial(n-1);
}

inline index_t
nchoosek( index_t n , index_t k )
{
  return factorial(n)/( factorial(k)*factorial(n-k) );
}

template<typename Basis>
class SimplexBase : public Master
{
public:
  SimplexBase( const coord_t number , const coord_t order ) :
    Master(number,order,"simplex")
  {}

  SimplexBase( const Topology<Simplex<Basis>>& topology , const coord_t order );

  void loadQuadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.

  template<typename BasisFrom_t,typename dof_t> void convert( const BasisFrom_t& masterFrom , const std::vector<dof_t>& A , std::vector<dof_t>& B ) const;

  index_t nb_poly() const { return phi_.m(); }
  index_t nb_quad() const { return phi_.n(); }

  index_t nb_facets( coord_t dim ) const
  {
    return nchoosek(number_+1,dim+1);
  }

  index_t nb_edges() const
  {
    return nb_facets(1);
  }

  index_t nb_vertices() const
  {
    return number_+1;
  }

  index_t nb_facets() const
  {
    // specialized for dim-1 facets
    return number_+1;
  }

  index_t nb_basis() const
  {
    return nb_basis(number_);
  }

  index_t nb_basis( coord_t dim ) const
  {
    index_t np = 1;
    for (coord_t d=1;d<=dim;d++)
      np *= (order_+d);
    return np/factorial(dim);
  }

  index_t nb_interior( coord_t dim ) const
  {
    index_t np = 1;
    for (coord_t d=1;d<=dim;d++)
      np *= (order_-d);
    return np/factorial(dim);
  }

  index_t nb_interior() const
  {
    return nb_interior(number_);
  }

  index_t get_index( index_t dim , index_t ifacet , index_t ilocal ) const;
  void get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , Element& f ) const;
  void get_edges( const index_t* v , index_t nv , std::vector<index_t>& edges ) const;
  void get_triangles( const index_t* v , index_t nv , std::vector<index_t>& triangles ) const;

protected:
  void get_edge( const index_t* v , index_t nv , index_t iedge , index_t* e ) const;
  void get_triangle( const index_t* v , index_t nv , index_t itriangle , index_t* t ) const;
  index_t get_vertex( const index_t* v , index_t nv , index_t ivertex ) const;
  void get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , std::vector<index_t>& f ) const;

private:

  numerics::MatrixD<real_t> phi_;
  std::vector< numerics::MatrixD<real_t> > dphi_;

  std::vector<real_t> xquad_;
  std::vector<real_t> wquad_;

protected:
  std::vector<index_t> edges_;
  std::vector<index_t> triangles_;
};

template<>
class Simplex<Lagrange> : public SimplexBase<Lagrange>
{
public:
  Simplex( coord_t number , coord_t order );

  template<typename BasisFrom_t,typename T> void transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const;
  template<typename BasisFrom_t,typename T> void transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T>& X , std::vector<T>& Y ) const;

  const real_t* get_reference_coordinate( index_t k ) const;
  const index_t* get_lattice_coordinate( index_t k ) const;
  void evaluate( const real_t* x , std::vector<real_t>& phi ) const;
  void evaluate( index_t k , std::vector<real_t>& phi ) const;

  void precalculate();

  void eval() const { printf("calling lagrange simplex eval\n"); }
  void eval() {}

private:

  std::vector<real_t> xunit_;
  std::vector<real_t> xorth_;

  std::vector<real_t>  xref_;
  std::vector<index_t> lref_;
};

template<>
class Simplex<Bezier> : public SimplexBase<Bezier>
{
public:
  Simplex( coord_t number , coord_t order );

  template<typename BasisFrom_t,typename T> void transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const;
  template<typename BasisFrom_t,typename T> void transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T>& X , std::vector<T>& Y ) const;

  void evaluate( const real_t* x , real_t* y ) const;
  void evaluate( index_t k , std::vector<real_t>& phi ) const;

  void eval() const { printf("calling bezier simplex eval\n"); }
  void eval() {}

private:
  // store the transformation matrix from a lagrange simplex
  numerics::MatrixD<real_t> B2L_;
};

template<>
class Simplex<Legendre> : public SimplexBase<Legendre>
{
public:
  Simplex( coord_t number , coord_t order );

  template<typename BasisFrom_t,typename T> void transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T*>& X , std::vector<T*>& Y , coord_t dim ) const;
  template<typename BasisFrom_t,typename T> void transfer( const Simplex<BasisFrom_t>& master , const std::vector<const T>& X , std::vector<T>& Y ) const;

  void evaluate( const real_t* x , real_t* y ) const;
  void evaluate( index_t k , std::vector<real_t>& phi ) const;

  void eval() const { printf("calling bezier simplex eval\n"); }
  void eval() {}

private:

};

} // ursa

#endif
