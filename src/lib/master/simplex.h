#ifndef LUNA_LIB_MASTER_SIMPLEX_H_
#define LUNA_LIB_MASTER_SIMPLEX_H_

#include "common/error.h"

#include "master/basis.h"
#include "master/master.h"

#include "numerics/functions.h"
#include "numerics/matrix.h"
#include "numerics/types.h"

#include <vector>

namespace luna
{

template<typename Shape> class Topology;
class Points;

class Simplex : public Master<Simplex>
{
public:
  Simplex( const coord_t number , const coord_t order ) :
    Master(number,order,"simplex")
  {}

  Simplex( const Topology<Simplex>& topology , const coord_t order );

  void precalculate();

  template<typename dof_t> void transfer( const Simplex& master , const std::vector<const dof_t*>& dof0 , std::vector<dof_t*>& dof1 , coord_t dim ) const;
  template<typename dof_t> void transfer( const Simplex& master , const std::vector<const dof_t>& dof0 , std::vector<dof_t>& dof1 ) const;

  template<typename dof_t>
  void convert( const Simplex& master_from , const std::vector<dof_t>& A , std::vector<dof_t>& B ) const;

  real_t volume( const Points& points , const index_t* v , index_t nv ) const;

  index_t nb_poly() const { return this->phi_.m(); }
  index_t nb_quad() const { return this->phi_.n(); }

  index_t nb_facets( coord_t dim ) const
  {
    return numerics::nchoosek(number_+1,dim+1);
  }

  index_t nb_edges() const
  {
    return nb_facets(1);
  }

  index_t nb_points() const
  {
    return number_+1;
  }

  index_t nb_facets() const
  {
    return number_+1; // specialized for dim-1 facets
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
    return np/numerics::factorial(dim);
  }

  index_t nb_interior( coord_t dim ) const
  {
    index_t np = 1;
    for (coord_t d=1;d<=dim;d++)
      np *= (order_-d);
    return np/numerics::factorial(dim);
  }

  index_t nb_interior() const
  {
    return nb_interior(number_);
  }

  index_t get_index( index_t dim , index_t ifacet , index_t ilocal ) const;
  void get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , Element& f ) const;
  void get_edges( const index_t* v , index_t nv , std::vector<index_t>& ek ) const;

  real_t closest( const Points& x , const index_t* v , const index_t nv , const real_t* p , std::vector<real_t>& y ) const;

  void facet( const index_t* v , index_t j , std::vector<index_t>& f ) const
    { luna_implement; }

  index_t edge( index_t k , index_t i ) const;

protected:
  void get_edge( const index_t* v , index_t nv , index_t iedge , index_t* e ) const;
  void get_triangle( const index_t* v , index_t nv , index_t itriangle , index_t* t ) const;
  index_t get_vertex( const index_t* v , index_t nv , index_t ivertex ) const;
  void get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , std::vector<index_t>& f ) const;

private:

  // these are all straight-sided
  std::vector<index_t> vertices_; // also high-order
  std::vector<index_t> edges_;
  std::vector<index_t> triangles_;

  // high-order edges
  std::vector<index_t> full_edges_;

  // store the transformation matrix from a lagrange simplex
  numerics::MatrixD<real_t> transformation_;

  real_t vunit_;
  real_t vorth_;

};

} // luna

#endif
