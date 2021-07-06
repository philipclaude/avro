//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_MASTER_SIMPLEX_H_
#define avro_LIB_MASTER_SIMPLEX_H_

#include "common/error.h"
#include "common/table.h"

#include "element/basis.h"
#include "element/element.h"
#include "element/reference.h"
#include "element/quadrature.h"

#include "numerics/functions.h"
#include "numerics/mat.h"

#include <map>
#include <vector>

namespace avro
{

template<typename type> class Topology;
template<typename type> class SimplicialDecomposition;
class Points;
class Entity;

index_t nb_simplex_basis( coord_t n , coord_t p );
index_t nb_simplex_basis_interior( coord_t n , coord_t p );

class Simplex : public Shape {
public:
  Simplex( const coord_t number , const coord_t order );
  Simplex( const Topology<Simplex>& topology , const coord_t order );

  static std::string type_name() { return "simplex"; }
  static TableLayoutCategory layout() { return TableLayout_Rectangular; }
  static index_t nb_basis( coord_t number , coord_t order ) { return nb_simplex_basis(number,order); }

  void precalculate();

  template<typename dof_t> void transfer( const Simplex&  , const std::vector<const dof_t*>& dof0 , std::vector<dof_t*>& dof1 , coord_t dim ) const;
  template<typename dof_t> void transfer( const Simplex&  , const std::vector<const dof_t>& dof0 , std::vector<dof_t>& dof1 ) const;

  template<typename dof_t>
  void convert( const Simplex& element_from , const std::vector<dof_t>& A , std::vector<dof_t>& B ) const;

  real_t volume( const Points& points , const index_t* v , index_t nv ) const;
  void edge_vector( const Points& points , index_t n0 , index_t n1 , real_t* edge , Entity* entity=nullptr ) const;

  index_t nb_facets( coord_t dim ) const {
    return numerics::nchoosek(number_+1,dim+1);
  }

  index_t nb_edges() const {
    return nb_facets(1);
  }

  index_t nb_points() const {
    return number_+1;
  }

  index_t nb_facets() const {
    return number_+1; // specialized for dim-1 facets
  }

  index_t nb_basis() const {
    return nb_basis(number_);
  }

  index_t nb_basis( coord_t dim ) const {
    index_t np = 1;
    for (coord_t d=1;d<=dim;d++)
      np *= (order_+d);
    return np/numerics::factorial(dim);
  }

  index_t nb_interior( coord_t dim ) const {
    index_t np = 1;
    for (coord_t d=1;d<=dim;d++)
      np *= (order_-d);
    return np/numerics::factorial(dim);
  }

  index_t nb_interior() const {
    return nb_interior(number_);
  }

  void get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , ElementIndices& f ) const;
  void get_edges( const index_t* v , index_t nv , std::vector<index_t>& ek ) const;
  void get_triangles( const index_t* v , index_t nv , std::vector<index_t>& tk ) const;

  real_t closest( const Points& x , const index_t* v , const index_t nv , const real_t* p , std::vector<real_t>& y , real_t distance2=-1 ) const;

  void get_canonical_indices( const index_t* v , index_t nv , const ElementIndices& f , std::vector<index_t>& canonical ) const;

  void facet( const index_t* v , index_t j , std::vector<index_t>& f ) const;

  index_t edge( index_t k , index_t i ) const;

  void triangulate( const index_t* v , index_t nv , SimplicialDecomposition<Simplex>& topology ) const;

  void set_entity( Entity* entity );
  Entity* entity() { return entity_; }

  void physical_to_reference( const Points& points , const index_t* v , index_t nv , const real_t* x , real_t* x0 ) const;

  real_t jacobian( const std::vector<const real_t*>& x , coord_t dim ) const;
  void   jacobian( const std::vector<const real_t*>& xk , matd<real_t>& J ) const;
  void   jacobian( const index_t* v , index_t nv , const Points& points , matd<real_t>& J ) const;

  const ReferenceElement<Simplex>& reference() const { return reference_; }
  void set_basis( BasisFunctionCategory category ) { reference_.set_basis(category); }

  void set_parameter( bool x ) { parameter_ = x; }
  bool parameter() const { return parameter_; }

  /*
  // REMOVE ME
  index_t nb_quad() const { return wquad_.size(); }
  void load_quadrature( Quadrature& quadrature ); // conical-product, grundmann-moeller, etc.
  real_t quad_weight(index_t k) const { return wquad_[k]; }
  const real_t* quad_point(index_t k) const { return &xquad_[number_*k]; }
  // END
  */

protected:
  void get_edge( const index_t* v , index_t nv , index_t iedge , index_t* e ) const;
  void get_triangle( const index_t* v , index_t nv , index_t itriangle , index_t* t ) const;
  index_t get_vertex( const index_t* v , index_t nv , index_t ivertex ) const;
  void get_facet_vertices( const index_t* v , index_t nv , index_t ifacet , std::vector<index_t>& f ) const;

private:

  // REMOVE ME
  //std::vector<real_t> xquad_;
  //std::vector<real_t> wquad_;
  // end

  ReferenceElement<Simplex> reference_;

  index_t triangle( const index_t itriangle , const index_t inode ) const;
  std::vector<index_t> facet( const index_t j , const index_t i ) const;

  // these are all straight-sided
  std::vector<index_t> vertices_; // also high-order
  std::vector<index_t> edges_;
  std::vector<index_t> triangles_;

  bool parameter_;
  Entity* entity_;
};

} // avro

#endif
