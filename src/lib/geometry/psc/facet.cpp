//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/psc/facet.h"
#include "geometry/psc/node.h"

#include "numerics/geometry.h"
#include "numerics/linear_algebra.h"

namespace avro
{

namespace PSC
{

Facet::Facet( Body* body , std::vector<std::shared_ptr<Entity>>& facets ) :
  Object(facets.size()/2,body->dim()),
  V_(dim_,number_),
  B_(number_,number_),
  x0_(dim_),
  dimension_(-1)
{
  avro_assert( number_>0 );
  for (index_t k=0;k<facets.size();k++)
    add_child(facets[k]);
  build_basis();
  body_ = body;
}

void
Facet::build_basis()
{
  // choose number+1 points to compute the basis
  std::vector<Entity*> children;
  get_children(children);

  // retrieve the nodes
  std::vector<Node*> nodes;
  for (index_t k=0;k<children.size();k++)
  {
    if (children[k]->number()==0)
      nodes.push_back( static_cast<Node*>(children[k]) );
  }
  // remove the following line (but keep the one after) when non-cube geometries are implemented
  avro_assert_msg( nodes.size() == pow(2,number_) , "|nodes| = %lu" , nodes.size() );
  avro_assert( nodes.size() > number_ );

  // store a set of vectors to compute a basis from
  matd<real_t> u(dim_,nodes.size()-1);
  for (index_t k=1;k<nodes.size();k++)
  for (index_t j=0;j<dim_;j++)
    u(j,k-1) = (*nodes[k])(j) - (*nodes[0])(j);

  for (index_t j=0;j<dim_;j++)
    x0_[j] = (*nodes[0])(j);

  // compute an orthornormal basis
  int result = numerics::range(u,basis_);
  avro_assert( result>=0 );
  avro_assert( basis_.n() == number_ );

  for (coord_t k=0;k<number_;k++)
  for (coord_t j=0;j<dim_;j++)
    V_(j,k) = basis_(j,k);

  matd<real_t> VtV = numerics::transpose(V_)*V_;
  for (coord_t i=0;i<number_;i++)
  for (coord_t j=0;j<number_;j++)
    B_(i,j) = VtV(i,j);

  real_t d = numerics::det(B_);
  if (d == 0.)
  {
    printf("bad basis!\n");
    printf("dim = %u\n",dim_);
    basis_.dump();
    V_.dump();
    B_.dump();
    avro_assert_not_reached;
  }

  // determine the dimension if any
  for (index_t k = 0; k < V_.m(); k++) {
    real_t sum = 0.0;
    for (index_t j = 0; j < V_.n(); j++)
      sum += V_(k,j);
    if (fabs(sum) < 1e-12) {
      dimension_ = k;
      break;
    }
  }
}

void
Facet::set_basis_by_name() {

  avro_assert( number_ == 3 );

  V_.zero();
  B_.zero();

  if (name_ == "tmin" || name_ == "tmax") {
    V_(0,0) = 1;
    V_(1,1) = 1;
    V_(2,2) = 1;
  }
  else if (name_ == "xmin" || name_ == "xmax") {
    V_(1,0) = 1;
    V_(2,1) = 1;
    V_(3,2) = 1;
  }
  else if (name_ == "ymin" || name_ == "ymax") {
    V_(2,0) = 1;
    V_(3,1) = 1;
    V_(0,2) = 1;
  }
  else if (name_ == "zmin" || name_ == "zmax") {
    V_(3,0) = 1;
    V_(0,1) = 1;
    V_(1,2) = 1;
  }

  matd<real_t> VtV = numerics::transpose(V_)*V_;
  for (coord_t i=0;i<number_;i++)
  for (coord_t j=0;j<number_;j++)
    B_(i,j) = VtV(i,j);
}

void
Facet::evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const
{
  avro_assert( x.size() == dim_ );
  avro_assert( basis_.m() == dim_ );
  avro_assert_msg( basis_.n() == number_ , "basis = %lu x %lu , number_ = %u, dim_ = %u" , basis_.m(),basis_.n(),number_,dim_ );

  // set the last barycentric coordinate
  vecd<real_t> u0(number_);
  for (coord_t d=0;d<number_;d++)
  {
    u0(d) = u[d];
  }

  // compute the linear combination of the basis using the barycentric coordinates
  for (coord_t j=0;j<dim_;j++)
  {
    x[j] = x0_(j);
    for (coord_t k=0;k<number_;k++)
      x[j] += V_(j,k)*u0[k];
  }
}

void
Facet::inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  // initialize parameter values to something useless
  std::fill( u.begin() , u.end() , 1e20 );

  // compute right-hand-side of minimization statement
  vecd<real_t> x0(x.size(),x.data());
  vecd<real_t> b = numerics::transpose(V_)*( x0 - x0_ );

  // solve VtV * alpha = Vt * (x - x0)
  vecd<real_t> alpha(number_);
  numerics::solveLUP(B_,b,alpha);

  // compute the closest point y = V*alpha + x0
  vecd<real_t> y = V_*alpha + x0_;

  real_t d = numerics::distance2( y.data() , x.data() , dim_ );
  if (d<1e-12)
  {
    // save the parameter values
    for (int j=0;j<number_;j++)
      u[j] = alpha[j];
    return;
  }

  // the projected point is outside the polytope
  // set it to something far away
  for (coord_t d=0;d<dim_;d++)
    x[d] = 1e22; // the 1e22 is a code to help look up possible errors
}

} // PSC

} // avro
