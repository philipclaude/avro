#include "geometry/psc/object.h"

#include "numerics/geometry.h"
#include "numerics/linear_algebra.h"

namespace avro
{

namespace PSC
{

Object::Object( coord_t number , coord_t dim ) :
  Entity(number),
  dim_(dim)
  //basis_(dim,number)
{
  tessellatable_ = true;
  interior_ = false;
}

void
Facet::inverse( std::vector<real_t>& x , std::vector<real_t>& u ) const
{
  numerics::VectorD<real_t> x0(x.size(),x.data());
  numerics::VectorD<real_t> b0 = numpack::Transpose(V_)*x0;
  numerics::VectorD<real_t> b( number_+2 );
  for (index_t j=0;j<x.size();j++)
    b[j] = b0[j];
  b[number_+1] = 1.0;

  // solve B*y = b
  numerics::VectorD<real_t> y = numpack::DLA::InverseLU::Solve(B_,b);

  // compute the closest point (z)
  numerics::VectorD<real_t> alpha(number_+1);
  for (index_t j=0;j<number_+1;j++)
    alpha[j] = y[j];
  numerics::VectorD<real_t> z = V_*alpha;

  for (index_t j=0;j<alpha.size();j++)
    u[j] = alpha[j];

  real_t d = numerics::distance2( z.data() , x.data() , dim_ );
  //printf("d = %g\n",d);
  if (d<1e-12) return;

  // the projected point is outside the polytope
  // set it to something far away
  for (coord_t d=0;d<dim_;d++)
    x[d] = 1e20;

}

void
Node::inverse( std::vector<real_t>& x , std::vector<real_t>& ) const
{
  avro_assert( x.size() == dim_ );
  for (index_t j=0;j<x.size();j++)
    x[j] = (*this)(j);
}

void
Facet::evaluate( const std::vector<real_t>& u , std::vector<real_t>& x ) const
{
  avro_assert( u.size() == dim_-1 );
  avro_assert( x.size() == dim_ );
  avro_assert( basis_.m() == dim_ );
  avro_assert_msg( basis_.n() == number_ , "basis = %d x %d , number_ = %u, dim_ = %u" , basis_.m(),basis_.n(),number_,dim_ );
  for (index_t j=0;j<dim_;j++)
  {
    x[j] = 0.0;
    for (index_t k=0;k<basis_.n();k++)
    {
      x[j] += u[k]*basis_(j,k);
    }
  }
}

void
Node::evaluate( const std::vector<real_t>& , std::vector<real_t>& x ) const
{
  avro_assert( x.size() == dim_ );
  for (index_t j=0;j<x.size();j++)
    x[j] = (*this)(j);
}

void
Object::build_hierarchy()
{
  // not necesssary but declared virtual...maybe should just
  // use in egads geometry
  avro_assert_not_reached;
}

void
Object::print( bool with_children ) const
{
  printf("PSC entity with number %hu at %p\n",number_,(void*)(this));
  if (!with_children) return;
  for (index_t k=0;k<nb_children();k++)
    child(k).print(with_children);
}

Facet::Facet( Body* body , std::vector<std::shared_ptr<Entity>>& facets ) :
  Object(facets.size()/2,body->dim()),
  V_(dim_,number_+1),
  B_(number_+2,number_+2)
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
  avro_assert( nodes.size() > number_ );

  // store a set of vectors to compute a basis from
  numerics::MatrixD<real_t> u(dim_,nodes.size()-1);
  for (index_t k=1;k<nodes.size();k++)
  for (index_t j=0;j<dim_;j++)
    u(j,k-1) = (*nodes[k])(j) - (*nodes[0])(j);

  // compute an orthornormal basis
  int result = numerics::range(u,basis_);
  avro_assert( result>=0 );
  avro_assert( basis_.n() == number_ );

  for (index_t j=0;j<dim_;j++)
    V_(j,0) = (*nodes[0])(j);
  for (index_t k=1;k<number_+1;k++)
  for (index_t j=0;j<dim_;j++)
    V_(j,k) = (*nodes[0])(j) + basis_(j,k-1);

  B_ = 1;
  numerics::MatrixD<real_t> VtV = numpack::Transpose(V_)*V_;
  for (index_t i=0;i<number_+1;i++)
  for (index_t j=0;j<number_+1;j++)
    B_(i,j) = VtV(i,j);
}

} // PSC

} // avro
