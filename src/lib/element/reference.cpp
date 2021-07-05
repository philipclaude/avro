//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"
#include "common/tools.h"

#include "element/reference.h"
#include "element/simplex.h"

#include "numerics/functions.h"
#include "numerics/geometry.h"
#include "numerics/mat.h"

#include <cmath>

namespace avro
{

/*
const real_t*
ReferenceElement<Simplex>::get_reference_coordinate( index_t k ) const
{
  return &xref_[k*(number_+1)];
}

const index_t*
ReferenceElement<Simplex>::get_lattice_coordinate( index_t k ) const
{
  return &lref_[k*(number_+1)];
}
int
ReferenceElement<Simplex>::find_index( const real_t* x ) const
{
  std::vector<index_t> y(number_);
  for (coord_t d=0;d<number_;d++)
    y[d] = x[d]*order_;
  avro_implement;
}
*/

static void
next_index( const int n , const int q , bool& more , std::vector<index_t>& x )
{
  int i,j;

  avro_assert( x.size()==index_t(n) );

  if (!more)
  {
    if (n<1)
      avro_assert(false);


    more = true;
    j = 1;

    x[0] = q;
    for (i=1;i<n;i++)
      x[i] = 0;

    if (n==1)
      more = false;
  }
  else
  {
    j = n -1;
    for (i=n-2;0<=i;i--)
    {
      if (0 < x[i])
      {
        j = i;
        break;
      }
    }

    x[j] = x[j] -1;
    x[j+1] = q;

    for (i=0;i<=j;i++)
      x[j+1] = x[j+1] -x[i];

    for (i=j+2;i<n;i++)
      x[i] = 0;

    if (x[n-1]==index_t(q))
      more = false;
  }
}

/*
void
ReferenceElement<Simplex>::precalculate()
{
  // set the unit (equilateral) coordinates
  if (number_==0)
  {
    xunit_.push_back(1.); // not sure?
  }
  else if (number_==1)
  {
    xunit_.push_back(1.);
    xunit_.push_back(0.);
  }
  else if (number_==2)
  {
    xunit_.push_back(1.); xunit_.push_back(0.);
    xunit_.push_back(-.5); xunit_.push_back(  std::sqrt(.25*3.) );
    xunit_.push_back(-.5); xunit_.push_back( -std::sqrt(.25*3.) );
  }
  else if (number_==3)
  {
    xunit_.push_back(1.); xunit_.push_back(0.); xunit_.push_back(0.);
    xunit_.push_back(-1./3); xunit_.push_back(  std::sqrt(8./9.) ); xunit_.push_back(  0.);
    xunit_.push_back(-1./3); xunit_.push_back( -std::sqrt(2./9.) ); xunit_.push_back(  std::sqrt(2./3.));
    xunit_.push_back(-1./3); xunit_.push_back( -std::sqrt(2./9.) ); xunit_.push_back( -std::sqrt(2./3.));
  }
  else if (number_==4)
  {
    // this is not actually unit...lengths are sqrt(5/2)
    xunit_.push_back(1.); xunit_.push_back(0.); xunit_.push_back(0.); xunit_.push_back(0.);
    xunit_.push_back(-.25); xunit_.push_back(  std::sqrt(15./16.) ); xunit_.push_back(0.); xunit_.push_back(0.);
    xunit_.push_back(-.25); xunit_.push_back( -std::sqrt( 5./48.) ); xunit_.push_back(  std::sqrt(5./6. ) ); xunit_.push_back(0.);
    xunit_.push_back(-.25); xunit_.push_back( -std::sqrt( 5./48.) ); xunit_.push_back( -std::sqrt(5./24.) ); xunit_.push_back(  std::sqrt(5./8.) );
    xunit_.push_back(-.25); xunit_.push_back( -std::sqrt( 5./48.) ); xunit_.push_back( -std::sqrt(5./24.) ); xunit_.push_back( -std::sqrt(5./8.) );

    for (index_t k=0;k<xunit_.size();k++)
      xunit_[k] *= sqrt(2./5.);
  }
  else
  {
    //printf("warning: this is untested\n");
    matd<real_t> X(number_,number_+1);
    X.zero();

    for (coord_t i=0;i<number_;i++)
      X(i,i) = 1.;
    real_t a = (1. -std::sqrt(1.+number_))/number_;
    for (coord_t i=0;i<number_;i++)
      X(i,number_) = a;

    std::vector<real_t> c(number_,0.);
    for (coord_t i=0;i<number_;i++)
    for (coord_t j=0;j<number_+1;j++)
      c[i] += X(i,j)/(number_+1);

    for (coord_t j=0;j<number_+1;j++)
    for (coord_t i=0;i<number_;i++)
      X(i,j) = X(i,j) -c[i];

    real_t s = 0.;
    for (coord_t i=0;i<number_;i++)
      s += X(i,0)*X(i,0);
    s = std::sqrt(s);
    for (coord_t i=0;i<number_;i++)
    for (coord_t j=0;j<number_+1;j++)
      X(i,j) = X(i,j)/s;

    for (coord_t i=0;i<number_+1;i++)
    for (coord_t j=0;j<number_;j++)
      xunit_.push_back( X(j,i) );
  }

  vunit_ = sqrt(number_+1)/(numerics::factorial(number_)*sqrt(pow(2.,number_)));
  vorth_ = 1./numerics::factorial(number_);

  // set the orthogonal coordinates
  xorth_.resize( number_*(number_+1) );
  std::fill( xorth_.begin() , xorth_.end() , 0. );
  for (index_t k=0;k<number_;k++)
    xorth_[ (k+1)*number_ +k ] = 1.;

  real_t length = order_;
  if (order_==0)
  {
    // it doesn't really make sense to have a zero-order lagrange element
    // but we can hack it by setting the coordinates to the centroid
    for (index_t j=0;j<index_t(number_+1);j++)
    {
      if (j<number_)
      {
        lref_.push_back( 1 );
        xref_.push_back( 1./(number_+1) );
      }
      else
      {
        // added on june 18th
        lref_.push_back( 0 );
        xref_.push_back( 0.);
      }
    }
    return;
  }
  avro_assert( length>0 );

  bool more = false;
  std::vector<index_t> xp(number_+1);
  for (index_t i=0;;i++)
  {
    next_index( number_+1 , order_ , more , xp );

    for (index_t j=0;j<xp.size();j++)
    {
      lref_.push_back(xp[j]);
      xref_.push_back( real_t(xp[j])/real_t(length));
    }

    if (!more) break;
  }

  for (index_t k=0;k<nb_basis();k++)
  {
    bool anyzero = false;
    for (coord_t j=0;j<number_+1;j++)
    {
      if (lref_[k*(number_+1)+j]==0)
      {
        anyzero = true;
        break;
      }
    }
    if (!anyzero)
    {
      //printf("lref = (%lu,%lu,%lu)\n",lref_[k*(number_+1)],lref_[k*(number_+1)+1],lref_[k*(number_+1)+2]);
      interior_.push_back(k);
    }
  }
}
*/

template<typename type,int N, int P>
struct GetLagrangeNodes {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,1,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,2,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,3,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type,int P>
struct GetLagrangeNodes<type,4,P> {
  static void get( std::vector<real_t>& nodes );
};

template<typename type, int P>
void
GetLagrangeNodes<type,1,P>::get( std::vector<real_t>& nodes ) {

  const std::vector<real_t>& s = LagrangeNodes<type,1,P>::coord_s_;
  for (index_t k = 0; k < s.size(); k++) {
    nodes.push_back( s[k] );
  }
}

template<typename type, int P>
void
GetLagrangeNodes<type,2,P>::get( std::vector<real_t>& nodes ) {

  const std::vector<real_t>& s = LagrangeNodes<type,2,P>::coord_s_;
  const std::vector<real_t>& t = LagrangeNodes<type,2,P>::coord_t_;

  avro_assert( s.size() == t.size() );

  for (index_t k = 0; k < s.size(); k++) {
    nodes.push_back( s[k] );
    nodes.push_back( t[k] );
  }
}

template<typename type, int P>
void
GetLagrangeNodes<type,3,P>::get( std::vector<real_t>& nodes ) {

  const std::vector<real_t>& s = LagrangeNodes<type,3,P>::coord_s_;
  const std::vector<real_t>& t = LagrangeNodes<type,3,P>::coord_t_;
  const std::vector<real_t>& u = LagrangeNodes<type,3,P>::coord_u_;

  avro_assert( s.size() == t.size() );
  avro_assert( u.size() == s.size() );

  for (index_t k = 0; k < s.size(); k++) {
    nodes.push_back( s[k] );
    nodes.push_back( t[k] );
    nodes.push_back( u[k] );
  }
}

template<typename type, int P>
void
GetLagrangeNodes<type,4,P>::get( std::vector<real_t>& nodes ) {

  const std::vector<real_t>& s = LagrangeNodes<type,4,P>::coord_s_;
  const std::vector<real_t>& t = LagrangeNodes<type,4,P>::coord_t_;
  const std::vector<real_t>& u = LagrangeNodes<type,4,P>::coord_u_;
  const std::vector<real_t>& v = LagrangeNodes<type,4,P>::coord_v_;

  avro_assert( s.size() == t.size() );
  avro_assert( u.size() == s.size() );
  avro_assert( v.size() == s.size() );

  for (index_t k = 0; k < s.size(); k++) {
    nodes.push_back( s[k] );
    nodes.push_back( t[k] );
    nodes.push_back( u[k] );
    nodes.push_back( v[k] );
  }
}

template<>
index_t
ReferenceElement<Simplex>::nb_basis() const {
  return nb_simplex_basis(number_,order_);
}

template<typename type>
ReferenceElement<type>::ReferenceElement( coord_t number , coord_t order ) :
  number_(number),
  order_(order) {

  if (number_ == 1) {
    static const int N = 1;
    if      (order_ == 1) GetLagrangeNodes<type,N,1>::get(nodes_);
    else if (order_ == 2) GetLagrangeNodes<type,N,2>::get(nodes_);
    else if (order_ == 3) GetLagrangeNodes<type,N,3>::get(nodes_);
    else if (order_ == 4) GetLagrangeNodes<type,N,4>::get(nodes_);
    else avro_implement;
  }
  if (number_ == 2) {
    static const int N = 2;
    if      (order_ == 1) GetLagrangeNodes<type,N,1>::get(nodes_);
    else if (order_ == 2) GetLagrangeNodes<type,N,2>::get(nodes_);
    else if (order_ == 3) GetLagrangeNodes<type,N,3>::get(nodes_);
    else if (order_ == 4) GetLagrangeNodes<type,N,4>::get(nodes_);
    else avro_implement;
  }
  if (number_ == 3) {
    static const int N = 3;
    if      (order_ == 1) GetLagrangeNodes<type,N,1>::get(nodes_);
    else if (order_ == 2) GetLagrangeNodes<type,N,2>::get(nodes_);
    else if (order_ == 3) GetLagrangeNodes<type,N,3>::get(nodes_);
    else if (order_ == 4) GetLagrangeNodes<type,N,4>::get(nodes_);
    else avro_implement;
  }
  if (number_ == 4) {
    static const int N = 4;
    if      (order_ == 1) GetLagrangeNodes<type,N,1>::get(nodes_);
    else if (order_ == 2) GetLagrangeNodes<type,N,2>::get(nodes_);
    else if (order_ == 3) GetLagrangeNodes<type,N,3>::get(nodes_);
    else if (order_ == 4) GetLagrangeNodes<type,N,4>::get(nodes_);
    else avro_implement;
  }

  real_t length = order_;

  bool more = false;
  std::vector<index_t> xp(number_+1);
  std::vector<real_t> xref_;

  for (index_t i=0;;i++)
  {
    next_index( number_+1 , order_ , more , xp );

    for (index_t j=0;j<xp.size();j++)
    {
      lattice_.push_back(xp[j]);
      xref_.push_back( real_t(xp[j])/real_t(length));
    }

    if (!more) break;
  }

  std::vector<index_t> lattice2canonical( nb_basis() );
  for (index_t k = 0; k < nb_basis(); k++) {
    int idx = find_index( &xref_[(number_+1)*k] );
    lattice2canonical[k] = idx;
  }

  for (index_t k=0;k<nb_basis();k++)
  {
    bool anyzero = false;
    for (coord_t j=0;j<number_+1;j++)
    {
      if (lattice_[k*(number_+1)+j]==0)
      {
        anyzero = true;
        break;
      }
    }
    if (!anyzero)
    {
      //printf("lref = (%lu,%lu,%lu)\n",lref_[k*(number_+1)],lref_[k*(number_+1)+1],lref_[k*(number_+1)+2]);
      interior_.push_back(k);
    }
  }
}

template<>
int
ReferenceElement<Simplex>::find_index( const index_t* x ) const
{
  for (index_t i=0;i<nb_basis();i++)
  {
    const index_t* coord0 = get_lattice_coordinate(i);
    index_t distance = 0;
    for (coord_t d=0;d<number_;d++)
    {
      if (coord0[d]-x[d]!=0)
      {
        distance = 1;
        break;
      }
    }
    if (distance==0)
    {
      return i;
    }
  }
  return -1;
}

template<>
int
ReferenceElement<Simplex>::find_index( const real_t* x ) const {
  for (index_t i=0;i<nb_basis();i++)
  {
    real_t d = numerics::distance( &nodes_[(number_)*i] , x , number_ );
    if (d < 1e-12) return i;
  }
  avro_assert_not_reached;
  return -1;
}

template<typename type>
const real_t*
ReferenceElement<type>::get_reference_coordinate( const index_t k ) const {
  avro_assert( k < nodes_.size()/number_ );
  return &nodes_[k*number_];
}

template<typename type>
const index_t*
ReferenceElement<type>::get_lattice_coordinate( const index_t k ) const {
  avro_assert( k < lattice_.size()/(number_+1) );
  return &lattice_[k*(number_+1)];
}

template class ReferenceElement<Simplex>;

} // avro
