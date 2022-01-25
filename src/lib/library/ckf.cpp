//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "library/ckf.h"
#include "mesh/facets.h"
#include "numerics/geometry.h"

namespace avro
{

class Cube : public Table<index_t>
{

public:
  Cube( coord_t number ) :
    Table<index_t>(TableLayout_Rectangular,number+1),
    number_(number)
  {
    allocate();
    calculate();
  }

  ~Cube()
  {
    delete [] x_;
  }

  void triangulate();

  index_t nb_points() const
  {
    return pow(2,number_);
  }

  real_t* get_coordinate( const index_t i )
  {
    avro_assert( i<nb_points() );
    return x_+i*number_;
  }

  void shift_points()
  {
    // shift from [-1,1] to [0,1]
    for (index_t k=0;k<nb_points();k++)
    {
      real_t* y = get_coordinate(k);
      for (coord_t d=0;d<number_;d++)
      {
        y[d] = .5*(y[d] +1);
      }
    }
  }

  void shift_to( const real_t* x0 , const real_t* dx , std::vector<real_t>& y )
  {
    y.resize( number_*nb_points() );
    for (index_t i=0;i<nb_points();i++)
    {
      const real_t* v = get_coordinate(i);
      for (coord_t d=0;d<number_;d++)
      {
        y[ number_*i + d ] = x0[d] +v[d]*dx[d];
      }
    }
  }

private:

  void allocate()
  {
    x_ = new real_t[ number_*nb_points() ];
  }

  void calculate();

  index_t find_point( const real_t* v )
  {
    for (index_t i=0;i<nb_points();i++)
    {
      real_t* y = get_coordinate(i);
      real_t  d = 0.;
      for (coord_t j=0;j<number_;j++)
        d += ( y[j] -v[j] )*( y[j] -v[j] );

      if (d<1e-6) return i;
    }
    printf("could not find vertex (%g,%g,%g,%g)\n",v[0],v[1],v[2],v[3]);
    avro_assert(false);
    return nb_points();
  }

  index_t number_;
  real_t* x_;

};

void
Cube::calculate()
{
  avro_assert(number_>0);

  // base case
  if (number_==1)
  {
    x_[0] = -1;
    x_[1] = 1;
    return;
  }

  // get the previous dimension cube
  Cube y(number_-1);
  index_t k=0;
  for (index_t i=0;i<y.nb_points();i++)
  {
    real_t* v1 = get_coordinate(k);
    real_t* v2 = get_coordinate(k+1);
    real_t* yi = y.get_coordinate(i);
    for (coord_t d=0;d<number_-1;d++)
    {
      v1[d] = yi[d];
      v2[d] = yi[d];
    }

    v1[number_-1] = -1;
    v2[number_-1] =  1;

    k = k +2;
  }
}

void
Cube::triangulate()
{
  // perform the Kuhn-Friedenthal triangulation such that we can
  // tile the generated simplices to create a mesh
  shift_points();

  std::vector<int>    p(number_);
  std::vector<real_t> v(number_,0);
  for (index_t i=0;i<number_;i++)
    p[i] = i;

  index_t idx0 = find_point(v.data());

  std::vector<index_t> s(number_+1);
  do
  {
    s[0] = idx0;

    std::fill( v.begin() , v.end() , 0 );

    for (index_t j=0;j<number_;j++)
    {
      std::vector<real_t> vj(number_);
      for (index_t d=0;d<number_;d++)
        vj[d] = v[d];

      vj[p[j]] = vj[p[j]] +1;
      s[j+1]   = find_point(vj.data());

      for (index_t d=0;d<number_;d++)
        v[d] = vj[d];
    }

    add( s.data() , s.size() );
  } while (std::next_permutation(p.begin(),p.end()));
}

CKF_Triangulation::CKF_Triangulation( const std::vector<index_t>& dims ) :
  Topology<Simplex>(points_,dims.size()),
  points_(dims.size()),
  p_(number_,0),
  u_(number_,0),
  dims_(dims),
  grid_(TableLayout_Rectangular,number_)
{
  for (index_t d=0;d<dims.size();d++)
    avro_assert_msg( dims[d] > 1 , "need more than 1 point in dimension %lu!" , d );
  generate();
}

void
CKF_Triangulation::ndgrid( const std::vector<real_t>& dx , coord_t current )
{

  for (index_t i=0;i<dims_[current];i++)
  {
    p_[current] = dx[current]*real_t(i);
    u_[current] = i;
    if (current==number_-1)
    {
      // add the vertex
      points_.create(p_.data());
      grid_.add( u_.data() , u_.size() );
    }
    else
    {
      // call the next dimension
      ndgrid(dx,current+1);
    }
  }
}

bool
check_skip( const real_t* x , coord_t dim )
{
  for (coord_t d=0;d<dim;d++)
  {
    if ( fabs(x[d]-1.0)<1e-6 )
      return true;
  }
  return false;
}

void
CKF_Triangulation::generate()
{
  avro_assert( dims_.size() == index_t(number_) );

  std::vector<real_t> dx(number_);
  for (coord_t d=0;d<number_;d++)
    dx[d] = 1./real_t(dims_[d] -1.);

  ndgrid( dx , 0 );

  std::vector< std::vector<index_t> > S0;

  // create the base cube triangulation
  // i.e. triangulation of a single n-cube with 2^n points
  Cube base(number_);
  base.triangulate();

  Table<index_t> table(TableLayout_Rectangular,number_);
  for (index_t k=0;k<base.nb_points();k++)
  {
    const real_t* x = base.get_coordinate(k);
    std::vector<index_t> idx(x,x+number_);
    table.add( idx.data() , idx.size() );
  }

  std::vector<index_t> simplex(number_+1);
  std::vector<index_t> idx( base.nb_points() );

  std::vector<index_t> offset(number_,1);
  for (coord_t d=0;d<number_;d++)
    for (coord_t i=d+1;i<number_;i++)
      offset[d] *= dims_[i];

  for (index_t i=0;i<points_.nb();i++)
  {

    real_t* v0 = points_[i];
    std::vector<index_t> u0 = grid_.get(i);

    // if the base point has any coordinate = 1, skip
    if (check_skip(v0,number_)) continue;

    // find the map between the base cube
    // indices and those resulting from the offset
    for (index_t j=0;j<base.nb_points();j++)
    {
      idx[j] = 0;
      for (coord_t d=0;d<number_;d++)
        idx[j] += (u0[d]+table(j,d))*offset[d];
    }

    for (index_t k=0;k<base.nb();k++)
    {
      // map the simplex indices
      for (index_t j=0;j<index_t(number_+1);j++)
        simplex[j] = idx[ base(k,j) ];

      // add the element to the topology
      this->add( simplex.data() , simplex.size() );
    }
  }

  // orient the topology so the volumes are positive
  this->orient();
}


template<>
CubeDomain<Simplex>::CubeDomain( coord_t number , coord_t dim , index_t n ) :
  Topology<Simplex>(points_,number),
  points_(dim)
{
  std::vector<index_t> dims(number,n);
  CKF_Triangulation ckf(dims);

  // copy the points and elements from the ckf triangulation
  avro_assert( dim >= number );
  avro_assert( ckf.points().dim() == number );
  std::vector<real_t> p(dim,0.0);
  for (index_t k = 0; k < ckf.points().nb(); k++)
  {
    // copy the points, any additional coordinates will be zero
    for (coord_t d = 0; d < ckf.points().dim(); d++)
      p[d] = ckf.points()[k][d];
    points_.create( p.data() );
  }
  TopologyBase::copy( ckf );

  // compute the neighbours
  neighbours().compute();
  inverse().build();
}


template<>
CubeDomain<Polytope>::CubeDomain( coord_t number , coord_t dim , index_t n0 ) :
  Topology<Polytope>(points_,number),
  points_(dim)
{
  // create a CKF triangulation with only 2 points in each direction
  std::vector<index_t> dims(number,2);
  CKF_Triangulation ckf(dims);

  // compute the facets of the mesh
  Facets facets(ckf);
  facets.compute();

  // determine which facets are on a common hyperplane
  std::vector<std::vector<int>> v2b( ckf.points().nb() );
  std::vector<index_t> f(number);
  std::vector<real_t*> x(number);
  std::vector<real_t> n(number);
  std::vector<real_t> c(number);
  for (index_t k=0;k<facets.nb();k++)
  {
    if (!facets.boundary(k)) continue;
    facets.retrieve(k,f);

    // retrieve the facet points
    for (index_t j=0;j<number;j++)
      x[j] = ckf.points()[ f[j] ];

    // compute the normal to the facet
    numerics::normal( x , n.data() , number );

    // which direction is this in?
    coord_t dir = 0;
    for (coord_t d=0;d<number;d++)
    {
      if (fabs( fabs(n[d])-1.0) < 1e-8)
      {
        dir = d;
        break;
      }
    }

    // are we at 0 or 1?
    numerics::centroid( f.data() , f.size() , ckf.points() , c );
    bool zero = true;
    if (fabs(c[dir]) > 1e-8)
    {
      avro_assert( fabs(c[dir] -1.0) < 1e-8 );
      zero = false;
    }

    // determine which plane this corresponds to
    int b;
    if (zero) b = - dir - 1;
    else b = - number - dir - 1;

    // let each vertex know it is on this bisector
    for (coord_t d=0;d<f.size();d++)
    {
      v2b[f[d]].push_back( b );
    }
  }

  // copy the points from the ckf triangulation
  avro_assert( dim >= number );
  avro_assert( ckf.points().dim() == number );
  std::vector<real_t> p(dim,0.0);
  for (index_t k = 0; k < ckf.points().nb(); k++)
  {
    // copy the points, any additional coordinates will be zero
    for (coord_t d = 0; d < ckf.points().dim(); d++)
      p[d] = ckf.points()[k][d];
    points_.create( p.data() );
  }

  // set the incidence relations
  for (index_t k=0;k<ckf.points().nb();k++)
  {
    uniquify( v2b[k] );
    points_.incidence().add( v2b[k].data() , v2b[k].size() );
  }

  std::vector<index_t> cube = linspace( std::pow(2,number) );
  add( cube.data() , cube.size() );
}

} // avro
