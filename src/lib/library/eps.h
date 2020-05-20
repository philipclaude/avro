#ifndef avro_LIB_LIBRARY_EPS_H_
#define avro_LIB_LIBRARY_EPS_H_

#include "common/types.h"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace avro
{

namespace library
{

class epsFile
{
public:
  epsFile();

  void add_triangles( const std::vector<real_t>& triangles , const std::vector<real_t>& colors );
  void add_edges( const std::vector<real_t>& edges , const std::vector<real_t>& colors );
  void add_points( const std::vector<real_t>& points );

  void write( const std::string& filename );

  void set_viewport( int* viewport );

private:
  FILE* fid_;
  void print_triangles() const;
  void print_edges() const;

  void print_primitives();

  std::vector<int> viewport_;
  std::vector<real_t> triangles_;
  std::vector<real_t> triangle_colors_;

  std::vector<real_t> edges_;
  std::vector<real_t> edge_colors_;

  std::vector<real_t> points_;

  class Primitive
  {
  public:
    Primitive( coord_t nv , coord_t nd , const real_t* x , const real_t* c , index_t id ) :
      nb_vertices_(nv),
      depth_(0),
      id_(id)
    {
      x_.resize( nd*nv , 0 ); // xy or xyz
      c_.resize( 3*nv , 0 );  // rgb
      for (coord_t i=0;i<nv;i++)
      {
        // save the nd-coordinates for the vertex
        for (coord_t j=0;j<nd;j++)
        {
          x_[i*nd+j] = 100*x[i*nd+j];
        }
        depth_ += x_[i*nd+2];

        // save the the (r,g,b) triplet for the vertex
        for (coord_t j=0;j<3;j++)
          c_[i*nd+j] = c[i*nd+j];
      }

      depth_ /= (real_t) nv;
    }

    real_t depth() const { return depth_; }
    index_t id() const { return id_; }

    bool coplanar( const Primitive& x , const Primitive& y ) const
    {

      real_t n[3];
      real_t v[3] = {0.,0.,0.};
      real_t t[3];
      real_t x1,y1,z1;
      real_t x2,y2,z2;
      real_t x3,y3,z3;

      if (x.nb_vertices()==y.nb_vertices()) return false;

      if (x.nb_vertices()==3 && y.nb_vertices()==2)
      {
        x1 = x.vertices()[0];
        y1 = x.vertices()[1];
        z1 = x.vertices()[2];
        x2 = x.vertices()[3];
        y2 = x.vertices()[4];
        z2 = x.vertices()[5];
        x3 = x.vertices()[6];
        y3 = x.vertices()[7];
        z3 = x.vertices()[8];
        v[0] = y.vertices()[3] - y.vertices()[0];
        v[1] = y.vertices()[4] - y.vertices()[1];
        v[0] = y.vertices()[5] - y.vertices()[2];
        t[0] = y.vertices()[0];
        t[1] = y.vertices()[1];
        t[2] = y.vertices()[2];
      }
      else if (x.nb_vertices()==2 && y.nb_vertices()==3)
      {
        x1 = y.vertices()[0];
        y1 = y.vertices()[1];
        z1 = y.vertices()[2];
        x2 = y.vertices()[3];
        y2 = y.vertices()[4];
        z2 = y.vertices()[5];
        x3 = y.vertices()[6];
        y3 = y.vertices()[7];
        z3 = y.vertices()[8];
        v[0] = x.vertices()[3] - x.vertices()[0];
        v[1] = x.vertices()[4] - x.vertices()[1];
        v[0] = x.vertices()[5] - x.vertices()[2];
        t[0] = x.vertices()[0];
        t[1] = x.vertices()[1];
        t[2] = x.vertices()[2];
      }
      else
        return false;

      n[0] = (y3-y1)*(z2-z1)-(y2-y1)*(z3-z1);
      n[1] = (z3-z1)*(x2-x1)-(x3-x1)*(z2-z1);
      n[2] = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);

      // check if n is perpendicular to v
      real_t dp = n[0]*v[0] + n[1]*v[1] + n[2]*v[2];

      if (fabs(dp)<1e-6)
      {
        dp = (t[0] - x1)*n[0] + (t[1] - y1)*n[1] + (t[2] -z1)*n[2];
        if (fabs(dp)<1e-6)
        {
          printf("coplanar!!\n");
          return true;
        }
      }

      return false;
    }


    bool operator< ( const Primitive& y ) const
    {
      if (coplanar(*this,y))
      {
        // this might not resolve z-fighting, but it's worth a try
        if (nb_vertices_==y.nb_vertices()) return id_ < y.id();

        // paint lower-dimensional entities last (lower depth)
        return nb_vertices_ < y.nb_vertices();
      }
      return depth_ < y.depth();
    }

    const real_t* color() const { return c_.data(); }
    const real_t* vertices() const { return x_.data(); }

    void print() const
    {
      for (index_t k=0;k<nb_vertices_;k++)
      {
      }
    }

    index_t nb_vertices() const { return nb_vertices_; }

  private:
    index_t nb_vertices_;
    std::vector<real_t> x_;
    std::vector<real_t> c_;
    real_t depth_;
    index_t id_;
  };

  std::vector<Primitive> primitives_;
};

} // library

} // avro


#endif
