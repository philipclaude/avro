#ifndef AVRO_LIB_GRAPHICS_BSP_H_
#define AVRO_LIB_GRAPHICS_BSP_H_

#include "common/error.h"

#include "element/simplex.h"

#include "graphics/math.h"
#include "graphics/plot.h"
#include "graphics/primitives.h"
#include "graphics/vao.h"

#include "mesh/topology.h"

#include "numerics/vec.hpp"

#include <memory>
#include <vector>

namespace avro
{

namespace graphics
{

#define BSP_TOL 1e-7

typedef struct {

  vec3 normal;
  vec3 center;

  int side( const vec3& p ) const {
    float dp = glm::dot( p - center , normal );
    if (fabs(dp) < BSP_TOL) return 0;
    if (dp > 0) return 1;
    return -1;
  }

  vec3 intersect( const vec3& p0 , const vec3& p1 ) const {
    if ( glm::norm(p1-p0) < BSP_TOL ) return p0;
    float num = glm::dot( p0 - center , normal );
    float den = glm::dot( normal , p0 - p1 );
    float d = num/den;
    if (d < -BSP_TOL || d > 1.0+BSP_TOL) {
      /*
      p0.print();
      p1.print();
      normal.print();
      center.print();
      printf("numerator = %g, denominator = %g\n",num,den);
      */
      return p0;
    }
    //avro_assert_msg( d >= -BSP_TOL && d <= (1.0+BSP_TOL) , "d = %g" , d );
    return p0 + d*(p1 - p0);
  }
} BSPPlane;

class BSPTriangle {

public:
  BSPTriangle( const vec3& p0 , const vec3& p1 , const vec3& p2 );
  const vec3& point( index_t j ) const { return points_[j]; }
  const BSPPlane& plane() const { return plane_; }
  bool ignore() const { return ignore_; }

private:
  vec3 points_[3];
  BSPPlane plane_;
  bool ignore_;
};

class BSPTriangles {

public:

  void build( const Plot& plot , const mat4& view_matrix , const mat4& projection_matrix , const mat4& screen_matrix );
  void add( std::shared_ptr<BSPTriangle> triangle ) {
    triangles_.push_back(triangle);
  }

  std::shared_ptr<BSPTriangle> operator[] (index_t k)       { return triangles_[k]; }
  std::shared_ptr<BSPTriangle> operator[] (index_t k) const { return triangles_[k]; }
  bool classify( index_t k , const BSPPlane& plane , BSPTriangles& front , BSPTriangles& back  ) const;
  index_t nb() const { return triangles_.size(); }

private:
  std::vector< std::shared_ptr<BSPTriangle> > triangles_;
};

class BSPTree {

public:
  BSPTree( index_t level = 0 ) : level_(level) {}

  void build( BSPTriangles& triangles );
  void get_triangles( const vec3& e , std::vector<BSPTriangle*>& triangles ) const;

private:
  BSPTriangles triangles_;
  std::shared_ptr<BSPTree> front_;
  std::shared_ptr<BSPTree> back_;
  index_t level_;
};

} // graphics

} // avro

#endif
