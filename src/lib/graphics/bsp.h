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

typedef struct {
  vec3 normal;
  vec3 center;

  int side( const vec3& p ) const {

    float dp = glm::dot( p - center , normal );
    if (fabs(dp) < 1e-12) return 0;
    if (dp > 0) return 1;
    return -1;
  }

  vec3 intersect( const vec3& p0 , const vec3& p1 ) const {
    float d = glm::dot( p0 - center , normal )/glm::dot(normal,p0-p1);
    return p0 + d*(p1-p0);
  }
} BSPPlane;

class BSPTriangle {
public:
  BSPTriangle( const vec3& p0 , const vec3& p1 , const vec3& p2 ) {
    points_[0] = p0;
    points_[1] = p1;
    points_[2] = p2;

    // initialize the plane of this triangle
    plane_.center = p0;
    plane_.normal = glm::cross( p1 - p0 , p2 - p0 );

    ignore_ = false;
    real_t l = glm::norm( plane_.normal );
    if (l < 1e-8) {
      ignore_ = true;
    }
    else {
      plane_.normal = glm::normalize(plane_.normal);
    }

    if (std::isnan(plane_.normal[0]) || std::isnan(plane_.normal[1]) || std::isnan(plane_.normal[2])) {
      ignore_ = true;
    }
  }

  const vec3& point( index_t j ) const { return points_[j]; }

  const BSPPlane& plane() const { return plane_; }

  bool ignore() const { return ignore_; }

private:
  vec3 points_[3];
  BSPPlane plane_;
  bool ignore_;
};

typedef struct {
  int side;
  bool swap;
  std::shared_ptr<BSPTriangle> triangle[3];
  int nb_triangles;
  int side0, side1, side2;
} BSPClassification;

class BSPTriangles {
public:

  BSPTriangles() {}

  void build( const Plot& plot , const mat4& view_matrix , const mat4& projection_matrix , const mat4& screen_matrix ) {

    const VertexAttributeObject& vao = plot.active_vao();
    const mat4& model_matrix = plot.model_matrix();

    // first map all the points
    std::vector<vec3> points;
    const std::vector<gl_float>& coordinates = vao.points().coordinates();
    mat4 transformation = /*screen_matrix * projection_matrix * */ view_matrix * model_matrix;
    for (index_t k = 0; k < coordinates.size()/3; k++) {

      vec4 p = { coordinates[3*k ] , coordinates[3*k+1] , coordinates[3*k+2] , 1.0f };
      vec4 q = transformation * p;

      vec3 v = { q[0] , q[1] , q[2] };
      points.push_back(v);
    }

    // loop through all triangle primitives
    for (index_t k = 0; k < vao.nb_triangles(); k++) {

      const TrianglePrimitive& prim = vao.triangles(k);

      const std::vector<gl_index>& indices = prim.indices();
      for (index_t j = 0; j < indices.size()/3; j++) {

        const vec3& p0 = points[ indices[3*j  ] ];
        const vec3& p1 = points[ indices[3*j+1] ];
        const vec3& p2 = points[ indices[3*j+2] ];

        triangles_.push_back( std::make_shared<BSPTriangle>(p0,p1,p2) );
      }
    }
  }

  void add( std::shared_ptr<BSPTriangle> triangle ) {
    triangles_.push_back(triangle);
  }

  std::shared_ptr<BSPTriangle> operator[] (index_t k) { return triangles_[k]; }

  bool classify( index_t k , const BSPPlane& plane , BSPTriangles& front , BSPTriangles& back  ) {

    // compute the side for all three vertices
    int s0 = plane.side( triangles_[k]->point(0) );
    int s1 = plane.side( triangles_[k]->point(1) );
    int s2 = plane.side( triangles_[k]->point(2) );

    if (s0 >= 0 && s1 >= 0 && s2 >= 0) {
      // triangle is entirely on the positive side
      front.add( triangles_[k] );
      return false;
    }
    if (s0 <= 0 && s1 <= 0 && s2 <= 0) {
      // triangle is entirely on the negative side
      back.add( triangles_[k] );
      return false;
    }

    // see if any vertices lie exactly on the plane
    int nb_coincident = 0;
    if (s0 == 0) nb_coincident++;
    if (s1 == 0) nb_coincident++;
    if (s2 == 0) nb_coincident++;

    // the entire triangle is coincident and should be added to the
    // node in the tree that called this function
    if (nb_coincident == 3) {
      return true;
    }

    #define ADD_ONE(s) if (s > 0) { front.add(triangles_[k]); } else { back.add(triangles_[k]); }
    #define ADD_TWO(s,t) if (s < 0) { avro_assert( t > 0 ); if (!t0->ignore()) back.add(t0); if (!t1->ignore()) front.add(t1); } else { if (!t1->ignore()) back.add(t1); if (!t0->ignore()) front.add(t0); }

    // one vertex is coincident to the plane
    if (nb_coincident == 1) {

      // determine if the other two are on opposite sides
      if (s0 == 0) {

        if (s1*s2 < 0) { // different signs for vertices 1 & 2

          // split the triangle in two
          vec3 q = plane.intersect( triangles_[k]->point(1) , triangles_[k]->point(2) );

          std::shared_ptr<BSPTriangle> t0 = std::make_shared<BSPTriangle>( triangles_[k]->point(0) , triangles_[k]->point(1) , q );
          std::shared_ptr<BSPTriangle> t1 = std::make_shared<BSPTriangle>( q , triangles_[k]->point(2) , triangles_[k]->point(0) );
          ADD_TWO(s1,s2);
        }
        else {
          ADD_ONE(s1); // vertices 1 & 2 are on the same side, so add to the side of vertex 1
        }
      }
      else if (s1 == 0) {
        if (s0*s2 < 0) { // different signs for vertices 0 & 2

          // split the triangle in two
          vec3 q = plane.intersect( triangles_[k]->point(0) , triangles_[k]->point(2) );

          std::shared_ptr<BSPTriangle> t0 = std::make_shared<BSPTriangle>( triangles_[k]->point(1) , triangles_[k]->point(2) , q );
          std::shared_ptr<BSPTriangle> t1 = std::make_shared<BSPTriangle>( q , triangles_[k]->point(0) , triangles_[k]->point(1) );
          ADD_TWO( s2 , s0 );
        }
        else {
          ADD_ONE( s0 ); // vertices 0 & 2 are on the same side, so add to the side of vertex 0
        }
      }
      else if (s2 == 0) {
        if (s0*s1 < 0) { // different signs for vertices 0 & 1

          // split the triangle in two
          vec3 q = plane.intersect( triangles_[k]->point(0) , triangles_[k]->point(1) );

          std::shared_ptr<BSPTriangle> t0 = std::make_shared<BSPTriangle>( triangles_[k]->point(2) , triangles_[k]->point(0) , q );
          std::shared_ptr<BSPTriangle> t1 = std::make_shared<BSPTriangle>( q , triangles_[k]->point(1) , triangles_[k]->point(2) );
          ADD_TWO(s0,s1);
        }
        else {
          ADD_ONE(s0); // vertice 0 & 1 are on the same side, so add to the side of vertex 0
        }
      }
      else
        avro_assert_not_reached;
    }
    else if (nb_coincident == 2) {
      // return the side of the one that is not coincident
      if (s0 != 0) {
        avro_assert( s1 == 0 && s2 == 0 );
        ADD_ONE( s0 );
      }
      else if (s1 != 0) {
        avro_assert( s0 == 0 && s2 == 0 );
        ADD_ONE( s1 );

      }
      else if (s2 != 0) {
        avro_assert( s0 == 0 && s1 == 0 );
        ADD_ONE( s2 );
      }
      else
        avro_assert_not_reached;
    }
    else {

      avro_assert( nb_coincident == 0 );

      int fa = s0;
      int fb = s1;
      int fc = s2;

      vec3 a = triangles_[k]->point(0);
      vec3 b = triangles_[k]->point(1);
      vec3 c = triangles_[k]->point(2);

      if (fa * fc >= 0) {
        std::swap( fb , fc );
        std::swap( b , c );
        std::swap( fa , fb );
        std::swap( a , b );
      }
      else if ( fb * fc >= 0 ) {
        std::swap( fa , fc );
        std::swap( a , c );
        std::swap( fa , fb );
        std::swap( a , b );
      }

      vec3 A = plane.intersect( a , b );
      vec3 B = plane.intersect( b , c );

      std::shared_ptr<BSPTriangle> T1 = std::make_shared<BSPTriangle>(a,b,A);
      std::shared_ptr<BSPTriangle> T2 = std::make_shared<BSPTriangle>(b,B,A);
      std::shared_ptr<BSPTriangle> T3 = std::make_shared<BSPTriangle>(A,B,c);

      if (fc >= 0) {
        if (!T1->ignore()) back.add(T1);
        if (!T2->ignore()) back.add(T2);
        if (!T3->ignore()) front.add(T3);
      }
      else {
        if (!T1->ignore()) front.add(T1);
        if (!T2->ignore()) front.add(T2);
        if (!T3->ignore()) back.add(T3);
      }
    }
    return false;
  }

  index_t nb() const { return triangles_.size(); }

private:

  std::vector< std::shared_ptr<BSPTriangle> > triangles_;

};

class BSPTree {

public:
  BSPTree( index_t level = 0 ) : level_(level) {}

  void build( BSPTriangles& triangles ) {

    if (triangles.nb() == 0) return;

    front_ = std::make_shared<BSPTree>(level_+1);
    back_  = std::make_shared<BSPTree>(level_+1);

    // find the triangle in the current list which would make a good node in the BSP tree
    // TODO something fancier...
    std::shared_ptr<BSPTriangle> root = triangles[0];
    index_t k = 0;
    while (root->ignore()) {
      root = triangles[k++];
      if (k == triangles.nb()) return;
    }
    triangles_.add(root);

    BSPTriangles front_triangles;
    BSPTriangles back_triangles;

    // loop through all the triangles and classify which are in front/back or split by the plane defined by the root
    for (index_t k = 0; k < triangles.nb(); k++) {

      if (triangles[k].get() == root.get()) continue;
      if (triangles[k]->ignore()) continue; // could have zero area

      bool coincident = triangles.classify( k , root->plane() , front_triangles , back_triangles );
      if (coincident) {
        // the triangle is exactly coincident with this plane, so we add it to this node of the tree
        triangles_.add( triangles[k] );
      }

    } // loop over BSP triangles

    // build the front and back trees
    front_->build(front_triangles);
    back_->build(back_triangles);

  }

private:
  BSPTriangles triangles_;
  std::shared_ptr<BSPTree> front_;
  std::shared_ptr<BSPTree> back_;
  index_t level_;
};

/*
static GLint gl2psFindRoot(GL2PSlist *primitives, GL2PSprimitive **root)
{
  GLint i, j, count, best = 1000000, idx = 0;
  GL2PSprimitive *prim1, *prim2;
  GL2PSplane plane;
  GLint maxp;

  if(!gl2psListNbr(primitives)){
    gl2psMsg(GL2PS_ERROR, "Cannot fint root in empty primitive list");
    return 0;
  }

  *root = *(GL2PSprimitive**)gl2psListPointer(primitives, 0);

  if(gl2ps->options & GL2PS_BEST_ROOT){
    maxp = gl2psListNbr(primitives);
    if(maxp > gl2ps->maxbestroot){
      maxp = gl2ps->maxbestroot;
    }
    for(i = 0; i < maxp; i++){
      prim1 = *(GL2PSprimitive**)gl2psListPointer(primitives, i);
      gl2psGetPlane(prim1, plane);
      count = 0;
      for(j = 0; j < gl2psListNbr(primitives); j++){
        if(j != i){
          prim2 = *(GL2PSprimitive**)gl2psListPointer(primitives, j);
          count += gl2psTestSplitPrimitive(prim2, plane);
        }
        if(count > best) break;
      }
      if(count < best){
        best = count;
        idx = i;
        *root = prim1;
        if(!count) return idx;
      }
    }
    // if(index) gl2psMsg(GL2PS_INFO, "GL2PS_BEST_ROOT was worth it: %d", index);
    return idx;
  }
  else{
    return 0;
  }
}

static void gl2psTraverseBspTree(GL2PSbsptree *tree, GL2PSxyz eye, GLfloat epsilon,
                                 GLboolean (*compare)(GLfloat f1, GLfloat f2),
                                 void (*action)(void *data), int inverse)
{
  GLfloat result;

  if(!tree) return;

  result = gl2psComparePointPlane(eye, tree->plane);

  if(GL_TRUE == compare(result, epsilon)){
    gl2psTraverseBspTree(tree->back, eye, epsilon, compare, action, inverse);
    if(inverse){
      gl2psListActionInverse(tree->primitives, action);
    }
    else{
      gl2psListAction(tree->primitives, action);
    }
    gl2psTraverseBspTree(tree->front, eye, epsilon, compare, action, inverse);
  }
  else if(GL_TRUE == compare(-epsilon, result)){
    gl2psTraverseBspTree(tree->front, eye, epsilon, compare, action, inverse);
    if(inverse){
      gl2psListActionInverse(tree->primitives, action);
    }
    else{
      gl2psListAction(tree->primitives, action);
    }
    gl2psTraverseBspTree(tree->back, eye, epsilon, compare, action, inverse);
  }
  else{
    gl2psTraverseBspTree(tree->front, eye, epsilon, compare, action, inverse);
    gl2psTraverseBspTree(tree->back, eye, epsilon, compare, action, inverse);
  }
}
*/

} // graphics

} // avro

#endif
