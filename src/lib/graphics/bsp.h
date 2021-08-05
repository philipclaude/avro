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

} BSPPlane;

class BSPTriangle {
public:
  BSPTriangle( const vec3& p0 , const vec3& p1 , const vec3& p2 ) {
    points_[0] = p0;
    points_[1] = p1;
    points_[2] = p2;

    // initialize the plane of this triangle
    plane_.center = p0;
    plane_.normal = glm::normalize( glm::cross( p1 - p0 , p2 - p0 ) );

    plane_.normal.print();
  }

  const vec3& point( index_t j ) const { return points_[j]; }

  const BSPPlane& plane() const { return plane_; }

private:
  vec3 points_[3];
  BSPPlane plane_;
};

typedef struct {

  int side;
  std::shared_ptr<BSPTriangle> triangle0;
  std::shared_ptr<BSPTriangle> triangle1;
  std::shared_ptr<BSPTriangle> triangle2;
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

  void add( std::shared_ptr<BSPTriangle>& triangle ) {
    triangles_.push_back(triangle);
  }

  std::shared_ptr<BSPTriangle>& operator[] (index_t k) { return triangles_[k]; }

  void split( index_t k , const BSPPlane& plane , BSPClassification& classification  ) {

    // compute the side for all three vertices
    int s0 = plane.side( triangles_[k]->point(0) );
    int s1 = plane.side( triangles_[k]->point(1) );
    int s2 = plane.side( triangles_[k]->point(2) );

    if (s0 > 0 && s1 > 0 && s2 > 0) {
      classification.side = 1;
      return;
    }
    if (s0 < 0 && s1 < 0 && s2 < 0) {
      classification.side = -1;
      return;
    }

    // see if any vertices lie exactly on the plane
    int nb_coincident = 0;
    if (s0 == 0) nb_coincident++;
    if (s1 == 0) nb_coincident++;
    if (s2 == 0) nb_coincident++;

    if (nb_coincident == 3) {
      classification.side = 2; // the entire triangle is coincident
      return;
    }

    // there is an intersection
    classification.side = 0;
    if (nb_coincident == 1) {
      // determine if the other two are on opposite sides

      if (s0 == 0) {

        if (s1*s2 < 0) {
          // intersection
          avro_implement;
        }
        else
          classification.side = s1;
      }
      else if (s1 == 0) {
        if (s0*s2 < 0) {
          // intersection
          avro_implement;
        }
        else
          classification.side = s0;
      }
      else if (s2 == 0) {
        if (s1*s2 < 0) {
          // intersection
          avro_implement;
        }
        else
          classification.side = s1;
      }
      else
        avro_assert_not_reached;
    }
    else if (nb_coincident == 2) {
      // return the side of the one that is not coincident

      if (s0 != 0) {
        avro_assert( s1 == 0 && s2 == 0 );
        classification.side = s0;
      }
      else if (s1 != 0) {
        avro_assert( s0 == 0 && s2 == 0 );
        classification.side = s1;
      }
      else if (s2 != 0) {
        avro_assert( s0 == 0 && s1 == 0 );
        classification.side = s2;
      }
      else
        avro_assert_not_reached;
    }
    else {

      avro_assert( nb_coincident == 0 );

      // none of the vertices are coincident, and there is a general intersection
      avro_implement;

    }
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
    std::shared_ptr<BSPTriangle>& root = triangles[0];
    triangles.add(root);

    BSPTriangles front_triangles;
    BSPTriangles back_triangles;

    // loop through all the triangles and classify which are in front/back or split by the plane defined by the root
    for (index_t k = 0; k < triangles.nb(); k++) {

      printf("level %lu, triangle %lu\n",level_,k);

      BSPClassification classification;
      triangles.split( k , root->plane() , classification );

      if (classification.side == -1) {

        // triangle is on negative side of plane
        back_triangles.add( triangles[k] );
      }
      else if (classification.side == 1) {
        // triangle is on positive side of plane
        front_triangles.add( triangles[k] );
      }
      else if (classification.side == 0) {

        avro_assert( classification.triangle0 != nullptr );
        avro_assert( classification.triangle1 != nullptr );
        avro_assert( classification.triangle2 != nullptr );

        // triangle is split by plane
        // the first two triangles are on the negative side
        front_triangles.add( classification.triangle0 );
        back_triangles.add( classification.triangle1 );

        // the third triangle is on the positive side
        front_triangles.add( classification.triangle2 );
      }
      else if (classification.side == 2) {
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

static void gl2psBuildBspTree(GL2PSbsptree *tree, GL2PSlist *primitives)
{
  GL2PSprimitive *prim, *frontprim = NULL, *backprim = NULL;
  GL2PSlist *frontlist, *backlist;
  GLint i, idx;

  tree->front = NULL;
  tree->back = NULL;
  tree->primitives = gl2psListCreate(1, 2, sizeof(GL2PSprimitive*));
  idx = gl2psFindRoot(primitives, &prim);
  gl2psGetPlane(prim, tree->plane);
  gl2psAddPrimitiveInList(prim, tree->primitives);

  frontlist = gl2psListCreate(1, 2, sizeof(GL2PSprimitive*));
  backlist = gl2psListCreate(1, 2, sizeof(GL2PSprimitive*));

  for(i = 0; i < gl2psListNbr(primitives); i++){
    if(i != idx){
      prim = *(GL2PSprimitive**)gl2psListPointer(primitives,i);
      switch(gl2psSplitPrimitive(prim, tree->plane, &frontprim, &backprim)){
      case GL2PS_COINCIDENT:
        gl2psAddPrimitiveInList(prim, tree->primitives);
        break;
      case GL2PS_IN_BACK_OF:
        gl2psAddPrimitiveInList(prim, backlist);
        break;
      case GL2PS_IN_FRONT_OF:
        gl2psAddPrimitiveInList(prim, frontlist);
        break;
      case GL2PS_SPANNING:
        gl2psAddPrimitiveInList(backprim, backlist);
        gl2psAddPrimitiveInList(frontprim, frontlist);
        gl2psFreePrimitive(&prim);
        break;
      }
    }
  }

  if(gl2psListNbr(tree->primitives)){
    gl2psListSort(tree->primitives, gl2psTrianglesFirst);
  }

  if(gl2psListNbr(frontlist)){
    gl2psListSort(frontlist, gl2psTrianglesFirst);
    tree->front = (GL2PSbsptree*)gl2psMalloc(sizeof(GL2PSbsptree));
    gl2psBuildBspTree(tree->front, frontlist);
  }
  else{
    gl2psListDelete(frontlist);
  }

  if(gl2psListNbr(backlist)){
    gl2psListSort(backlist, gl2psTrianglesFirst);
    tree->back = (GL2PSbsptree*)gl2psMalloc(sizeof(GL2PSbsptree));
    gl2psBuildBspTree(tree->back, backlist);
  }
  else{
    gl2psListDelete(backlist);
  }

  gl2psListDelete(primitives);
}

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


static GLint gl2psSplitPrimitive(GL2PSprimitive *prim, GL2PSplane plane,
                                 GL2PSprimitive **front, GL2PSprimitive **back)
{
  GLshort i, j, in = 0, out = 0, in0[5], in1[5], out0[5], out1[5];
  GLint type;
  GLfloat d[5];

  type = GL2PS_COINCIDENT;

  for(i = 0; i < prim->numverts; i++){
    d[i] = gl2psComparePointPlane(prim->verts[i].xyz, plane);
  }

  switch(prim->type){
  case GL2PS_POINT :
    if(d[0] > GL2PS_EPSILON)       type = GL2PS_IN_BACK_OF;
    else if(d[0] < -GL2PS_EPSILON) type = GL2PS_IN_FRONT_OF;
    else                           type = GL2PS_COINCIDENT;
    break;
  default :
    for(i = 0; i < prim->numverts; i++){
      j = gl2psGetIndex(i, prim->numverts);
      if(d[j] > GL2PS_EPSILON){
        if(type == GL2PS_COINCIDENT)      type = GL2PS_IN_BACK_OF;
        else if(type != GL2PS_IN_BACK_OF) type = GL2PS_SPANNING;
        if(d[i] < -GL2PS_EPSILON){
          gl2psAddIndex(in0, in1, &in, i, j);
          gl2psAddIndex(out0, out1, &out, i, j);
          type = GL2PS_SPANNING;
        }
        gl2psAddIndex(out0, out1, &out, j, -1);
      }
      else if(d[j] < -GL2PS_EPSILON){
        if(type == GL2PS_COINCIDENT)       type = GL2PS_IN_FRONT_OF;
        else if(type != GL2PS_IN_FRONT_OF) type = GL2PS_SPANNING;
        if(d[i] > GL2PS_EPSILON){
          gl2psAddIndex(in0, in1, &in, i, j);
          gl2psAddIndex(out0, out1, &out, i, j);
          type = GL2PS_SPANNING;
        }
        gl2psAddIndex(in0, in1, &in, j, -1);
      }
      else{
        gl2psAddIndex(in0, in1, &in, j, -1);
        gl2psAddIndex(out0, out1, &out, j, -1);
      }
    }
    break;
  }

  if(type == GL2PS_SPANNING){
    *back = (GL2PSprimitive*)gl2psMalloc(sizeof(GL2PSprimitive));
    *front = (GL2PSprimitive*)gl2psMalloc(sizeof(GL2PSprimitive));
    gl2psCreateSplitPrimitive(prim, plane, *back, out, out0, out1);
    gl2psCreateSplitPrimitive(prim, plane, *front, in, in0, in1);
  }

  return type;
}



*/

} // graphics

} // avro

#endif
