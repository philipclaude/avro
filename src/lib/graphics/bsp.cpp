#include "bsp.h"

namespace avro
{

namespace graphics
{

// this helped while debugging
#define BSP_CALL(X) try {X} catch( ... ) { printf("error calling %s on line %d of %s\n",#X,__LINE__,__FILE__); avro_assert_not_reached; }

BSPTriangle::BSPTriangle( const vec3& p0 , const vec3& p1 , const vec3& p2 ) {

  // save the vertices defining this triangle
  points_[0] = p0;
  points_[1] = p1;
  points_[2] = p2;

  // initialize edges to be non-ghost
  ghost_[0] = ghost_[1] = ghost_[2] = false;

  // initialize the plane of this triangle
  plane_.center = p0;
  plane_.normal = glm::cross( p1 - p0 , p2 - p0 );

  // check if the normal has a zero length (zero area triangle)
  ignore_ = false;
  real_t l = glm::norm( plane_.normal );
  if (l < BSP_TOL) {
    ignore_ = true;
  }
  else {
    plane_.normal = glm::normalize(plane_.normal);
  }

  if (std::isnan(plane_.normal[0]) || std::isnan(plane_.normal[1]) || std::isnan(plane_.normal[2])) {
    p0.print();
    p1.print();
    p2.print();
    ignore_ = true;
    //avro_assert_not_reached;
  }
}

void
BSPTriangles::build( const Plot& plot , const mat4& view_matrix , const mat4& projection_matrix , const mat4& screen_matrix ) {

  const VertexArrayObject& vao = plot.active_vao();
  const mat4& model_matrix = plot.model_matrix();

  // first map all the points using the transformation that would typically be done in the vertex shader
  // so we lose a bit of performance here (since it doesn't run on the GPU), but we don't need to use transform feedback
  std::vector<vec3> points;
  const std::vector<gl_float>& coordinates = vao.points().coordinates();
  mat4 transformation = /*screen_matrix * projection_matrix * view_matrix * */ model_matrix;
  for (index_t k = 0; k < coordinates.size()/3; k++) {

    vec4 p = { coordinates[3*k ] , coordinates[3*k+1] , coordinates[3*k+2] , 1.0f };
    vec4 q = transformation * p;

    vec3 v = { q[0] , q[1] , q[2] };
    points.push_back(v);
  }

  // add every triangle in the vao
  for (index_t k = 0; k < vao.nb_triangles(); k++) {

    const TrianglePrimitive& prim = vao.triangles(k);
    if (!prim.visible()) continue;

    const std::vector<gl_index>& indices = prim.indices();
    for (index_t j = 0; j < indices.size()/3; j++) {

      const vec3& p0 = points[ indices[3*j  ] ];
      const vec3& p1 = points[ indices[3*j+1] ];
      const vec3& p2 = points[ indices[3*j+2] ];

      triangles_.push_back( std::make_shared<BSPTriangle>(p0,p1,p2) );
    }
  }
}

bool
BSPTriangles::classify( index_t k , const BSPPlane& plane , BSPTriangles& front , BSPTriangles& back  ) const {

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

  // count how many vertices lie exactly on the plane
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

  if (nb_coincident == 2) {
    if (s0 != 0) {
      // add to the side of s0
      avro_assert( s1 == 0 && s2 == 0 );
      ADD_ONE( s0 );
    }
    else if (s1 != 0) {
      // add to the side of s1
      avro_assert( s0 == 0 && s2 == 0 );
      ADD_ONE( s1 );

    }
    else if (s2 != 0) {
      // add to the side of s2
      avro_assert( s0 == 0 && s1 == 0 );
      ADD_ONE( s2 );
    }
    else
      avro_assert_not_reached;
  }
  else if (nb_coincident == 1) {

    // one vertex is coincident to the plane
    // coordinates of the intersection point
    vec3 q;

    // determine if the other two are on opposite sides
    if (s0 == 0) {

      avro_assert( s1*s2 != 0 );
      if (s1*s2 < 0) { // different signs for vertices 1 & 2

        // split the triangle in two
        BSP_CALL( q = plane.intersect( triangles_[k]->point(1) , triangles_[k]->point(2) ); );

        std::shared_ptr<BSPTriangle> t0 = std::make_shared<BSPTriangle>( triangles_[k]->point(0) , triangles_[k]->point(1) , q );
        std::shared_ptr<BSPTriangle> t1 = std::make_shared<BSPTriangle>( q , triangles_[k]->point(2) , triangles_[k]->point(0) );
        ADD_TWO(s1,s2);
      }
      else {
        ADD_ONE(s1); // vertices 1 & 2 are on the same side, so add to the side of vertex 1
      }
    }
    else if (s1 == 0) {

      avro_assert( s0*s2 != 0 );
      if (s0*s2 < 0) { // different signs for vertices 0 & 2

        // split the triangle in two
        BSP_CALL( q = plane.intersect( triangles_[k]->point(0) , triangles_[k]->point(2) ); );

        std::shared_ptr<BSPTriangle> t0 = std::make_shared<BSPTriangle>( triangles_[k]->point(1) , triangles_[k]->point(2) , q );
        std::shared_ptr<BSPTriangle> t1 = std::make_shared<BSPTriangle>( q , triangles_[k]->point(0) , triangles_[k]->point(1) );
        ADD_TWO( s2 , s0 );
      }
      else {
        ADD_ONE( s0 ); // vertices 0 & 2 are on the same side, so add to the side of vertex 0
      }
    }
    else if (s2 == 0) {

      avro_assert( s0*s1 != 0 );
      if (s0*s1 < 0) { // different signs for vertices 0 & 1

        // split the triangle in two
        BSP_CALL( q = plane.intersect( triangles_[k]->point(0) , triangles_[k]->point(1) ); );

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
  else {

    avro_assert( nb_coincident == 0 );

    vec3 p0 = triangles_[k]->point(0);
    vec3 p1 = triangles_[k]->point(1);
    vec3 p2 = triangles_[k]->point(2);

    int fc;
    vec3 a, b, c;
    if (s0 * s1 >= 0) {
      a  = p0;
      b  = p1;
      c  = p2;
      fc = s2;
    }
    else if (s0 * s2 >= 0) {
      a  = p2;
      b  = p0;
      c  = p1;
      fc = s1;
    }
    else if (s1 * s2 >= 0) {
      a  = p1;
      b  = p2;
      c  = p0;
      fc = s0;
    }
    else {
      avro_assert_not_reached;
    }

    vec3 A, B;
    BSP_CALL( A = plane.intersect( a , c ); );
    BSP_CALL( B = plane.intersect( b , c ); );

    std::shared_ptr<BSPTriangle> T1, T2, T3;

    BSP_CALL( T1 = std::make_shared<BSPTriangle>(a,b,A); )
    BSP_CALL( T2 = std::make_shared<BSPTriangle>(b,B,A); )
    BSP_CALL( T3 = std::make_shared<BSPTriangle>(A,B,c); )

    // set the ghost edges
    T1->set_ghost(1,true);
    T2->set_ghost(1,true);
    T2->set_ghost(2,true);
    T3->set_ghost(0,true);
    T3->set_ghost(0,true);

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

void
BSPTree::build( BSPTriangles& triangles ) {

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

    // classify the triangle into the front or back subtree
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

void
BSPTree::get_triangles( const vec3& e , std::vector<BSPTriangle*>& triangles ) const {

  if (triangles_.nb() == 0) return;

  // orders the triangles from back to front
  BSPTriangle* root = triangles_[0].get();

  int side = root->plane().side(e);

  if (side < 0) {
    front_->get_triangles(e,triangles);
    for (index_t j = 0; j < triangles_.nb(); j++)
      triangles.push_back(triangles_[j].get());
    //triangles.push_back( root );
    back_->get_triangles(e,triangles);
  }
  else {
    back_->get_triangles(e,triangles);
    for (index_t j = 0; j < triangles_.nb(); j++)
      triangles.push_back(triangles_[j].get());
    //triangles.push_back(root);
    front_->get_triangles(e,triangles);
  }

}

} // avro

} // graphics
