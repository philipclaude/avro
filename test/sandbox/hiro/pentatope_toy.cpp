//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "element/simplex.h"

#include "library/ckf.h"
#include "library/tesseract.h"

#include "numerics/geometry.h"

using namespace avro;

UT_TEST_SUITE( pentatope_test_suite )

void
calculate_normal( const Points& points , const std::vector<index_t>& tet , std::vector<real_t>& normal ) {

  // see https://en.wikipedia.org/wiki/Cross_product#Generalizations, section on Multilinear algebra

  // retrieve the coordinates of the tetrahedron
  real_t x0,y0,z0,t0;
  real_t x1,y1,z1,t1;
  real_t x2,y2,z2,t2;
  real_t x3,y3,z3,t3;
  x0 = points[tet[0]][0]; y0 = points[tet[0]][1]; z0 = points[tet[0]][2]; t0 = points[tet[0]][3];
  x1 = points[tet[1]][0]; y1 = points[tet[1]][1]; z1 = points[tet[1]][2]; t1 = points[tet[1]][3];
  x2 = points[tet[2]][0]; y2 = points[tet[2]][1]; z2 = points[tet[2]][2]; t2 = points[tet[2]][3];
  x3 = points[tet[3]][0]; y3 = points[tet[3]][1]; z3 = points[tet[3]][2]; t3 = points[tet[3]][3];

  // compute the vectors from the first vertex to all other vertices
  real_t u1,u2,u3,u4;
  real_t v1,v2,v3,v4;
  real_t w1,w2,w3,w4;
  u1 = x1 - x0; u2 = y1 - y0; u3 = z1 - z0; u4 = t1 - t0;
  v1 = x2 - x0; v2 = y2 - y0; v3 = z2 - z0; v4 = t2 - t0;
  w1 = x3 - x0; w2 = y3 - y0; w3 = z3 - z0; w4 = t3 - t0;

  // calculate the product (derived using MATLAB's symbolic toolbox)
  /*
    syms u1 u2 u3 u4
    syms v1 v2 v3 v4
    syms w1 w2 w3 w4
    syms e1 e2 e3 e4

    A = [u1,u2,u3,u4;v1,v2,v3,v4;w1,w2,w3,w4;e1,e2,e3,e4];
    detA = simplify(det(A));

    dA_de1 = diff(detA,e1)
    dA_de2 = diff(detA,e2)
    dA_de3 = diff(detA,e3)
    dA_de4 = diff(detA,e4)
  */
  normal[0] = u2*v4*w3 - u2*v3*w4 + u3*v2*w4 - u3*v4*w2 - u4*v2*w3 + u4*v3*w2;
  normal[1] = u1*v3*w4 - u1*v4*w3 - u3*v1*w4 + u3*v4*w1 + u4*v1*w3 - u4*v3*w1;
  normal[2] = u1*v4*w2 - u1*v2*w4 + u2*v1*w4 - u2*v4*w1 - u4*v1*w2 + u4*v2*w1;
  normal[3] = u1*v2*w3 - u1*v3*w2 - u2*v1*w3 + u2*v3*w1 + u3*v1*w2 - u3*v2*w1;

  numerics::normalize( normal.data() , normal.size() );
}

UT_TEST_CASE( test1 )
{
  const real_t tol = 1e-10; // tolerance for floating point unit test assertions

  coord_t number = 4;   // topological number of the mesh
  coord_t dim = number; // ambient dimension of the mesh

  // this vector defines the number of points in each dimension of the mesh
  // (feel free to change these)
  std::vector<index_t> sizes(number);
  sizes[0] = 2; // number of points in x-direction
  sizes[1] = 2; // number of points in y-direction
  sizes[2] = 2; // number of points in z-direction
  sizes[3] = 2; // number of points in t-direction

  CKF_Triangulation mesh(sizes);
  mesh.orient();
  mesh.build_structures(); // builds the internal data structures such as pentatope-pentatope and vertex-pentatope adjacency information

  // loop through the elements of the mesh and print them out
  for (index_t k = 0; k < mesh.nb(); k++) {

    // note: mesh.nv(k) returns the number of vertices for element k
    UT_ASSERT( mesh.nv(k) == number+1 );

    printf("pentatope %lu = (%lu,%lu,%lu,%lu,%lu)\n", k , mesh(k,0), mesh(k,1) , mesh(k,2), mesh(k,3), mesh(k,4) );
  }

  // there are 24 pentatopes in each cube (division)
  index_t nb_pentatopes_per_cube = 24;
  index_t nb_cubes = 1;
  for (coord_t d = 0; d < number; d++)
    nb_cubes *= (sizes[d]-1);
  UT_ASSERT_EQUALS( nb_cubes*nb_pentatopes_per_cube , mesh.nb() ); // mesh.nb() is the number of elements in the mesh

  // extract the edges of the mesh
  std::vector<index_t> edges;
  mesh.get_edges( edges );
  index_t nb_edges = edges.size()/2; // two vertices per edge
  for (index_t k = 0; k < nb_edges; k++) {
    printf("edge %lu: (%lu,%lu)\n",k,edges[2*k],edges[2*k+1]);
  }

  // extract the boundary facets by traversing the mesh neighbours and looking for null neighbours
  std::vector<real_t> volume( dim , 0.0 );
  real_t S = 0.0;
  for (index_t k = 0; k < mesh.nb(); k++) {

    for (coord_t j = 0; j < number+1; j++) {

      int n = mesh.neighbours()(k,j);

      // extract the facet
      std::vector<index_t> tet;
      for (index_t i = 0; i < number+1; i++) {
        if (i == j) continue;
        tet.push_back(mesh(k,i));
      }

      // compute the normal
      std::vector<real_t> normal(dim);
      calculate_normal( mesh.points() , tet , normal );

      // flip the sign based on the ordering of the extracted facet
      real_t sign = (j % 2 == 0) ? -1.0 : 1.0;
      for (coord_t d = 0; d < dim; d++)
        normal[d] *= sign;

      // compute the absolute value of the volume of the tetrahedron
      std::vector<const real_t*> X(number);
      for (coord_t i = 0; i < number; i++)
        X[i] = mesh.points()[tet[i]];
      real_t v = numerics::volume_nd( X , dim );

      // add the contribution to the total volume
      for (coord_t d = 0; d < dim; d++)
        volume[d] += v * normal[d];

      // only perform the following checks for boundary tetrahedra
      if (n >= 0) continue; // not a boundary tetrahedron

      // add the contribution to the boundary volume
      S += v;

      print_inline(tet,"boundary tet: ");
      print_inline(normal,"\t normal: ");

      // the mesh is defined within [0,1]^4
      //  so normal vectors that point in -x should be at x = 0
      // and normal vectors that point in +x should be at x = 1
      // (the same idea for other dimensions as well)
      coord_t D;
      real_t  value = 0.0;
      if      (fabs(normal[0] + 1.0) < tol) D = 0, value = 0.0;
      else if (fabs(normal[0] - 1.0) < tol) D = 0, value = 1.0;
      else if (fabs(normal[1] + 1.0) < tol) D = 1, value = 0.0;
      else if (fabs(normal[1] - 1.0) < tol) D = 1, value = 1.0;
      else if (fabs(normal[2] + 1.0) < tol) D = 2, value = 0.0;
      else if (fabs(normal[2] - 1.0) < tol) D = 2, value = 1.0;
      else if (fabs(normal[3] + 1.0) < tol) D = 3, value = 0.0;
      else if (fabs(normal[3] - 1.0) < tol) D = 3, value = 1.0;
      else avro_assert_not_reached;

      for (coord_t i = 0; i < number; i++)
        UT_ASSERT_NEAR( mesh.points()[tet[i]][D] , value , tol );
    }
  }

  // the total signed volume should be zero
  for (coord_t d = 0; d < dim; d++)
    UT_ASSERT_NEAR( volume[d] , 0.0 , tol );

  // the boundary volume should be 8 (8 cubes of unit volume bound this tesseract domain)
  printf("boundary volume = %g\n",S);
  UT_ASSERT_NEAR( S , 8.0 , tol );


  // demo of how to compute the elements attached to an edge (shell)
  std::vector<index_t> shell;
  for (index_t k = 0; k < nb_edges; k++) {

    // edge endpoint vertices
    index_t p = edges[2*k];
    index_t q = edges[2*k+1];

    shell.clear();
    mesh.intersect( {p,q} , shell );
    print_inline( shell , "elements attached to edge " + std::to_string(k) );
  }

  // demo of how to compute the elements attached to a vertex (ball)
  std::vector<index_t> ball;
  for (index_t k = 0; k < mesh.points().nb(); k++) {

    ball.clear();
    mesh.intersect( {k} , ball );
    print_inline( ball , "elements attached to vertex " + std::to_string(k) );
  }

}
UT_TEST_CASE_END( test1 )


UT_TEST_SUITE_END( pentatope_test_suite )
