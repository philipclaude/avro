//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "numerics/linear_algebra.h"
#include "numerics/mat.h"
#include "numerics/predicates.h"
#include "avro_types.h"

#include "mesh/points.h"

#include <triangle/predicates.h>

#include <climits>
#include <algorithm>

using namespace GEO;
using namespace GEO::PCK;

static coord_t dim1 = 10;

using namespace avro;

class SimplexPoints : public Points
{
public:
  SimplexPoints( coord_t ns , coord_t nd=0 ) :
    Points(nd)
  {
    if (nd==0 && ns>0) nd = ns;

    for (coord_t k=0;k<ns+1;k++)
    {
      std::vector<real_t> x( nd , 0. );

      if (k>0)
        x[k-1] = 1.;

      create(x);
    }
  }
};

UT_TEST_SUITE(predicates_suite)

UT_TEST_CASE(side1_test)
{
  GEO::PCK::initialize();

  GEO::Sign answers[6] = {GEO::POSITIVE,GEO::NEGATIVE};

  for (coord_t nd=2;nd<=dim1;nd++)
  {

    // create the sites
    index_t nsites = 2;
    std::vector<real_t> zi(nd,0.);
    std::vector<real_t> zj(nd,0.);

    zi[0] =  -.5;

    zj[0] =  .5;

    // loop through allocation permutations for the points to test SOS
    std::vector<int> perm(nsites);
    for (index_t k=0;k<nsites;k++)
      perm[k] = k;

    // simplex points (point)
    Points v(nd);
    std::vector<real_t> x0(nd,0.);
    v.create(x0);

    real_t* p[2];
    for (index_t k=0;k<nsites;k++)
      p[k] = new double[nd];

    std::sort(p,p+nsites);

    index_t permutation = 0;

    do
    {

      // create the delaunay points in the allocation order
      Points delaunay(nd);

      for (coord_t i=0;i<nd;i++)
      {
        p[ perm[0] ][i] = zi[i];
        p[ perm[1] ][i] = zj[i];
      }

      const real_t* q0 = v[0];

      GEO::Sign sign = side1_SOS( p[perm[0]],p[perm[1]],q0,nd );

      UT_ASSERT( sign!=GEO::ZERO );
      UT_ASSERT_EQUALS( sign , answers[permutation] );

      permutation++;

    } while( std::next_permutation( perm.begin(),perm.begin()+perm.size() ) );

    for (index_t k=0;k<nsites;k++)
      delete [] p[k];
  }
}
UT_TEST_CASE_END(side1_test)

UT_TEST_CASE(side2_test)
{

  GEO::PCK::initialize();

  GEO::Sign answers[6] = {GEO::POSITIVE,GEO::POSITIVE,GEO::POSITIVE,
                          GEO::NEGATIVE,GEO::NEGATIVE,GEO::NEGATIVE};

  for (coord_t nd=2;nd<=dim1;nd++)
  {

    // create the sites
    index_t nsites = 3;
    std::vector<real_t> zi(nd,0.);
    std::vector<real_t> zj(nd,0.);
    std::vector<real_t> zk(nd,0.);

    zi[0] =  .5;
    zi[1] = -.5;

    zj[0] = -.5;
    zj[1] =  .5;

    zk[0] =  .5;
    zk[1] =  .5;

    // loop through allocation permutations for the points to test SOS
    std::vector<int> perm(nsites);
    for (index_t k=0;k<nsites;k++)
      perm[k] = k;

    // simplex points (edge)
    SimplexPoints v(1,nd);

    real_t* p[3];
    for (index_t k=0;k<nsites;k++)
      p[k] = new double[nd];

    std::sort(p,p+nsites);

    index_t permutation = 0;

    do
    {

      // create the delaunay points in the allocation order
      Points delaunay(nd);

      for (coord_t i=0;i<nd;i++)
      {
        p[ perm[0] ][i] = zi[i];
        p[ perm[1] ][i] = zj[i];
        p[ perm[2] ][i] = zk[i];
      }

      const real_t* q0 = v[0];
      const real_t* q1 = v[1];

      GEO::Sign sign = side2_SOS( p[perm[0]],p[perm[1]],p[perm[2]],q0,q1,nd );

      UT_ASSERT( sign!=GEO::ZERO );
      UT_ASSERT_EQUALS( sign , answers[permutation] );

      permutation++;

    } while( std::next_permutation( perm.begin(),perm.begin()+perm.size() ) );

    for (index_t k=0;k<nsites;k++)
      delete [] p[k];

  }

}
UT_TEST_CASE_END(side2_test)

UT_TEST_CASE(side3_test)
{
  GEO::PCK::initialize();

  for (coord_t nd=3;nd<=dim1;nd++)
  {

    // create the sites
    index_t nsites = 4;
    std::vector<real_t> zi(nd,0.);
    std::vector<real_t> zj(nd,0.);
    std::vector<real_t> zk(nd,0.);
    std::vector<real_t> zl(nd,0.);

    zi[0] =  .5;
    zi[1] =  .5;
    zi[2] = -.5;

    zj[0] = -.5;
    zj[1] =  .3;
    zj[2] = 1.5;

    zk[0] = -.5;
    zk[1] = 1.6;
    zk[2] = -.5;

    zl[0] =  .5;
    zl[1] =  .21;
    zl[2] =  .5;

    matd<real_t> A(3,3);
    for (coord_t ii=0;ii<3;ii++)
    {
      A(ii,0) = zj[ii] - zi[ii];
      A(ii,1) = zk[ii] - zi[ii];
      A(ii,2) = zl[ii] - zi[ii];
    }
    avro_assert( numerics::det(A) != 0 );

    // loop through allocation permutations for the points to test SOS
    std::vector<int> perm(nsites);
    for (index_t k=0;k<nsites;k++)
      perm[k] = k;

    // simplex points (edge)
    SimplexPoints v(2,nd);

    real_t* p[4];
    for (index_t k=0;k<nsites;k++)
      p[k] = new double[nd];

    std::sort(p,p+nsites);

    index_t permutation = 0;

    do
    {

      // create the delaunay points in the allocation order
      Points delaunay(nd);

      for (coord_t i=0;i<nd;i++)
      {
        p[ perm[0] ][i] = zi[i];
        p[ perm[1] ][i] = zj[i];
        p[ perm[2] ][i] = zk[i];
        p[ perm[3] ][i] = zl[i];
      }

      const real_t* q0 = v[0];
      const real_t* q1 = v[1];
      const real_t* q2 = v[2];

      GEO::Sign sign = GEO::ZERO;
      try
      {
        sign = side3_SOS( p[perm[0]],p[perm[1]],p[perm[2]],p[perm[3]],q0,q1,q2,nd );
      }
      catch (...)
      {
        printf("failure for d = %u\n",nd);
      }

      UT_ASSERT( sign!=GEO::ZERO );
      UT_ASSERT_EQUALS( sign , GEO::NEGATIVE ); // all the signs should be negative for this geometry

      permutation++;

    } while( std::next_permutation( perm.begin(),perm.begin()+perm.size() ) );

    for (index_t k=0;k<nsites;k++)
      delete [] p[k];

  }
}
UT_TEST_CASE_END(side3_test)

UT_TEST_CASE(side4_test)
{
  GEO::PCK::initialize();

  for (coord_t nd=4;nd<=dim1;nd++)
  {

    // create the sites
    index_t nsites = 5;
    std::vector<real_t> zi(nd,0.);
    std::vector<real_t> zj(nd,0.);
    std::vector<real_t> zk(nd,0.);
    std::vector<real_t> zl(nd,0.);
    std::vector<real_t> zm(nd,0.);

    zi[0] =  .5;
    zi[1] =  .5;
    zi[2] =  .5;
    zi[3] =  -.5;

    zj[0] = -.5;
    zj[1] =  .5;
    zj[2] =  .5;
    zj[3] =  .5;

    zk[0] = -.5;
    zk[1] =  .5;
    zk[2] = -.5;
    zk[3] =  .5;

    zl[0] =  .5;
    zl[1] = -.5;
    zl[2] =  .5;
    zl[3] =  .5;

    zm[0] =  .5;
    zm[1] = -.5;
    zm[2] =  .5;
    zm[3] =  .5;

    // this is the only one that's important relative to zi
    zm[0] =  .5;
    zm[1] =  .5;
    zm[2] =  .5;
    zm[3] =  .5;

    matd<real_t> A(4,4);
    for (coord_t ii=0;ii<4;ii++)
    {
      A(ii,0) = zj[ii] - zi[ii];
      A(ii,1) = zk[ii] - zi[ii];
      A(ii,2) = zl[ii] - zi[ii];
      A(ii,3) = zm[ii] - zi[ii];
    }
    avro_assert( numerics::det(A) != 0 );

    // loop through allocation permutations for the points to test SOS
    std::vector<int> perm(nsites);
    for (index_t k=0;k<nsites;k++)
      perm[k] = k;

    // simplex points (edge)
    SimplexPoints v(3,nd);

    real_t* p[5];
    for (index_t k=0;k<nsites;k++)
      p[k] = new double[nd];

    std::sort(p,p+nsites);

    index_t permutation = 0;

    do
    {

      // create the delaunay points in the allocation order
      Points delaunay(nd);

      for (coord_t i=0;i<nd;i++)
      {
        p[ perm[0] ][i] = zi[i];
        p[ perm[1] ][i] = zj[i];
        p[ perm[2] ][i] = zk[i];
        p[ perm[3] ][i] = zl[i];
        p[ perm[4] ][i] = zm[i];
      }

      const real_t* q0 = v[0];
      const real_t* q1 = v[1];
      const real_t* q2 = v[2];
      const real_t* q3 = v[3];

      GEO::Sign sign = GEO::ZERO;
      try
      {
        sign = side4_SOS( p[perm[0]],p[perm[1]],p[perm[2]],p[perm[3]],p[perm[4]], q0,q1,q2,q3,nd );
      }
      catch(...)
      {
        printf("failure for d = %u\n",nd);
      }

      UT_ASSERT( sign!=GEO::ZERO );
      //UT_ASSERT_EQUALS( sign , GEO::NEGATIVE ); // all the signs should be negative for this geometry

      permutation++;

    } while( std::next_permutation( perm.begin(),perm.begin()+perm.size() ) );

    for (index_t k=0;k<nsites;k++)
      delete [] p[k];

  }
}
UT_TEST_CASE_END(side4_test)

UT_TEST_CASE(side5_test)
{
  GEO::PCK::initialize();

  for (coord_t nd=5;nd<=dim1;nd++)
  {

    // create the sites
    index_t nsites = 6;
    std::vector<real_t> zi(nd,0.);
    std::vector<real_t> zj(nd,0.);
    std::vector<real_t> zk(nd,0.);
    std::vector<real_t> zl(nd,0.);
    std::vector<real_t> zm(nd,0.);
    std::vector<real_t> zo(nd,0.);

    zi[0] =  .5;
    zi[1] =  .5;
    zi[2] =  .5;
    zi[3] =  .5;
    zi[4] = -.5;

    zj[0] = -.5;
    zj[1] =  .5;
    zj[2] =  .5;
    zj[3] =  .5;
    zj[4] =  .5;

    zk[0] = -.5;
    zk[1] =  .5;
    zk[2] =  .5;
    zk[3] = -.5;
    zk[4] =  .5;

    zl[0] =  .5;
    zl[1] = -.5;
    zl[2] =  .5;
    zl[3] =  .5;
    zl[4] =  .5;

    zm[0] =  .5;
    zm[1] =  .5;
    zm[2] = -.5;
    zm[3] =  .5;
    zm[4] =  .5;

    // this is the only one that's important relative to zi
    zo[0] =  .5;
    zo[1] =  .5;
    zo[2] =  .5;
    zo[3] =  .5;
    zo[4] =  .5;

    matd<real_t> A(5,5);
    for (coord_t ii=0;ii<5;ii++)
    {
      A(ii,0) = zj[ii] - zi[ii];
      A(ii,1) = zk[ii] - zi[ii];
      A(ii,2) = zl[ii] - zi[ii];
      A(ii,3) = zm[ii] - zi[ii];
      A(ii,4) = zo[ii] - zi[ii];
    }
    avro_assert( numerics::det(A) > 0 );

    // loop through allocation permutations for the points to test SOS
    std::vector<int> perm(nsites);
    for (index_t k=0;k<nsites;k++)
      perm[k] = k;

    // simplex points (edge)
    SimplexPoints v(4,nd);

    real_t* p[6];
    for (index_t k=0;k<nsites;k++)
      p[k] = new double[nd];

    std::sort(p,p+nsites);

    index_t permutation = 0;

    do
    {

      // create the delaunay points in the allocation order
      Points delaunay(nd);

      for (coord_t i=0;i<nd;i++)
      {
        p[ perm[0] ][i] = zi[i];
        p[ perm[1] ][i] = zj[i];
        p[ perm[2] ][i] = zk[i];
        p[ perm[3] ][i] = zl[i];
        p[ perm[4] ][i] = zm[i];
        p[ perm[5] ][i] = zo[i];
      }

      const real_t* q0 = v[0];
      const real_t* q1 = v[1];
      const real_t* q2 = v[2];
      const real_t* q3 = v[3];
      const real_t* q4 = v[4];

      GEO::Sign sign = GEO::ZERO;
      try
      {
        sign = side5_SOS( p[perm[0]],p[perm[1]],p[perm[2]],p[perm[3]],p[perm[4]],p[perm[5]],q0,q1,q2,q3,q4,nd );
      }
      catch(...)
      {
        printf("failure for d = %u\n",nd);
      }

      UT_ASSERT( sign!=GEO::ZERO );
      //UT_ASSERT_EQUALS( sign , GEO::NEGATIVE ); // all the signs should be negative for this geometry

      permutation++;

    } while( std::next_permutation( perm.begin(),perm.begin()+perm.size() ) );

    for (index_t k=0;k<nsites;k++)
      delete [] p[k];

  }
}
UT_TEST_CASE_END(side5_test)

UT_TEST_SUITE_END(predicates_suite)
