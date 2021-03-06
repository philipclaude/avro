//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "element/element.h"
#include "element/lagrange_nodes.h"
#include "element/reference.h"

namespace avro
{

#define O3 1./3. // One over 3
#define T3 2./3. // Two over 3

#define Q4 1./4. // Quarter
#define H2 1./2. // Half
#define T4 3./4. // Three quarters

// line
template<> const std::vector<real_t> LagrangeNodes<Simplex,1,1>::coord_s_ = {0.0,1.0};
template<> const std::vector<real_t> LagrangeNodes<Simplex,1,2>::coord_s_ = {0.0,0.5,1.0};
template<> const std::vector<real_t> LagrangeNodes<Simplex,1,3>::coord_s_ = {0.0,1./3,2./3,1.0};
template<> const std::vector<real_t> LagrangeNodes<Simplex,1,4>::coord_s_ = {0.0,0.25,0.5,0.75,1.0};
template<> const std::vector<real_t> LagrangeNodes<Simplex,1,5>::coord_s_ = {0.0,0.2,0.4,0.6,0.8,1.0};

// triangle
template<> const std::vector<real_t> LagrangeNodes<Simplex,2,1>::coord_s_ = {0.0, 1.0, 0.0};
template<> const std::vector<real_t> LagrangeNodes<Simplex,2,1>::coord_t_ = {0.0, 0.0, 1.0};

template<> const std::vector<real_t> LagrangeNodes<Simplex,2,2>::coord_s_ = {0.0, 1.0, 0.0, 0.5, 0.0, 0.5};
template<> const std::vector<real_t> LagrangeNodes<Simplex,2,2>::coord_t_ = {0.0, 0.0, 1.0, 0.5, 0.5, 0.0};

template<> const std::vector<real_t> LagrangeNodes<Simplex,2,3>::coord_s_ = {0.0, 1.0, 0.0, 2./3., 1./3., 0.0, 0.0, 1./3., 2./3., 1./3.};
template<> const std::vector<real_t> LagrangeNodes<Simplex,2,3>::coord_t_ = {0.0, 0.0, 1.0, 1./3., 2./3., 2./3., 1./3., 0.0, 0.0, 1./3.};

template<> const std::vector<real_t> LagrangeNodes<Simplex,2,4>::coord_s_ = {0.0, 1.0, 0.0,
                                                                             0.75, 0.5, 0.25,
                                                                             0.0, 0.0, 0.0,
                                                                             0.25, 0.5, 0.75,
                                                                             0.25, 0.5, 0.25
                                                                            };
template<> const std::vector<real_t> LagrangeNodes<Simplex,2,4>::coord_t_ = {0.0, 0.0, 1.0,
                                                                             0.25, 0.5, 0.75,
                                                                             0.75, 0.5, 0.25,
                                                                             0.0, 0.0, 0.0,
                                                                             0.25, 0.25, 0.5
                                                                            };

// tetrahedron
#define COORD(P, c) template<> const std::vector<real_t> LagrangeNodes<Simplex,3,P>:: coord_##c## _

COORD(1, s) = {0, 1, 0, 0};
COORD(1, t) = {0, 0, 1, 0};
COORD(1, u) = {0, 0, 0, 1};

COORD(2, s) = {0, 1, 0, 0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5};
COORD(2, t) = {0, 0, 1, 0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0};
COORD(2, u) = {0, 0, 0, 1, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0};

COORD(3, s) = {0, 1, 0, 0, 0 , 0 , O3, T3, T3, O3, 0 , 0 , 0 , 0 , O3, T3, O3, 0 , O3, O3};
COORD(3, t) = {0, 0, 1, 0, T3, O3, 0 , 0 , O3, T3, T3, O3, 0 , 0 , 0 , 0 , O3, O3, 0 , O3};
COORD(3, u) = {0, 0, 0, 1, O3, T3, T3, O3, 0 , 0 , 0 , 0 , O3, T3, 0 , 0 , O3, O3, O3, 0 };
COORD(4,s) = {0, 1, 0, 0, 0 , 0 , 0 , Q4, H2, T4, T4, H2, Q4, 0 , 0 , 0 , 0 , 0 , 0 , Q4, H2, T4, H2, Q4, Q4, 0 , 0 , 0 , Q4, H2, Q4, Q4, Q4, H2, Q4};
COORD(4,t) = {0, 0, 1, 0, T4, H2, Q4, 0 , 0 , 0 , Q4, H2, T4, T4, H2, Q4, 0 , 0 , 0 , 0 , 0 , 0 , Q4, H2, Q4, Q4, Q4, H2, 0 , 0 , 0 , Q4, H2, Q4, Q4};
COORD(4,u) = {0, 0, 0, 1, Q4, H2, T4, T4, H2, Q4, 0 , 0 , 0 , 0 , 0 , 0 , Q4, H2, T4, 0 , 0 , 0 , Q4, Q4, H2, Q4, H2, Q4, Q4, Q4, H2, 0 , 0 , 0 , Q4};
#undef COORD

// pentatope
#define COORD(P, c) template<> const std::vector<real_t> LagrangeNodes<Simplex,4,P>:: coord_##c## _

// P = 1
COORD(1, s) = {0.0,  1.0,  0.0,  0.0,  0.0};
COORD(1, t) = {0.0,  0.0,  1.0,  0.0,  0.0};
COORD(1, u) = {0.0,  0.0,  0.0,  1.0,  0.0};
COORD(1, v) = {0.0,  0.0,  0.0,  0.0,  1.0};

// P = 2
COORD(2, s) = {0.0,  1.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.0,  0.5,  0.5,  0.5,  0.0,  0.0,  0.0};
COORD(2, t) = {0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5,  0.5,  0.0};
COORD(2, u) = {0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5,  0.0,  0.5,  0.0,  0.5};
COORD(2, v) = {0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.5,  0.0,  0.0,  0.5,  0.0,  0.5,  0.5};

// P = 3
COORD(3, s) = {0.0,    1.0,    0.0,    0.0,    0.0,    1./3.,  2./3.,  0.0,    0.0,    0.0,
//             0       1       2       3       4       5       6       7       8       9
               0.0,    0.0,    0.0,    2./3.,  1./3.,  2./3.,  1./3.,  2./3.,  1./3.,  0.0,
//             10      11      12      13      14      15      16      17      18      19
               0.0,    0.0,    0.0,    0.0,    0.0,    1./3.,  1./3.,  1./3.,  0.0,    0.0,
//             20      21      22      23      24      25      26      27      28      29
               0.0,    1./3.,  1./3.,  1./3.,  0.0};
//             30      31      32      33      34
COORD(3, t) = {0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    1./3.,  2./3.,  0.0,
//             0       1       2       3       4       5       6       7       8       9
               0.0,    0.0,    0.0,    1./3.,  2./3.,  0.0,    0.0,    0.0,    0.0,    2./3.,
//             10      11      12      13      14      15      16      17      18      19
               1./3.,  2./3.,  1./3.,  0.0,    0.0,    1./3.,  0.0,    0.0,    1./3.,  1./3.,
//             20      21      22      23      24      25      26      27      28      29
               0.0,    1./3.,  1./3.,  0.0,    1./3.};
//             30      31      32      33      34
COORD(3, u) = {0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,    1./3.,
//             0       1       2       3       4       5       6       7       8       9
               2./3.,  0.0,    0.0,    0.0,    0.0,    1./3.,  2./3.,  0.0,    0.0,    1./3.,
//             10      11      12      13      14      15      16      17      18      19
               2./3.,  0.0,    0.0,    2./3.,  1./3.,  0.0,    1./3.,  0.0,    1./3.,  0.0,
//             20      21      22      23      24      25      26      27      28      29
               1./3.,  1./3.,  0.0,    1./3.,  1./3.};
//             30      31      32      33      34
COORD(3, v) = {0.0,    0.0,    0.0,    0.0,    1.0,    0.0,    0.0,    0.0,    0.0,    0.0,
//             0       1       2       3       4       5       6       7       8       9
               0.0,    1./3.,  2./3.,  0.0,    0.0,    0.0,    0.0,    1./3.,  2./3.,  0.0,
//             10      11      12      13      14      15      16      17      18      19
               0.0,    1./3.,  2./3.,  1./3.,  2./3.,  0.0,    0.0,    1./3.,  0.0,    1./3.,
//             20      21      22      23      24      25      26      27      28      29
               1./3.,  0.0,    1./3.,  1./3.,  1./3.};
//             30      31      32      33      34

// P = 4
COORD(4, s) = {0.0,  1.0,  0.0,  0.0,  0.0,  0.25, 0.5,  0.75, 0.0,  0.0,
//             0     1     2     3     4     5     6     7     8     9
               0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.75, 0.5,  0.25,
//             10    11    12    13    14    15    16    17    18    19
               0.75, 0.5,  0.25, 0.75, 0.5,  0.25, 0.0,  0.0,  0.0,  0.0,
//             20    21    22    23    24    25    26    27    28    29
               0.0,  0.0,  0.0,  0.0,  0.0,  0.25, 0.5,  0.25, 0.25, 0.5,
//             30    31    32    33    34    35    36    37    38    39
               0.25, 0.25, 0.5,  0.25, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
//             40    41    42    43    44    45    46    47    48    49
               0.0,  0.0,  0.0,  0.5,  0.25, 0.25, 0.5,  0.25, 0.25, 0.5,
//             50    51    52    53    54    55    56    57    58    59
               0.25, 0.25, 0.0,  0.0,  0.0,  0.25, 0.0,  0.25, 0.25, 0.25};
//             60    61    62    63    64    65    66    67    68    69
COORD(4, t) = {0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.25, 0.5,
//             0     1     2     3     4     5     6     7     8     9
               0.75, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.25, 0.5,  0.75,
//             10    11    12    13    14    15    16    17    18    19
               0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.75, 0.5,  0.25, 0.75,
//             20    21    22    23    24    25    26    27    28    29
               0.5,  0.25, 0.0,  0.0,  0.0,  0.25, 0.25, 0.5,  0.0,  0.0,
//             30    31    32    33    34    35    36    37    38    39
               0.0,  0.0,  0.0,  0.0,  0.25, 0.5,  0.25, 0.25, 0.5,  0.25,
//             40    41    42    43    44    45    46    47    48    49
               0.0,  0.0,  0.0,  0.25, 0.5,  0.25, 0.25, 0.5,  0.25, 0.0,
//            50    51    52    53    54    55    56    57    58    59
               0.0,  0.0,  0.5,  0.25, 0.25, 0.25, 0.25, 0.0,  0.25, 0.25};
//            60    61    62    63    64    65    66    67    68    69
COORD(4, u) = {0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
//             0     1     2     3     4     5     6     7     8     9
              0.0,  0.25, 0.5,  0.75, 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
//             10    11    12    13    14    15    16    17    18    19
               0.25, 0.5,  0.75, 0.0,  0.0,  0.0,  0.25, 0.5,  0.75, 0.0,
//             20    21    22    23    24    25    26    27    28    29
               0.0,  0.0,  0.75, 0.5,  0.25, 0.0,  0.0,  0.0,  0.25, 0.25,
//             30    31    32    33    34    35    36    37    38    39
               0.5,  0.0,  0.0,  0.0,  0.25, 0.25, 0.5,  0.0,  0.0,  0.0,
//             40    41    42    43    44    45    46    47    48    49
               0.25, 0.5,  0.25, 0.25, 0.25, 0.5,  0.0,  0.0,  0.0,  0.25,
//             50    51    52    53    54    55    56    57    58    59
               0.5,  0.25, 0.25, 0.5,  0.25, 0.25, 0.25, 0.25, 0.0,  0.25};
//             60    61    62    63    64    65    66    67    68    69
COORD(4, v) = {0.0,  0.0,  0.0,  0.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
//             0     1     2     3     4     5     6     7     8     9
               0.0,  0.0,  0.0,  0.0,  0.25, 0.5,  0.75, 0.0,  0.0,  0.0,
//             10    11    12    13    14    15    16    17    18    19
               0.0,  0.0,  0.0,  0.25, 0.5,  0.75, 0.0,  0.0,  0.0,  0.25,
//             20    21    22    23    24    25    26    27    28    29
               0.5,  0.75, 0.25, 0.5,  0.75, 0.0,  0.0,  0.0,  0.0,  0.0,
//             30    31    32    33    34    35    36    37    38    39
               0.0,  0.25, 0.25, 0.5,  0.0,  0.0,  0.0,  0.25, 0.25, 0.5,
//             40    41    42    43    44    45    46    47    48    49
               0.25, 0.25, 0.5,  0.0,  0.0,  0.0,  0.25, 0.25, 0.5,  0.25,
//             50    51    52    53    54    55    56    57    58    59
               0.25, 0.5,  0.25, 0.25, 0.5,  0.25, 0.25, 0.25, 0.25, 0.0};
//             60    61    62    63    64    65    66    67    68    69

#undef COORD

template<typename type, int P>
void
GetLagrangeNodes<type,1,P>::get( std::vector<real_t>& nodes ) {

  const std::vector<real_t>& s = LagrangeNodes<type,1,P>::coord_s_;
  for (index_t k = 0; k < s.size(); k++) {
    nodes.push_back( s[k] );
    nodes.push_back( 1 - s[k] );
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
    nodes.push_back( 1 - s[k] - t[k] );
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
    nodes.push_back( 1 - s[k] - t[k] - u[k] );
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
    nodes.push_back( 1 - s[k] - t[k] - u[k] - v[k] );
  }
}

template<typename type>
void
get_lagrange_nodes( coord_t number , coord_t order , std::vector<real_t>& nodes ) {
  nodes.clear();
  if (number == 0) {
    nodes.push_back(1.0);
  }
  else if (number == 1) {
    static const int N = 1;
    if      (order == 1) GetLagrangeNodes<type,N,1>::get(nodes);
    else if (order == 2) GetLagrangeNodes<type,N,2>::get(nodes);
    else if (order == 3) GetLagrangeNodes<type,N,3>::get(nodes);
    else if (order == 4) GetLagrangeNodes<type,N,4>::get(nodes);
    else avro_implement;
  }
  else if (number == 2) {
    static const int N = 2;
    if      (order == 1) GetLagrangeNodes<type,N,1>::get(nodes);
    else if (order == 2) GetLagrangeNodes<type,N,2>::get(nodes);
    else if (order == 3) GetLagrangeNodes<type,N,3>::get(nodes);
    else if (order == 4) GetLagrangeNodes<type,N,4>::get(nodes);
    else avro_implement;
  }
  else if (number == 3) {
    static const int N = 3;
    if      (order == 1) GetLagrangeNodes<type,N,1>::get(nodes);
    else if (order == 2) GetLagrangeNodes<type,N,2>::get(nodes);
    else if (order == 3) GetLagrangeNodes<type,N,3>::get(nodes);
    else if (order == 4) GetLagrangeNodes<type,N,4>::get(nodes);
    else avro_implement;
  }
  else if (number == 4) {
    static const int N = 4;
    if      (order == 1) GetLagrangeNodes<type,N,1>::get(nodes);
    else if (order == 2) GetLagrangeNodes<type,N,2>::get(nodes);
    else if (order == 3) GetLagrangeNodes<type,N,3>::get(nodes);
    else if (order == 4) GetLagrangeNodes<type,N,4>::get(nodes);
    else avro_implement;
  }
  else
    avro_implement;
}

template void get_lagrange_nodes<Simplex>( coord_t , coord_t , std::vector<real_t>& );

} // avro
