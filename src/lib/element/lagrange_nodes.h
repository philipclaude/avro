#ifndef AVRO_LIB_ELEMENT_LAGRANGE_NODES_H_
#define AVRO_LIB_ELEMENT_LAGRANGE_NODES_H_

#include "common/error.h"
#include "avro_types.h"

namespace avro
{

template<typename type, int N, int P> struct LagrangeNodes;

template<typename type,int P>
struct LagrangeNodes<type,1,P> {
  static const std::vector<real_t> coord_s_;
};

template<typename type,int P>
struct LagrangeNodes<type,2,P> {
  static const std::vector<real_t> coord_s_;
  static const std::vector<real_t> coord_t_;
};

template<typename type,int P>
struct LagrangeNodes<type,3,P> {
  static const std::vector<real_t> coord_s_;
  static const std::vector<real_t> coord_t_;
  static const std::vector<real_t> coord_u_;
};

template<typename type,int P>
struct LagrangeNodes<type,4,P> {
  static const std::vector<real_t> coord_s_;
  static const std::vector<real_t> coord_t_;
  static const std::vector<real_t> coord_u_;
  static const std::vector<real_t> coord_v_;
};

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

}
#endif
