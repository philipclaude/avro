// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef URSA_COMMON_STRINGIFY_H_
#define URSA_COMMON_STRINGIFY_H_

#include "common/error.h"
#include "common/types.h"

#include <algorithm>
#include <sstream>
#include <vector>

namespace ursa
{

template <typename T>
inline std::string
stringify( const T& x )
{
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

template <typename type> inline type unstringify(const std::string &txt);

template <>
inline int
unstringify<int>(const std::string &txt)
{
  return atoi(txt.c_str());
}

template <>
inline double
unstringify<double>(const std::string &txt)
{
  return atof(txt.c_str());
}

template <>
inline std::string
unstringify<std::string>(const std::string &txt)
{
  return txt;
}

template<>
inline index_t
unstringify<index_t>( const std::string &txt)
{
  return (index_t) atoi(txt.c_str());
}

template<>
inline coord_t
unstringify<coord_t>( const std::string &txt)
{
  return (coord_t) atoi(txt.c_str());
}

template<>
inline bool
unstringify<bool>( const std::string &txt)
{
  if (txt=="true" || txt=="True") return true;
  return false;
}

inline std::vector<std::string>
split(const std::string &txt,const std::string &delim)
{
  std::vector<std::string> strings;
  std::string::const_iterator start,end;
  if (delim.empty()) {
    strings.push_back(txt);
    return strings;
  }

  start = txt.begin();
  while (true) {
    end = search(start,txt.end(),delim.begin(),delim.end());
    std::string s(start,end);
    if (s.size()>0)
      strings.push_back(s);
    if (end==txt.end())
      break;
    start = end +delim.size();
  }
  return strings;
}

inline std::string
lowercase( const std::string& s )
{
  std::string S = s;
  std::transform( S.begin() , S.end() , S.begin() , ::tolower );
  return S;
}

} // ursa

#endif
