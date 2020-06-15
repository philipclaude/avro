// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef STRINGIFY_H
#define STRINGIFY_H

#include <string>
#include <sstream>
#include <iomanip>

namespace numpack 
{

  //Casts type T to string
  template<typename T>
  inline std::string stringify(const T& k)
  {
    std::stringstream ss;
    ss << std::boolalpha << k;
    return ss.str();
  }

  //A special case where a string is stringified
  template<>
  inline std::string stringify(const std::string& k) { return k; }

  template<typename T>
  inline std::string stringify(const T& k, const int& width, const char fill)
  {
    std::stringstream ss;
    ss << std::setfill(fill) << std::setw(width) << k;
    return ss.str();
  }

}

#endif //STRINGIFY_H
