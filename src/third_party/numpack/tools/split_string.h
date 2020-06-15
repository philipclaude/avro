// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef STRING_UTIL_H
#define STRING_UTIL_H

#include <string>
#include <vector>

namespace numpack 
{

// utility function for separating a string by a delimiter
std::vector<std::string>
split_string(const std::string &txt, const std::string &delim);


}

#endif //STRING_UTIL_H
