// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef FILE_SKIPLINES_H
#define FILE_SKIPLINES_H

#include <fstream>
#include <limits>

namespace numpack 
{

inline void
file_skiplines( std::ifstream& input, const int nlines )
{
  // skip lines in the ASCII file
  for (int i = 0; i < nlines; i++)
    input.ignore(std::numeric_limits<std::streamsize>::max(), input.widen('\n'));
}

} //namespace numpack 

#endif //FILE_SKIPLINES_H
