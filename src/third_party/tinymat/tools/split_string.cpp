// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include <string>
#include <algorithm>

#include "tools/split_string.h"

namespace tinymat 
{

// utility function for separating a string by a delimiter
std::vector<std::string>
split_string(const std::string &txt, const std::string &delim)
{
  std::vector<std::string> strings;
  std::string::const_iterator start,end;

  if (delim.empty())
  {
    // just return the original text if there is no delimiter
    strings.push_back(txt);
    return strings;
  }

  start = txt.begin();
  while (true)
  {
    end = std::search( start , txt.end() , delim.begin() , delim.end() );
    std::string s(start,end);
    if (s.size() > 0)
      strings.push_back(s);
    if (end == txt.end())
      break;
    start = end + delim.size();
  }

  return strings;
}

}
