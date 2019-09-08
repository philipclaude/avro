// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef FILE_EXTENSION_H_
#define FILE_EXTENSION_H_

// we don't currently require this additional dependency and it's probably not worth it
#define BOOST_WITH_FILESYSTEM 0

#if BOOST_WITH_FILYESYSTEM
#include <boost/filesystem.hpp>
#endif

#include <string>

namespace SANS
{

#if BOOST_WITH_FILESYSTEM

// use the boost version
inline std::string
get_file_extension( const std::string& filename )
{
  return boost::filesystem::extension(filename);
}

#else

inline std::string
get_file_extension( const std::string& filename )
{
  std::string::size_type idx;
  idx = filename.rfind('.'); // find the '.' in reverse order
  if (idx!=std::string::npos)
    return filename.substr(idx+1);
  std::string dummy;
  return dummy; // empty
}

#endif

}

#endif
