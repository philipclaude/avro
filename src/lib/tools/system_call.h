// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include <string>

namespace SANS
{

// Functions to for more robust system calls by trying multiple times
void system_call(const std::string& name, const std::string& exec, const std::string& args);
void wait_for_file(const std::string& filename);

}
