// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef PROGRESSBAR_H
#define PROGRESSBAR_H

#include <iostream>
#include "tools/SANSnumerics.h"     // Real

namespace numpack 
{

class ProgressBar
{
public:

  ProgressBar(const std::string& pretext, int max_levels = 10) :
    pretext_(pretext), max_levels_(max_levels), level_(0)
  {}

  void update(const Real fraction)
  {
    if (fraction >= ((Real)level_/(Real)max_levels_))
    {
      print(fraction);
      level_++;
    }
  }

  void end()
  {
    std::cout<<std::endl<<std::endl<<std::flush;
  }

protected:

  void print(const Real fraction)
  {
    const int max_length = 50;
    int bar_length = (int)(fraction*max_length);

    std::cout << "\r" << pretext_ << " [";
    for (int i=0; i < bar_length-1; i++)
        std::cout << "=";

    if (bar_length > 0)
      std::cout << ">";

    for (int i=0; i < (max_length - bar_length); i++)
        std::cout << " ";
    std::cout << "] " << (int)(fraction*100) << "%";
    std::cout << std::flush;
  }

  const std::string pretext_;
  int max_levels_;
  int level_;
};

}

#endif //PROGRESSBAR_H
