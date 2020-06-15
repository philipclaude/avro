// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// Specialization for adding vectors to gether.

#ifndef PLUS_STD_VECTOR_H
#define PLUS_STD_VECTOR_H

#include "SANSException.h"

#include <vector>
#include <functional>

namespace std
{
  template<typename _Tp>
  struct plus<std::vector<_Tp>> : public std::binary_function<std::vector<_Tp>, std::vector<_Tp>, std::vector<_Tp>>
  {
    std::vector<_Tp>
    operator()(const std::vector<_Tp>& __x, const std::vector<_Tp>& __y) const
    {
      SANS_ASSERT(__x.size() == __y.size());

      std::plus<_Tp> Op;
      std::vector<_Tp> tmp(__x.size());

      for (std::size_t i = 0; i < __x.size(); i++)
        tmp[i] = Op(__x[i], __y[i]);

      return tmp;
    }
  };
}

#endif //PLUS_STD_VECTOR_H
