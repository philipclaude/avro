// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

/*
  This is essentially a copy-paste of the make_unique<T> class implemented within c++-14
  That it was not included in c++-11 is a travesty as it has a very useful role of allowing unique_ptrs to be assigned easily.
  There are a large number of places in the code where shared_ptr are used unnnecessarily, and unique_ptr would have been the better choice.
  This should go some way to helping to address this going forward. In short, unless you're passing the pointer around a lot,
  it should generally be a unique_ptr rather than a shared_ptr. This helps prevent resource leaks.

  Version by Stephan T. Lavavej (also known by STL) who originally proposed adding this function to C++14

  https://stackoverflow.com/questions/17902405/how-to-implement-make-unique-function-in-c11
  -- Hugh
*/

#ifndef MAKE_UNIQUE_H
#define MAKE_UNIQUE_H

#include <memory>
#include <type_traits>
#include <utility>

namespace tinymat  // This namespace is used so that it looks as similar as possible to the std::make_shared choice.
{
  template <typename T, typename... Args>
  std::unique_ptr<T> make_unique_helper(std::false_type, Args&&... args) { return std::unique_ptr<T>(new T(std::forward<Args>(args)...)); }

  template <typename T, typename... Args>
  std::unique_ptr<T> make_unique_helper(std::true_type, Args&&... args)
  {
    static_assert(std::extent<T>::value == 0,
    "make_unique<T[N]>() is forbidden, please use make_unique<T[]>().");

    typedef typename std::remove_extent<T>::type U;
    return std::unique_ptr<T>(new U[sizeof...(Args)]{std::forward<Args>(args)...});
  }

  template <typename T, typename... Args>
  std::unique_ptr<T> make_unique(Args&&... args) { return make_unique_helper<T>(std::is_array<T>(), std::forward<Args>(args)...); }

} // namespace std

#endif // MAKE_UNIQUE_H
