// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

//----------------------------------------------------------------------------//
//  overloaded derivative operator
//  ref: derivify.h

#include <cmath>
#include <iostream>

#include "tools/SANSnumerics.h"     // Real
#include "SurrealD.h"


// I/O
/*
std::istream&
operator>>( std::istream& is, SurrealD& z )
{
  Real v = 0;
  Real d[10] = {0};
  char c = 0;
  int n = 0;
  bool done;

  is >> c;
  if (c == '(')
  {
    is >> v;

    is >> c;
    done = false;
    while (! done)
    {
      if (c != ')') is.clear(std::ios::badbit);
      if (c == ',')
      {
        is >> d[n]; n++;
      }
      else if (c == ')')
      {
        done = true;
      }
    }
  }
  else
  {
    is.putback(c);
    is >> v;
  }

  if (is) z = SurrealD(v, d, n);
  return is;
}
*/


// debug dump of private data
void
SurrealD::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "SurrealD: v_ = " << v_;
  if (N_ == 0)
  {
    if ( d_ == NULL )
      out << "  d_ = " << 0;
    else
      out << "  d_ = " << d_;
    out << "  N_ = " << N_ << std::endl;
  }
  else
  {
    out << "  d_[" << N_ << "] = (";
    for (unsigned int n = 0; n < N_-1; n++)
      out << d_[n] << ",";
    out << d_[N_-1] << ")" << std::endl;
  }
}
