// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef SRC_LINEARALGEBRA_SPARSELINALG_WRITEPLAINVECTOR_H_
#define SRC_LINEARALGEBRA_SPARSELINALG_WRITEPLAINVECTOR_H_

#include <fstream>
#include <iostream>
#include <iomanip>

#include "ScalarVector.h"

namespace numpack 
{

namespace SLA
{

// ------------------------------------------------------------------------ //
// I/O
//
// write a plain (i.e. unformatted) vector to output stream
std::ostream& WritePlainVector( const ScalarVector& vec_plain, std::ostream& out );

// write formatted vector (e.g. SANS custom vector types) into a plain (i.e. unformatted) vector
template<class SparseVector>
std::ostream& WritePlainVector( const SparseVector& v, std::ostream& out )
{
  ScalarVector v_plain(v);
  WritePlainVector(v_plain, out);

  return out;
}
// write vector given filename
template<class SparseVector>
void
WritePlainVector( const SparseVector& v, const std::string& filename )
{
  std::fstream file(filename, std::ios::out);
  WritePlainVector(v, file);
}

} // namespace SLA

} // namespace numpack 

#endif /* SRC_LINEARALGEBRA_SPARSELINALG_WRITEPLAINVECTOR_H_ */
