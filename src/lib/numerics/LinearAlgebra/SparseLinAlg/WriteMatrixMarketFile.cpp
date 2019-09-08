// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "WriteMatrixMarketFile.h"

namespace SANS
{
namespace SLA
{

void WriteMatrixMarketFile( const DLA::MatrixDView<Real>& A, std::ostream& file )
{
  //Write the banner
  file << "%%MatrixMarket matrix coordinate real general" << std::endl;
  file << A.m() << " " << A.n() << " " << A.m()*A.n() << std::endl;

  //Write out the matrix data
  //Add one to the row and column index as the file format is 1-based
  for ( int row = 0; row < A.m(); row++ )
    for ( int col = 0; col < A.n(); col++ )
      file << std::setprecision( 16 ) << row+1 << " " << col+1 << " " << A(row,col) << std::endl;
}

}
}
