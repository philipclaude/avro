// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MatrixD.h"
#include "numpack/dense/static/MatrixS.h"

namespace numpack 
{
namespace DLA
{

template <>
void
MatrixDView<Real>::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "MatrixD<" << m_ << "," << n_ << ",Real>:" << std::endl;

  out << indent << "  data = ";
  out << std::endl << "{";
  for (int i = 0; i < m_; i++)
  {
    out << "{";
    for (int j = 0; j < n_; j++)
    {
      out << v_[i*n_+j];
      if (j < n_-1) out << ", ";
    }
    out << "}";
    if (i < m_-1) out << std::endl << " ";
  }
  out << "}" << std::endl;

}

template <>
void
MatrixDView< MatrixS<1,1,Real> >::dump( int indentSize, std::ostream& out ) const
{
  std::string indent(indentSize, ' ');
  out << indent << "MatrixD<" << m_ << "," << n_ << ",MatrixS<1,1,Real> >:" << std::endl;

  out << indent << "  data = ";
  out << std::endl << "{";
  for (int i = 0; i < m_; i++)
  {
    out << "{";
    for (int j = 0; j < n_; j++)
    {
      out << v_[i*n_+j](0,0);
      if (j < n_-1) out << ", ";
    }
    out << "}";
    if (i < m_-1) out << std::endl << " ";
  }
  out << "}" << std::endl;
}

}
}
