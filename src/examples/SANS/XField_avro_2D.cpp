// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#define XFIELD_AVRO_INSTANTIATE
#include "XField_avro_impl.h"

#include "Field/XFieldArea.h"

namespace SANS
{

template<>
void
XField_avro<PhysD2,TopoD2>::fill( const std::vector<double>& params )
{
  SANS_DEVELOPER_EXCEPTION("not implemented");
}

// explicit instantiation
template class XFieldBase_avro<XField<PhysD2,TopoD2>>;
template class XField_avro<PhysD2,TopoD2>;

}
