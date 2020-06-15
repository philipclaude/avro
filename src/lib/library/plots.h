#ifndef avro_LIB_LIBRARY_PLOTS_H_
#define avro_LIB_LIBRARY_PLOTS_H_

#include "mesh/topology.h"

namespace avro
{

class Points;

namespace library
{

template<typename type>
class Plot : public Topology<type>
{
public:
  Plot( Points& points );
  Plot( Topology<type>& topology );
};

} // library

} // avro

#endif
