#include "library/plots.h"

#include "element/simplex.h"

#include "mesh/points.h"

namespace avro
{

namespace library
{

template<typename type>
Plot<type>::Plot( Points& points ) :
  Topology<type>(points,0)
{
  for (index_t k=0;k<points.nb();k++)
    Table<index_t>::add( &k , 1 );
}

template class Plot<Simplex>;

} // library

} // avro
