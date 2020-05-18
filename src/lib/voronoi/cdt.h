#ifndef avro_LIB_VORONOI_CDT_H_
#define avro_LIB_VORONOI_CDT_H_

#include "master/simplex.h"

#include "mesh/topology.h"

namespace avro
{

class PSC;

class ConstrainedDelaunayTriangulation : public Topology<Simplex>
{
public:
  ConstrainedDelaunayTriangulation( const Topology<Simplex>& psc );
  ConstrainedDelaunayTriangulation( const PSC& psc ); // will need to be broken into points, segments, triangles, etc.

  void compute();

  void recover_segments();
  void recover_facets();

private:
  const Topology<Simplex>& psc_;
};

} // avro

#endif
