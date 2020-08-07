#ifndef AVRO_LIB_VORONOI_CELL_H_
#define AVRO_LIB_VORONOI_CELL_H_

#include "element/polytope.h"

#include "mesh/topology.h"

class VoronoiCell : public Topology<Polytope>
{
public:
  VoronoiCell( const index_t k , const NearestNeighbours& neighbours , const Topology<Polytope>& domain );

  void initialize();
  void compute();

private:
  void clip_by_bisector( index_t bj );
  void clip_edge( index_t e0 , index_t e1 );

};

class LaguerreDiagram : public Topology<Polytope>
{
public:
  LaguerreDiagram( const Topology<Polytope>& domain , Delaunay& delaunay );

  void compute();
};

#endif
