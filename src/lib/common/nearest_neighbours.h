#ifndef AVRO_GRAPH_NEIGHBOURS_H_
#define AVRO_GRAPH_NEIGHBOURS_H_

#include "common/table.h"
#include "common/types.h"

namespace avro
{

class Points;

class NearestNeighbours
{
public:
  NearestNeighbours( Points& points , const index_t knear=0 );

  void compute();

  index_t operator() ( index_t k , index_t j ) const
  {
    return neighbours_(k,j);
  }

  index_t nb( const index_t k ) const
  {
    return neighbours_.nv(k);
  }

  void print() const;

private:
  Points& points_;
  index_t knear_;
  Table<index_t> neighbours_;

};

} // avro

#endif