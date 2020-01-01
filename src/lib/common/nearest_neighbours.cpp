#include "common/nearest_neighbours.h"
#include "common/tools.h"

#include "mesh/points.h"

#include "numerics/geometry.h"

#include <algorithm>

namespace avro
{

NearestNeighbours::NearestNeighbours( Points& _points , const index_t _knear ) :
  points_(_points) , knear_(_knear),
  neighbours_(TableLayout_Rectangular)
{
  if (_knear==0) knear_ = points_.nb();
  neighbours_.set_rank(knear_);
  compute();
}

void
NearestNeighbours::compute()
{
  neighbours_.clear();

  index_t nv = points_.nb();
  coord_t nd = points_.dim();

  for (index_t k=0;k<nv;k++)
  {
    // compute all the distances
    std::vector<real_t> d( nv , 0. );
    std::vector<index_t> nn( nv );
    for (index_t j=0;j<nv;j++)
    {
      d[j]  = numerics::distance2( points_[k] , points_[j] , nd );
      nn[j] = j;
    }

    // sort by distance
    std::sort( nn.begin() , nn.end() , SortBy<real_t>(d) );

    // add the result
    neighbours_.add( nn.data() , knear_ );
  }
}

void
NearestNeighbours::print() const
{
  for (index_t k=0;k<points_.nb();k++)
  {
    printf("neighbours[%d] = (",int(k));
    for (index_t j=0;j<knear_;j++)
      printf(" %d",int(neighbours_(k,j)));
    printf(" )\n");
  }
}

} // avro
