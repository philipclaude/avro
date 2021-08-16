#include "nn_search.h"
#include "kd_tree.h"

/****************************************************************************/

namespace GEO {

    NearestNeighborSearch::NearestNeighborSearch(
        coord_index_t dimension ) :
        dimension_(dimension),
        nb_points_(0),
        stride_(0),
        points_(nullptr),
        exact_(true)
    {}

    void NearestNeighborSearch::get_nearest_neighbors(
      	index_t nb_neighbors,
      	const double* query_point,
      	index_t* neighbors,
      	double* neighbors_sq_dist,
      	KeepInitialValues ) const
    {
    	get_nearest_neighbors(
    	    nb_neighbors,
    	    query_point,
    	    neighbors,
    	    neighbors_sq_dist
    	);
    }

    void NearestNeighborSearch::get_nearest_neighbors(
        index_t nb_neighbors,
        index_t query_point,
        index_t* neighbors,
        double* neighbors_sq_dist ) const {
        get_nearest_neighbors(
            nb_neighbors,
            point_ptr(query_point),
            neighbors,
            neighbors_sq_dist
        );
    }

    void NearestNeighborSearch::set_points(
        index_t nb_points, const double* points )
    {
        nb_points_ = nb_points;
        points_ = points;
        stride_ = dimension_;
    }

    bool NearestNeighborSearch::stride_supported() const { return false; }

    void NearestNeighborSearch::set_points(index_t nb_points, const double* points, index_t stride)
    {
      if (stride == index_t(dimension()))
      {
        set_points(nb_points, points);
        return;
      }
      avro_assert(stride_supported());
      nb_points_ = nb_points;
      points_ = points;
      stride_ = stride;
    }

    void NearestNeighborSearch::set_exact(bool x) { exact_ = x; }

    NearestNeighborSearch::~NearestNeighborSearch() {}

    NearestNeighborSearch* NearestNeighborSearch::create( unsigned short dimension, const std::string& name_in )
    {
      std::string name = name_in;
      if (name == "ANN") return new AdaptiveKdTree(dimension);
      if (name == "BNN") return new BalancedKdTree(dimension);

      // default
      return new BalancedKdTree(dimension);
    }

} // GEO
