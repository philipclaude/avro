#ifndef GEOGRAM_POINTS_NN_SEARCH
#define GEOGRAM_POINTS_NN_SEARCH

#include "common/error.h"

#include "common.h"
#include "defs.h"
#include "numeric.h"

#include <memory>
#include <vector>

/**
 * \file geogram/points/nn_search.h
 * \brief Abstract interface for nearest neighbor searching
 */

namespace GEO {

    /*typedef unsigned long index_t;
    typedef unsigned short coord_t;
    typedef coord_t coord_index_t;
    typedef int signed_index_t;*/

    /**
     * \brief Abstract interface for nearest neighbor search algorithms.
     * \details
     * Given a point set in arbitrary dimension, creates a data
     * structure for efficient nearest neighbor queries.
     *
     * NearestNeighborSearch objects are created using method create() which
     * uses the Factory service. New search algorithms can be implemented and
     * registered to the factory using
     * geo_register_NearestNeighborSearch_creator().
     * \see NearestNeighborSearchFactory
     * \see geo_register_NearestNeighborSearch_creator
     */
    class NearestNeighborSearch {
    public:
        /**
         * \brief Creates a new search algorithm
         * \param[in] dimension dimension of the points (e.g., 3 for 3d)
         * \param[in] name name of the search algorithm to create:
         * - "ANN" - uses the standard ANN algorithm
         * - "BNN" - uses the optimized KdTree
         * - "default", uses the command line argument "algo:nn_search"
         * \retval nullptr if \p name is not a valid search algorithm name
         * \retval otherwise, a pointer to a search algorithm object. The
         * returned pointer must be stored in a NearestNeighborSearch_var that
         * does automatic destruction:
         * \code
         * NearestNeighborSearch_var nnsearch =
         *     NearestNeighborSearch::create(3, "ANN");
         * \endcode
         */
        static NearestNeighborSearch* create( unsigned short dimension, const std::string& name = "default" );

        /**
         * \brief Sets the points and create the search data structure.
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension()
         */
        virtual void set_points(index_t nb_points, const double* points);

        /**
         * \brief Tests whether the stride variant of set_points() is supported
         * \return true if stride different from dimension can be used
         *  in set_points(), false otherwise
         */
        virtual bool stride_supported() const;

        /**
         * \brief Sets the points and create the search data structure.
         * \details This variant has a stride parameter. Call
         * stride_supported() before to check whether it is supported.
         *
         * \param[in] nb_points number of points
         * \param[in] points an array of nb_points * dimension()
         * \param[in] stride number of doubles between two consecutive
         *  points (stride=dimension() by default).
         */
        virtual void set_points(
            index_t nb_points, const double* points, index_t stride
        );

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  coordinates.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as an array of dimension() doubles
         * \param[out] neighbors array of nb_neighbors index_t
         * \param[out] neighbors_sq_dist array of nb_neighbors doubles
         */
        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            const double* query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const = 0;


      	/**
      	 * \brief A structure to discriminate between the two
      	 *  versions of get_nearest_neighbors()
      	 */
      	struct KeepInitialValues {};

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  coordinates. Uses input neighbors and squared distance as
      	 *  an initialization.
      	 * \details Default implementation ignores the input values.
      	 *  Derived classes may have more efficient implementations.
               * \param[in] nb_neighbors number of neighbors to be searched.
               *  Should be smaller or equal to nb_points() (else it triggers
               *  an assertion)
               * \param[in] query_point as an array of dimension() doubles
               * \param[in,out] neighbors array of nb_neighbors index_t
               * \param[in,out] neighbors_sq_dist array of nb_neighbors doubles
      	 * \param[in] dummy a dummy parameter to discriminate between the
      	 *  two forms of get_nearest_neighbors()
         */
        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            const double* query_point,
            index_t* neighbors,
            double* neighbors_sq_dist,
	          KeepInitialValues dummy
        ) const;

        /**
         * \brief Finds the nearest neighbors of a point given by
         *  its index.
         * \details For some implementation, may be faster than
         *  nearest neighbor search by point coordinates.
         * \param[in] nb_neighbors number of neighbors to be searched.
         *  Should be smaller or equal to nb_points() (else it triggers
         *  an assertion)
         * \param[in] query_point as the index of one of the points that
         *  was inserted in this NearestNeighborSearch
         * \param[out] neighbors array of nb_neighbors index_t
         * \param[out] neighbors_sq_dist array of nb_neighbors doubles
         */
        virtual void get_nearest_neighbors(
            index_t nb_neighbors,
            index_t query_point,
            index_t* neighbors,
            double* neighbors_sq_dist
        ) const;

        /**
         * \brief Nearest neighbor search.
         * \param[in] query_point array of dimension() doubles
         * \return the index of the nearest neighbor from \p query_point
         */

        index_t get_nearest_neighbor(
            const double* query_point
        ) const {
            index_t result;
            double sq_dist;
            get_nearest_neighbors(1, query_point, &result, &sq_dist);
            avro_assert(signed_index_t(result) >= 0);
            return index_t(result);
        }

        /**
         * \brief Gets the dimension of the points.
         * \return the dimension
         */
        coord_index_t dimension() const {
            return dimension_;
        }

        /**
         * \brief Gets the number of points.
         * \return the number of points
         */
        index_t nb_points() const {
            return nb_points_;
        }

        /**
         * \brief Gets a point by its index
         * \param[in] i index of the point
         * \return a const pointer to the coordinates of the point
         */
        const double* point_ptr(index_t i) const {
            avro_assert(i < nb_points());
            return points_ + i * stride_;
        }

        /**
         * \brief Search can be exact or approximate. Approximate
         *  search may be faster.
         * \return true if nearest neighbor search is exact,
         *   false otherwise
         */
        bool exact() const {
            return exact_;
        }

        /**
         * \brief Search can be exact or approximate. Approximate
         *  search may be faster.
         * \param[in] x true if nearest neighbor search is exact,
         *   false otherwise. Default mode is exact.
         */
        virtual void set_exact(bool x);

    protected:
        /**
         * \brief Constructs a NearestNeighborSearch.
         * \param[in] dimension dimension of the points
         */
        NearestNeighborSearch(coord_index_t dimension);

        /**
         * \brief NearestNeighborSearch destructor
         */
        virtual ~NearestNeighborSearch();

    protected:
        coord_index_t dimension_;
        index_t nb_points_;
        index_t stride_;
        const double* points_;
        bool exact_;
    };

    /**
     * \brief A smart pointer that contains a NearestNeighborSearch object.
     * \relates NearestNeighborSearch
     */
    typedef std::shared_ptr<NearestNeighborSearch> NearestNeighborSearch_var;
}

#endif
