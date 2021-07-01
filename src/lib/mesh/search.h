//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_MESH_SEARCH_H_
#define avro_MESH_SEARCH_H_

#include "avro_types.h"

#include "mesh/boundary.h"

#include <memory> // shared_ptr
#include <time.h>
#include <vector>

namespace avro
{

template<typename type> class Neighbours;
template<typename type> class Topology;
class Points;

class PointCloud;
class KdTreeNd;

class PointSearch
{
public:
  PointSearch( Points& points );

  // insert a coordinate into the kd tree
  void insert( const real_t* x );

  // obtain the N closest points to x and store their indices in "points"
  void closestPoints( const real_t* x , std::vector<index_t>& points ) const;

private:
  Points& points_;

  std::shared_ptr<PointCloud> cloud_;
  std::shared_ptr<KdTreeNd>   kdtree_;

};

class BoundarySearch
{
public:
  BoundarySearch( const Points& points );

  index_t nearest( const real_t* p ) const;

  const Points& points() const { return points_; }

private:
  BoundaryPoints points_;
  PointSearch searcher_;
};

template<typename type>
class ElementSearch
{
public:
    ElementSearch( const Topology<type>& _topology );

    int find( const real_t* x , const index_t start );
    int step( const real_t* x , const index_t start );
    int brute( const real_t* x );

    index_t closest( const real_t* x , std::vector<real_t>& alpha ) const;

    void print() const;

    index_t nb_search() const { return nb_search_; }

    real_t timing() const { return real_t(time_)/real_t(CLOCKS_PER_SEC); }

    bool& brute() { return brute_; }

    const std::vector<index_t>& path() const { return path_; }

private:
    const Topology<type>& topology_;
    const Neighbours<type>& neighbours_;
    BoundarySearch boundary_;
    std::vector<bool> visited_;
    std::vector<index_t> path_;

    index_t depth_;
    index_t nb_brute_force_;
    index_t nb_search_;
    clock_t time_;
    index_t nb_steps_;
    bool brute_;

};

} // avro

#endif
