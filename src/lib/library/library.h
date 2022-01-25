//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_SRC_LIBRARY_LIBRARY_H_
#define AVRO_SRC_LIBRARY_LIBRARY_H_

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace avro
{

class Mesh;

class Library
{
public:
  static Library* instance;

  static Library* get()
  {
    if (!instance) instance = new Library;
    return instance;
  }

  void add_mesh( const std::string& mesh , const std::string& path=std::string() )
  {
    meshes_.push_back(mesh);
    if (!path.empty()) meshname2file_.insert({mesh,path});
  }

  void add_geometry( const std::string& geometry , const std::string& path=std::string() )
  {
    geometries_.push_back(geometry);
    if (!path.empty()) geometryname2file_.insert({geometry,path});
  }

  void add_metric( const std::string& metric )
  { metrics_.push_back(metric); }

  void add_points( const std::string& points )
  { points_.push_back(points); }

  void add_mesh_ptr( std::shared_ptr<Mesh> mesh )
  {
    mesh_ptr_.push_back(mesh);
  }

  std::string meshname2file( const std::string& mesh ) const
  {
    if (meshname2file_.find(mesh)==meshname2file_.end()) return "n/a";
    return meshname2file_.at(mesh)+"/"+mesh;
  }

  void reset()
  {
    meshes_ = meshes0_;
    metrics_ = metrics0_;
    geometries_ = geometries0_;
  }

  const std::vector<std::string>& meshes() const { return meshes_; }
  const std::vector<std::string>& metrics() const { return metrics_; }
  const std::vector<std::string>& geometries() const { return geometries_; }

private:
  const std::vector<std::string> meshes0_ = {};
  const std::vector<std::string> metrics0_ = {"Linear-3d","Polar1","Polar2","Linear-4d","Wave-4d"};
  const std::vector<std::string> geometries0_ = {"square","box","tesseract"};
  const std::vector<std::string> points0_ = {"Random"};

  std::vector<std::string> meshes_ = meshes0_;
  std::vector<std::string> metrics_ = metrics0_;
  std::vector<std::string> geometries_ = geometries0_;
  std::vector<std::string> points_ = points0_;

  std::map<std::string,std::string> meshname2file_;
  std::map<std::string,std::string> geometryname2file_;

  std::vector< std::shared_ptr<Mesh> > mesh_ptr_;

  Library() {}
};

} // avro

#endif
