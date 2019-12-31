#ifndef avro_LIB_LIBRARY_EPS_H_
#define avro_LIB_LIBRARY_EPS_H_

#include "common/types.h"

#include <cstdlib>
#include <string>
#include <vector>

namespace avro
{

namespace library
{

class epsFile
{
public:
  epsFile();

  void add_triangles( const std::vector<real_t>& triangles , const std::vector<real_t>& colors );
  void add_edges( const std::vector<real_t>& edges , const std::vector<real_t>& colors );
  void add_points( const std::vector<real_t>& points );

  void write( const std::string& filename );

  void set_viewport( int* viewport );

private:
  FILE* fid_;
  void print_triangles() const;

  std::vector<int> viewport_;
  std::vector<real_t> triangles_;
  std::vector<real_t> triangle_colors_;

  std::vector<real_t> edges_;
  std::vector<real_t> edge_colors_;

  std::vector<real_t> points_;
};

} // library

} // avro


#endif
