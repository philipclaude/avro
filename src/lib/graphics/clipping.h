#ifndef avro_LIB_GRAPHICS_CLIPPING_H_
#define avro_LIB_GRAPHICS_CLIPPING_H_

#include "avro_types.h"

#include "graphics/math.h"

#include "mesh/points.h"

namespace avro
{

namespace graphics
{

class GraphicsManager;

class ClippingPlane
{
public:
  ClippingPlane( coord_t dim , real_t length=1.0 );

  void set_coordinates( const real_t* bounding_box=nullptr );

  void initialize();

  void reset();

  void update(); // for opengl
  void update( real_t distance , real_t* angles , int dir ); // for wv

  bool visible( const Points& points , const index_t* v , index_t nv ) const;

  const Points& points() const { return transformed_points_; }

  void plot( GraphicsManager& manager , const real_t* focus , const mat4& transformation );
  void hide( GraphicsManager& manager );

  void flip_normal() { sign_ *= -1; }
  int sign() const { return sign_; }


  void append_transformation( const mat4& mr , const mat4& mt );

private:
  std::vector<real_t> normal_;
  std::vector<real_t> center_;

  Points points_;
  Points transformed_points_;
  real_t length_;
  int sign_;

  mat4 transformation_;

  std::vector<real_t>  plot_points_;
  std::vector<index_t> plot_triangles_;
  std::vector<real_t>  plot_colors_;

  real_t distance_;
  std::vector<real_t> u_;
  std::vector<real_t> v_;
  real_t angles_[2];
};


} // graphics

} // avro

#endif
