#ifndef AVRO_LIB_GRAPHICS_PRIMITIVES_H_
#define AVRO_LIB_GRAPHICS_PRIMITIVES_H_

#include "common/error.h"

#include "graphics/gl.h"
#include "graphics/math.h"

#include <map>
#include <memory>
#include <vector>

namespace avro
{

class Points;

namespace graphics
{

class PrimitiveBase {
public:
  gl_index& buffer()       { return buffer_; }
  gl_index  buffer() const { return buffer_; }

  index_t  memory() const { return memory_; }
  index_t& memory()       { return memory_; }

protected:
  gl_index buffer_;
  index_t  memory_;
};

class PointPrimitive : public PrimitiveBase {

public:

	PointPrimitive( const Points& points );

	void write();

  void draw( ShaderProgram& program );

  index_t nb() const { return coordinates_.size() / 3; }
  const std::vector<gl_float>& coordinates() const { return coordinates_; }

  void print() const;

private:
  const Points& points_;
	std::vector<gl_float> coordinates_;
	bool  visible_;
};

class EdgePrimitive : public PrimitiveBase {

public:
  EdgePrimitive( coord_t order );

  coord_t order() const { return order_; }
  void add( const index_t* v , index_t nv );

  index_t nb() const { return indices_.size() / nb_basis_; }
  const std::vector<gl_index>& indices() const { return indices_; }

  void write();

  void draw(bool with_tess=true);

  void print() const;

  bool visible() const { return visible_; }
  void set_visible(bool x) { visible_ = x; }
  bool& visible() { return visible_; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_index> indices_;
  bool visible_;
};

class TrianglePrimitive : public PrimitiveBase {
public:
  TrianglePrimitive( coord_t order );

  void add( const index_t* v , index_t nv );

  void print() const;

  index_t nb() const { return indices_.size() / nb_basis_; }
  const std::vector<gl_index>& indices() const { return indices_; }

  void write();

  void draw(bool with_tess=true);

  const vec3& color() const { return color_; }
  void set_color( vec3 color ) { color_ = color; }
  bool& visible() { return visible_; }
  coord_t order() const { return order_; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_index> indices_;
  vec3 color_;
  bool visible_;
};

class FieldData {
public:
  FieldData( coord_t order );

  void add( real_t* f , index_t ndof );

  const std::vector<gl_float>& data() const { return data_; }
  coord_t order() const { return order_; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_float> data_;
};

class FieldPrimitive : public PrimitiveBase {
public:
  FieldPrimitive() {}

  void set_active( const std::string& active , index_t rank = 0 ) {
    active_field_ = active;
    active_rank_  = rank;
  }

  void set_active_rank( index_t rank ) { active_rank_ = rank; }

  std::string field_name( const std::string& name , index_t rank ) const {
     return name + "-" + std::to_string(rank);
  }
  std::string active_name() const { return field_name(active_field_,active_rank_); }

  void add( const std::string& name , index_t rank , std::shared_ptr<FieldData> data );

  void write();

  void activate( ShaderProgram& shader );

  const std::map<std::string,std::shared_ptr<FieldData>>& data() const { return data_; }

  const FieldData& active() const { return *data_.at(active_name()).get(); }

  gl_index& texture()       { return texture_; }
  gl_index  texture() const { return texture_; }

  float umin() const { return umin_; }
  float umax() const { return umax_; }

private:
  std::string active_field_;
  index_t active_rank_;
  std::map<std::string, std::shared_ptr<FieldData> > data_;
  gl_index texture_;

  float umin_, umax_;
};

} // graphics

} // avro

#endif
