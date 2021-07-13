#ifndef AVRO_LIB_GRAPHICS_PRIMITIVES_H_
#define AVRO_LIB_GRAPHICS_PRIMITIVES_H_

#include "common/error.h"

#include "graphics/gl.h"

namespace avro
{

class Points;

namespace graphics
{

class PointPrimitive {

public:

	PointPrimitive( const Points& points );

  ~PointPrimitive();

	void write();

  void draw( ShaderProgram& program );

	gl_index buffer() const { return buffer_; }

  index_t nb() const { return coordinates_.size() / 3; }
  const std::vector<gl_float>& coordinates() const { return coordinates_; }

  void print() const;

private:
  const Points& points_;
	std::vector<gl_float> coordinates_;
	gl_index buffer_;
	bool  visible_;
};

class EdgePrimitive {

public:
  EdgePrimitive( coord_t order );
  void add( const index_t* v , index_t nv );

  index_t nb() const { return indices_.size() / nb_basis_; }
  const std::vector<gl_index>& indices() const { return indices_; }

  void write();

  void draw();

  void print() const;

private:
  coord_t order_;
  index_t nb_basis_;
  gl_index buffer_;
  std::vector<gl_index> indices_;
};

class FieldPrimitive;

class TrianglePrimitive {
public:
  TrianglePrimitive( coord_t order );

  ~TrianglePrimitive();

  void add( const index_t* v , index_t nv );

  void print() const;

  index_t nb() const { return indices_.size() / nb_basis_; }
  const std::vector<gl_index>& indices() const { return indices_; }

  void write();

  void draw();

  const vec3& color() const { return color_; }
  void set_color( vec3 color ) { color_ = color; }
  bool& visible() { return visible_; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_index> indices_;
  gl_index buffer_;
  FieldPrimitive* field_;
  vec3 color_;
  bool visible_;
};

class FieldData {
public:
  FieldData( coord_t order );

  void add( real_t* f , index_t ndof );

  const std::vector<gl_float>& data() const { return data_; }

private:
  coord_t order_;
  index_t nb_basis_;
  std::vector<gl_float> data_;
};

class FieldPrimitive {
public:
  FieldPrimitive(bool nowrite=false);

  ~FieldPrimitive();

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

private:
  std::string active_field_;
  index_t active_rank_;
  std::map<std::string, std::shared_ptr<FieldData> > data_;
  gl_index texture_;
  gl_index buffer_;
};

} // graphics

} // avro

#endif
