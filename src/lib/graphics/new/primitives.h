#ifndef AVRO_LIB_GRAPHICS_PRIMITIVES_H_
#define AVRO_LIB_GRAPHICS_PRIMITIVES_H_

#include "graphics/new/gl.h"

namespace avro
{

namespace graphics
{

#if 0

class SolutionPrimitive {

private:
	SolutionPrimitive( std::string name , coord_t order ) {
		nb_basis_ = (order+1)*(order+1)/2;
		sampler_index_ = 1; // colormap is always the first sampler in the shader
	}

	index_t nb_basis() const { return nb_basis_; }
	coord_t order() const { return order_; }

	index_t nb() const { return solution_.size() / nb_basis_; }

	void write( gl_index texture_buffer , gl_index texture ) {

		// bind the solution values to the texture buffer
	  GL_CALL( glBindBuffer( GL_TEXTURE_BUFFER , texture_buffer) );
	  GL_CALL( glBufferData( GL_TEXTURE_BUFFER , sizeof(float) * solution_.size() , solution_.data() , GL_STATIC_DRAW) );

		// generate a texture to hold the solution buffer
		GL_CALL( glGenTextures( 1 , &texture) );
		GL_CALL( glActiveTexture( GL_TEXTURE0 + 1 ) );
		GL_CALL( glBindTexture( GL_TEXTURE_BUFFER , texture) );
		GL_CALL( glTexBuffer( GL_TEXTURE_BUFFER , GL_R32F , texture_buffer ) );
	}

	void set_active( Shader& shader ) {
    glActiveTexture(GL_TEXTURE0 + 1 ); // colormap always comes first
		shader.setUniform( name_ , int(sampler_index_) );
	}

private:
	coord_t order_;
	index_t nb_basis_;
	gl_index buffer_;
	gl_index texture_;
	gl_index texture_index_;
	gl_index sampler_index_;
	std::string name_;

	std::vector<gl_float> solution_; // solution control points
};

class TrianglePrimitive {

public:

	TrianglePrimitive( Plot& plot , coord_t order ) :
	 	order_(order),
		plot_(plot),
		active_solution_("") {
		nb_basis_ = (order+1)*(order+2)/2;
	}

	~TrianglePrimitive() {
		// free the buffers
		glDeleteBuffers(1,&element_buffer_);
		glDeleteBuffers(1,&texture_buffer_);
		glDeleteTextures(1,&texture_);
	}

	void add_triangle( const index_t* p , index_t np );

	void write() {
		// make sure we are bound to the correct vao

		GL_CALL( glBindVertexArray( plot_.vao() ) );
	}

	void render( const std::string& solution_name ) {

		if (!visible_) return;

		// use the shader
		shader_.use();

		// set the appropriate solution to visualize
		std::map<std::string,std::shared_ptr<SolutionPrimitive>> it;
		it = solution_.find(solution_name);
		if (it != solution_.end()) {

			// tell the shader we want to plot a solution
			shader.setUniform("with_solution",1);

			// retrieve the solution primitive
			SolutionPrimitive& solution = *it->second.get();

			// if the solution is not currently active, we need to overwrite the data in the buffer
			if (solution_name != active_solution_) solution.write(texture_buffer_,texture);

			// set the solution texture to be active in the shader
			solution.set_active(shader_);
		}
		else {
			// tell the shader to use a constant color
			shader.setUniform("with_solution",0);

			// which color do we want to assign to this triangle patch?
			shader.setUniform("color", vec3(0.9,0.9,0.2) );
		}

		// set any uniforms for this patch of triangles
		// for example, maybe we want geometry triangles to glow, etc.

		// bind the buffer associated with the triangle patch
		GL_CALL( glBindBuffer( GL_ELEMENT_ARRAY_BUFFER , element_buffer_ ) );

		// set the appropriate patch parameter
		if (order_ > 1) {
			GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );

			// draw the elements
			GL_CALL( glDrawElements( GL_PATCHES , indices_.size() ) );
		}
		else {
			GL_CALL( glDrawElements( GL_TRIANGLES , indices_.size() ) );
		}
	}

private:
	Plot& plot_;
	index_t order_;
	index_t nb_basis_;
	std::vector<gl_index> indices_;
	std::map<std::string,SolutionPrimitive> solution_;
	std::string active_solution_;

	gl_index element_buffer_;
	gl_index texture_buffer_;
	gl_index texture_;
	bool  visible_;
};

class EdgePrimitive {
public:
	EdgePrimitive( coord_t qorder );

	void write();
	void render();

	std::vector<gl_index> indices_;
	gl_index buffer_;
	bool  visible_;
};

class PointPrimitive {

public:

	PointPrimitive( const Points& points );

	void write();
	void render() {
		// draw the points
	}

	gl_index buffer() const { return buffer_; }

private:
	std::vector<gl_float> coordinates_;
	gl_index buffer_;
	bool  visible_;
};
#endif

template<typename type> class Topology;

template<typename type>
class PrimitiveDecomposition {

public:
	PrimitiveDecomposition( const Topology<type>& topology ) :
		topology_(topology)
	{}

	void compute() {

		// count how many geometry entities there are
		std::set<index_t> entities;
		for (index_t k = 0; k < topology_.points().nb(); k++) {
			if (topology_.points().entity(k) == nullptr) continue;
			entities.insert( topology_.points().entity(k)->identifier() );
		}
		entities_.assign( entities.begin() , entities.end() );

		compute_triangles();
		compute_edges();
	}

	void compute_triangles();
	void compute_edges();

private:
  const Topology<type>& topology_;
	std::vector<real_t>  points_;
	std::vector<index_t> edges_;
	std::vector<index_t> triangles_;
	std::vector<int>     entities_;
};

} // graphics

} // avro

#endif
