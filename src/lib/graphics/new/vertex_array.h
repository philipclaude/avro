#ifndef AVRO_LIB_GRAPHICS_VERTEX_ARRAY_
#define AVRO_LIB_GRAPHICS_VERTEX_ARRAY_

namespace avro
{

class TopologyBase;

namespace graphics
{

class VertexArray {};

class OpenGL_VertexArray : public VertexArray {

public:
	OpenGL_VertexArray( const TopologyBase& topology ) :
    topology_(topology)
  {}

	void create() {

		// a 4d should have been decomposed into 3d meshes before this
		// it could be a mesh of simplices for the boundary,
		// or a mesh of polyhedra for a slice
		avro_assert( topology_.number() <= 3);

		std::shared_ptr<PrimitiveDecomposition> decomposition = nullptr;
		if (topology.element().type() == "simplex")
			decomposition = std::make_shared<PrimitiveDecomposition<type>>( static_cast<Topology<Simplex>&>(topology_) );
		else
			avro_implement;

		const Fields& fields = topology.fields();

		decomposition->compute();

		// split up the facets into triangle primitives

		// define a solution texture on each triangle patch

		// compute all the edges of the topology

		// split up the edges into edge primitives
	}

	void render() {

		// bind the vertex array
		// assign the coordinate buffer to the vertex attribute
		GL_CALL( glBindVertexArray(vertex_array_) );
		GL_CALL( glBindBuffer( GL_ARRAY_BUFFER, points_.buffer()  ) );
		GL_CALL( glVertexAttribPointer( 0, 3, GL_FLOAT, GL_FALSE, 0, 0 ) ); // enable attribute in position 0 which holds coordinates
		GL_CALL( glEnableVertexAttribArray(0) );
		GL_CALL( glBindBuffer( GL_ARRAY_BUFFER , 0 ) );

		// set the matrices


		// render the triangle primitives
		for (index_t k = 0; k < triangles_.size(); k++) {
			triangles_[k].render();
		}

		// loop through the edge buffers
		for (index_t k = 0; k < edge_buffers_.size(); k++) {
			edges_[k].render();
		}

		// do we want to render points?
		if (points_visible_) {
			points_.render();
		}

	}

private:

  const TopologyBase& topology_;
	gl_index vertex_array_;

	PointPrimitive points_;
	std::vector<TrianglePrimitive> triangles_;
	std::vector<EdgePrimitive> edges_;
};

} // graphics

} // avro

#endif
