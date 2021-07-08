#include "graphics/new/primitives.h"

#include "graphics/new/vao.h"

namespace avro
{

namespace graphics
{

void
TrianglePrimitive::draw() {

    //if (field_ != nullptr) field_->activate();

    // bind the buffer for the indices we want to draw
    GL_CALL( glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffer_ ) );
    GL_CALL( glPatchParameteri( GL_PATCH_VERTICES , nb_basis_ ) );
    GL_CALL( glDrawElements(GL_PATCHES, indices_.size() , GL_UNSIGNED_INT , 0 ) );
}

} // graphics

} // avro
