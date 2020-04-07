#ifndef avro_LIB_GRAPHICS_MANAGER_H_
#define avro_LIB_GRAPHICS_MANAGER_H_

#include "common/error.h"
#include "common/types.h"

#include "graphics/listener.h"

namespace avro
{

namespace graphics
{

class Primitive;
class SceneGraph;
class TransformFeedbackResult;

class GraphicsManager
{
public:
  virtual ~GraphicsManager() {}
  virtual void write( Primitive& primitive ) = 0;
  virtual void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr ) = 0;

  Listener& listener() { return listener_; }

private:
  Listener listener_;
};

class Vulkan_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  Vulkan_Manager() {}

  void write( Primitive& primitive ) { avro_implement; }
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_implement; }

  void create_shaders()
  { avro_implement; }

};

} // graphics

} // avro

#endif
