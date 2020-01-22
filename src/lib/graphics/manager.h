#ifndef avro_LIB_GRAPHICS_MANAGER_H_
#define avro_LIB_GRAPHICS_MANAGER_H_

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

class WV_Manager : public GraphicsManager
{
  template<typename type> friend class Application;

private:
  WV_Manager();

  void write( Primitive& primitive ) { avro_implement; }
  void draw( SceneGraph& scene , TransformFeedbackResult* feedback=nullptr )
  { avro_assert_not_reached; }

};

} // graphics

} // avro

#endif
