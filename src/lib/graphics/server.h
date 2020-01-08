#ifndef avro_LIB_GRAPHICS_SERVER_H_
#define avro_LIB_GRAPHICS_SERVER_H_

namespace avro
{

namespace graphics
{

#if 0
class Server
{
public:
  virtual void run() = 0;
  virtual void receive( const std::string& msg );
  virtual void send( const std::string& msg );
};

namespace OpenGL
{
  class Plotter : public graphics::Server
  {
  public:
    void run();
    void receive( const std::string& msg );
    void send( const std::string& msg );
  };
}

namespace WebGL
{
  class Plotter : public graphics::Server
  {
  public:
    Server( index_t port=7681 );

    void run();
    void receive( const std::string& msg );
    void send( const std::string& msg );

  private:
    index_t port_;
  };
}
#endif

} // graphics

} // server

#endif
