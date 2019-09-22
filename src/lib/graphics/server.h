#ifndef URSA_LIB_GRAPHICS_SERVER_H_
#define URSA_LIB_GRAPHICS_SERVER_H_

namespace ursa
{

namespace graphics
{

class Message;

template<typename Framework>
class Server
{
public:
  Server();

  void send( int client , Message& ) const;
  void receive( const Message& );

private:
  Framework& derived() { return static_cast<Framework>(*this); }
  const Framework& derived() const { return static_cast<const Framework>(*this); }
  
};

class WebServer : public Server<WebServer>
{

};

class LocalServer : public Server<LocalServer>
{

};

} // graphics

} // server

#endif
