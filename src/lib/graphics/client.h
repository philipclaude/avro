#ifndef URSA_LIB_GRAPHICS_CLIENT_H_
#define URSA_LIB_GRAPHICS_CLIENT_H_

namespace ursa
{

namespace graphics
{

template<typename Framework>
class Client
{

private:
  Framework& derived() { return static_cast<Framework>(*this); }
  const Framework& derived() { return static_cast<const Framework>(*this); }

};

class WebClient : public Client<WebClient>
{
};

class LocalClient : public Client<LocalClient>
{
};

} // graphics

} // ursa

#endif
