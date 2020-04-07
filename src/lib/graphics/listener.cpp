#include "graphics/listener.h"

#include <json/json.hpp>

namespace avro
{

namespace graphics
{

void
Listener::send( const json& request , json& response )
{
  std::string command = request.at("command");

  if (command=="ls")
  {
    std::vector<json> items;
    directory_.ls(items);
    response["ls-response"] = items;
  }
  if (command=="cd")
  {
    directory_.cd( request["data"] );
  }
}

} // graphics

} // avro
