#ifndef avro_LIB_GRAPHICS_LISTENER_H_
#define avro_LIB_GRAPHICS_LISTENER_H_

#include "common/directory.h"
#include "common/json.h"

namespace avro
{

namespace graphics
{


class Listener
{
public:
  Listener() {}

  void send( const json& request , json& response );

  std::string pwd() const { return directory_.pwd(); }

private:
  Directory directory_;
};

} // graphics

} // avro

#endif
