#include "graphics/application.h"

#include <wsss.h>
#include <wsserver.h>

#include <sstream>
#include <unistd.h>

namespace avro
{

namespace graphics
{

void
Application<Web_Interface>::connect_client()
{
  json msg;
  msg["commands"] = {"make_tree"};
  msg["make_tree-data"] = scene_.menu();

  wv_broadcastText( const_cast<char*>(msg.dump().c_str()) );
}

void
Application<Web_Interface>::run()
{
  write();

  if (wv_startServer( 7681 , NULL , NULL , NULL , 0 , manager_.context() ) == 0)
  {
    printf("waiting for client...\n");
    while (wv_nClientServer(0)==0) {}
    usleep(500000);

    connect_client();

    while( wv_statusServer(0) ) usleep(500000);
  }
}

} // graphics

} // avro
