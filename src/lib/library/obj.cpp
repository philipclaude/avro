#include "common/error.h"

#include "library/obj.h"

#include <fstream>

namespace ursa
{

namespace library
{

objFile::objFile( const std::string& filename ) :
  Topology<Simplex<Lagrange>>(vertices_,2),
  vertices_(3),
  filename_(filename)
{
  read();
}

void
objFile::read()
{
  FILE * file = fopen(filename_.c_str() , "r");
  if ( file == NULL )
  {
    printf("cannot open the file! %s\n",filename_.c_str());
    ursa_assert_not_reached;
  }

  real_t x[3],u[2],n[3];
  index_t t[3],nt[3],ut[3];
  char lineHeader[128];

  while (true)
  {


    // read the first word of the line
    int res = fscanf(file, "%s", lineHeader);
    if (res == EOF)
      break; // EOF = End Of File. Quit the loop.

    // else : parse lineHeader

    if ( strcmp( lineHeader, "v" ) == 0 )
    {
      fscanf(file, "%lg %lg %lg\n", &x[0] , &x[1] , &x[2] );
      vertices_.create(x);
    }
    else if ( strcmp( lineHeader, "vt" ) == 0 )
    {
      fscanf(file, "%lg %lg\n", &u[0] , &u[1] );
      u[1] = -u[1]; // invert V coordinate since we will only use DDS texture, which are inverted. Remove if you want to use TGA or BMP loaders.

      //TODO
    }
    else if ( strcmp( lineHeader, "vn" ) == 0 )
    {
      fscanf(file, "%lg %lg %lg\n", &n[0] , &n[1] , &n[2] );
      // TODO
    }
    else if ( strcmp( lineHeader, "f" ) == 0 )
    {
      int matches = fscanf(file, "%lu/%lu/%lu %lu/%lu/%lu %lu/%lu/%lu\n",
          &t[0], &ut[0], &nt[0], &t[1], &ut[1], &nt[1], &t[2], &ut[2], &nt[2] );
      if (matches == 9)
      {
        this->add( t , 3 );
        // TODO track other stuff....
      }
      else
        ursa_implement;
    }
    else
    {
      // probably a comment, eat up the rest of the line
      char comment[1024];
      fgets(comment, 1024, file);
    }

  }

  fclose(file);

}

} // library

} // ursa
