#include "programs.h"

#include "common/directory.h"

#include "geometry/egads/context.h"
#include "geometry/egads/model.h"

#include "graphics/application.h"

#include "library/meshb.h"

#include <stdio.h>

#include <egads.h>

extern "C"
{
int EG_exportModel( egObject* model , size_t* nbytes , char* stream[] );
}

namespace avro
{

namespace programs
{

int
convert( int nb_input , const char** inputs )
{
  int nb_expected = 2;
  if (nb_input<nb_expected || nb_input==-1)
  {
    printf("\t\tconvert [input] [output]\n");
    printf("\t\toptions:\n");
    printf("\t\t--> none (so far)\n");
    return 0;
  }

  std::string filename_input( inputs[0] );
  std::string filename_output( inputs[1] );

  // get the file extensions
  std::string ext_input = get_file_ext( filename_input );
  std::string ext_output = get_file_ext( filename_output );

  printf("ext's are %s -> %s\n",ext_input.c_str(),ext_output.c_str());

  if (ext_input == "egads" && ext_output == "legads") {
    // convert to EGADS lite
    #ifndef AVRO_NO_ESP
    EGADS::Context context;
    EGADS::Model model(&context,filename_input);

    size_t nbytes;
    char *stream;
    EG_exportModel( model.object() , &nbytes , &stream );

    FILE *fp;
    fp = fopen(filename_output.c_str(),"wb");
    avro_assert( fp!=nullptr );
    fwrite(stream,sizeof(char),nbytes,fp);
    fclose(fp);
    EG_free(stream);
    #else
    printf("cannot convert .egads to .legads without full ESP");
    #endif
  }
  else if ((ext_input == "mesh" || ext_input == "meshb") && (ext_output == "mesh" || ext_output == "meshb")) {
    library::meshb meshb_in(filename_input);
    meshb_in.write(meshb_in,filename_output,false);
  }
  else {
    avro_assert_not_reached;
  }


  return 0;
}

} // programs

} // avro
