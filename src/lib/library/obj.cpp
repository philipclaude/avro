#include "common/error.h"

#include "library/field.h"
#include "library/obj.h"

#include <fstream>
#include <cstring>

namespace ursa
{

namespace library
{

objFile::objFile( const std::string& filename ) :
  Topology<Simplex<Lagrange>>(vertices_,2),
  vertices_(3),
  filename_(filename)
{
  this->setSorted(false);
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
  index_t t[3],ut[3];
  char lineHeader[128];

  Data<real_t> normals,uv;
  std::vector< std::vector<index_t> > nt_vals,ut_vals;

  typedef Simplex<Lagrange> Shape_t;
  typedef Field<Shape_t,Shape_t,std::vector<real_t>> FieldType;
  typedef Field<Shape_t,Shape_t,std::vector<index_t>> FieldType_idx;

  std::shared_ptr<FieldType> normal_fld = std::make_shared<FieldType>(*this,1,CONTINUOUS);
  std::shared_ptr<FieldType> uv_fld = std::make_shared<FieldType>(*this,1,CONTINUOUS);

  printf("making discontinuous fields..\n");

  std::shared_ptr<FieldType_idx> nt_fld = std::make_shared<FieldType_idx>(*this,0,DISCONTINUOUS);
  std::shared_ptr<FieldType_idx> ut_fld = std::make_shared<FieldType_idx>(*this,0,DISCONTINUOUS);

  printf("reading file..\n");

  index_t line = 0;
  int matches = -1;
  while (true)
  {
    // read the first word of the line
    int res = fscanf(file, "%s", lineHeader);
    if (res == EOF)
      break; // we reached the end of the file

    if ( strcmp( lineHeader, "v" ) == 0 )
    {
      matches = fscanf(file, "%lg %lg %lg\n", &x[0] , &x[1] , &x[2] );
      ursa_assert( matches==3 );
      //printf("read vertex (%g,%g,%g)\n",x[0],x[1],x[2]);
      vertices_.create(x);
    }
    else if ( strcmp( lineHeader, "vt" ) == 0 )
    {
      matches = fscanf(file, "%lg %lg\n", &u[0] , &u[1] );
      ursa_assert( matches==2 );
      u[1] = -u[1]; // invert V coordinate since we will only use DDS texture, which are inverted. Remove if you want to use TGA or BMP loaders.
      uv.add( u , 2 );
      //printf("read vt = (%g,%g)\n",u[0],u[1]);
    }
    else if ( strcmp( lineHeader, "vn" ) == 0 )
    {
      matches = fscanf(file, "%lg %lg %lg\n", &n[0] , &n[1] , &n[2] );
      ursa_assert( matches==3 );
      normals.add( n , 3 );
    }
    else if ( strcmp( lineHeader, "f" ) == 0 )
    {
      //printf("read facet %s\n",lineHeader);
      int matches; /*= fscanf(file, "%lu/%lu/%lu %lu/%lu/%lu %lu/%lu/%lu\n",
          &t[0], &ut[0], &nt[0], &t[1], &ut[1], &nt[1], &t[2], &ut[2], &nt[2] );
      if (matches == 9)
      {
        for (coord_t d=0;d<3;d++)
          t[d]--;
        this->add( t , 3 );

        std::vector<index_t> NT(nt,nt+3);
        std::vector<index_t> UT(ut,ut+3);

        nt_vals.push_back(NT);
        ut_vals.push_back(UT);
        continue;
      }*/
      matches = fscanf(file, "%lu/%lu %lu/%lu %lu/%lu\n",
          &t[0], &ut[0], &t[1], &ut[1], &t[2], &ut[2] );
      if (matches==6)
      {
        for (coord_t d=0;d<3;d++)
          t[d]--;
        this->add( t , 3 );
        std::vector<index_t> UT(ut,ut+3);
        ut_vals.push_back(UT);
        continue;
      }
      else
      {
        printf("error reading line %lu\n",line);
        printf("t = (%lu,%lu,%lu)\n",t[0],t[1],t[2]);
        ursa_assert_msg( false , "matches = %d " ,matches );
      }
    }
    else
    {
      // probably a comment, eat up the rest of the line
      char comment[1024];
      fgets(comment, 1024, file);
    }
    line++;

  }
  fclose(file);

  printf("read file..skipping normals, uv and texture...\n");

  return;

  // build the field with this connectivity and the same linear master element
  normal_fld->build();
  uv_fld->build();

  // store the values
  for (index_t k=0;k<normals.nb();k++)
  for (index_t j=0;j<3;j++)
  {
    index_t idx = (*this)(k,j);
    std::vector<real_t> n( normals(idx) , normals(idx)+3 );
    (*normal_fld)(k,j) = n;

    std::vector<real_t> U( uv(idx) , uv(idx)+3 );
    (*uv_fld)(k,j) = U;

  }

  ursa_assert_msg( normal_fld->nb_data()==vertices_.nb() ,
    "|fld| = %lu , |vertices| = %lu", normal_fld->nb(),vertices_.nb() );

  nt_fld->build();
  ut_fld->build();

  //ursa_assert( nt_vals.size()==ut_vals.size() );
  for (index_t k=0;k<nt_vals.size();k++)
  {
    //(*nt_fld)(k,0) = nt_vals[k];
    (*ut_fld)(k,0) = ut_vals[k];
  }

  // add the field
  fields_.make( "vertex_normals" , normal_fld );
  fields_.make( "vertex_uv" , uv_fld );
  fields_.make( "triangle_normals" , nt_fld );
  fields_.make( "triangle_uv" , ut_fld );

}

} // library

} // ursa
