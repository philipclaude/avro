//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/error.h"
#include "common/tools.h"

#include "library/field.h"
#include "library/obj.h"

#include "mesh/dof.h"

#include <fstream>
#include <cstring>

namespace avro
{

namespace library
{

objFile::objFile( const std::string& filename ) :
  Topology<Simplex>(points_,2),
  points_(3),
  filename_(filename)
{
  this->set_sorted(false);
  read();
}

void
get_vertex( std::string& s , int& vertex , int& normal , int& texture )
{
  std::vector<std::string> pieces = split(s,"/");
  avro_assert(pieces.size()>0);

  vertex = normal = texture = -1;
  vertex = unstringify<int>(pieces[0]);

  if (pieces.size()>1)
    normal = unstringify<int>(pieces[1]);
  if (pieces.size()>2)
    texture = unstringify<int>(pieces[2]);
}

void
objFile::read()
{
  real_t x[3],u[2],n[3];
  index_t t[3];

  DOF<real_t> normals(3),uv(2);
  DOF<real_t> texture(2);
  std::vector< std::vector<index_t> > nt_vals,ut_vals;

  typedef Simplex Element_t;
  typedef Field<Element_t,std::vector<real_t>> FieldType;
  typedef Field<Element_t,std::vector<index_t>> FieldType_idx;

  std::shared_ptr<FieldType> normal_fld = std::make_shared<FieldType>(*this,1,CONTINUOUS);
  std::shared_ptr<FieldType> uv_fld = std::make_shared<FieldType>(*this,1,CONTINUOUS);

  std::shared_ptr<FieldType_idx> nt_fld = std::make_shared<FieldType_idx>(*this,0,DISCONTINUOUS);
  std::shared_ptr<FieldType_idx> ut_fld = std::make_shared<FieldType_idx>(*this,0,DISCONTINUOUS);

  this->set_sorted(false);

  // read the file
  std::ifstream file( filename_.c_str() , std::ios::in );
  if (!file.good()) {
    printf("could not open file %s\n",filename_.c_str());
  }
  avro_assert( file.good() );
  std::string line,token;
  getline(file,line);
  while (!file.eof())
  {
    size_t pos = line.find_first_of("#");
    if (pos!=std::string::npos)
      line = line.substr(0,pos);
    trim_string(line);

    if (line.length()>0)
    {
      std::istringstream line_stream( line );

      line_stream >> token;

      if (token == "v")
      {
        float X,Y,Z;
        line_stream >> X >> Y >> Z;
        x[0] = X;
        x[1] = Y;
        x[2] = Z;
        points_.create(x);
      }
      else if (token == "vt")
      {
        float S,T;
        line_stream >> S >> T;
        u[0] = S;
        u[1] = T;
        texture.add( u , 2 );
      }
      else if (token == "vn")
      {
        float NX,NY,NZ;
        line_stream >> NX >> NY >> NZ;
        n[0] = NX;
        n[1] = NY;
        n[2] = NZ;
        normals.add( n , 3 );
      }
      else if (token == "f")
      {
        std::vector<std::string> parts;
        while (line_stream.good())
        {
          std::string s;
          line_stream >> s;
          parts.push_back(s);
        }
        avro_assert_msg( parts.size()==3 , "only triangles supported" );

        int vertex0[3],vertex1[3],vertex2[3];
        get_vertex( parts[0] , vertex0[0] , vertex0[1] , vertex0[2] );
        get_vertex( parts[1] , vertex1[0] , vertex1[1] , vertex1[2] );
        get_vertex( parts[2] , vertex2[0] , vertex2[1] , vertex2[2] );

        t[0] = vertex0[0] -1;
        t[1] = vertex1[0] -1;
        t[2] = vertex2[0] -1;

        add(t,3);
      }
      else if (token=="usemtl")
      {
        printf("[warning] usemtl not supported!\n");
      }
      else if (token=="s")
      {
        printf("[warning] s token not supported!\n");
      }
      else
      {
        printf("bad token %s\n",token.c_str());
        //avro_assert_not_reached;
      }
    }
    getline(file,line);
  }
  file.close();

  printf("[warning] read facets...skipping normals, uv and texture...\n");

  return;

  // build the field with this connectivity and the same linear  element
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

  avro_assert_msg( normal_fld->nb_data()==points_.nb() ,
    "|fld| = %lu , |vertices| = %lu", normal_fld->nb(),points_.nb() );

  nt_fld->build();
  ut_fld->build();

  //avro_assert( nt_vals.size()==ut_vals.size() );
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

} // avro
