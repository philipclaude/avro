#include "geometry/egads/object.h"
#include "geometry/egads/model.h"

#include "library/meshb.h"

#include <libmeshb/libmeshb7.h>

namespace luma
{

namespace library
{

template<class... Args>
void
get_line( int fid , int keyword , Args... args )
{
  int ret = GmfGetLin(fid,keyword,args...);
  luma_assert_msg( ret > 0 , "ret = %d", ret );
}

meshb::meshb( const std::string& filename , const EGADS::Model* model ) :
  Mesh(0),
  filename_(filename),
  model_(model)
{
  read();
}

template<typename type>
void
meshb::read_elements( int GmfType )
{
  int status=1;
  index_t nb;
  int indices[6] = {-1,-1,-1,-1,-1};
  real_t dvalues[4] = {0.,0.,0.};
  float fvalues[4] = {0.,0.,0.};
  index_t simplex[6] = {0,0,0,0,0};

  nb = index_t( GmfStatKwd( fid_ , GmfType ) );
  printf("--> reading %lu %lu-elements\n",nb,nv(GmfType)-1);
  for (index_t k=0;k<nb;k++)
  {
    // read the element type
    int ref = -1;
    if (GmfType==GmfTetrahedra)
    {
      status = GmfGetLin( fid_ , GmfTetrahedra , &indices[0] , &indices[1] , &indices[2] , &indices[3] , &indices[4] );
      ref = indices[4];
    }
    else if (GmfType==GmfTriangles)
    {
      status = GmfGetLin( fid_ , GmfTriangles , &indices[0] , &indices[1] , &indices[2] , &indices[3] );
      ref = indices[3];
    }
    else if (GmfType==GmfEdges)
    {
      status = GmfGetLin( fid_ , GmfEdges , &indices[0] , &indices[1] , &indices[2] );
      ref = indices[2];
    }
    else if (GmfType==GmfCorners)
    {
      status = GmfGetLin( fid_ , GmfCorners , &indices[0] );
      ref = 0;
    }
    else if (GmfType==GmfVertices)
    {
      int domain = -1;
      if (points_.dim()==2)
      {
        if (version_==1)
          status = GmfGetLin( fid_ , GmfVertices , &fvalues[0] , &fvalues[1] , &domain );
        else
          status = GmfGetLin( fid_ , GmfVertices , &dvalues[0] , &dvalues[1] , &domain  );
      }
      else if (points_.dim()==3)
      {
        if (version_==1)
          GmfGetLin( fid_ , GmfVertices , &fvalues[0] , &fvalues[1] , &fvalues[2] , &domain );
        else
          GmfGetLin( fid_ , GmfVertices , &dvalues[0] , &dvalues[1] , &dvalues[2] , &domain  );
      }

      if (version_==1)
      {
        for (coord_t d=0;d<3;d++)
          dvalues[d] = real_t(fvalues[d]);
      }
      points_.create( dvalues );
      continue; // skip the topology stuff below
    }

    luma_assert( status==1 );

    for (index_t j=0;j<nv(GmfType);j++)
      simplex[j] = index_t(indices[j]-1); // 1-indexed

    std::map<int,index_t>::iterator it = ref_index_.find(ref);
    if (it==ref_index_.end())
    {
      // we need to create a new topology
      std::shared_ptr<Topology<type>> topology =
            std::make_shared<Topology<type>>(points_,nv(GmfType)-1);

      ref_index_.insert( { ref , nb_topologies() } );
      add( topology );

      topology->add(simplex, nv(GmfType) );
    }
    else
    {
      retrieve<type>( it->second ).add( simplex , nv(GmfType) );
    }
  }

}

index_t
meshb::nv( const int GmfType ) const
{
  if (GmfType==GmfVertices) return 1;
  if (GmfType==GmfTetrahedra) return 4;
  if (GmfType==GmfTriangles) return 3;
  if (GmfType==GmfEdges) return 2;
  if (GmfType==GmfCorners) return 1;
  luma_assert_not_reached;
  return 0;
}

void
meshb::read()
{
  // open the file
  int dim;
  fid_ = GmfOpenMesh(filename_.c_str(),GmfRead,&version_,&dim);
  luma_assert_msg( fid_ , "could not open mesh file %s ",filename_.c_str() );

  printf("--> vertices dimension = %d\n",dim);

  points_.set_dim(coord_t(dim));
  points_.set_parameter_dim(coord_t(dim-1));

  // determine the simplex number of the mesh
  int number = -1;
  if ( GmfStatKwd( fid_ , GmfTetrahedra ) > 0 )
    number = 3;
  else if ( GmfStatKwd( fid_ , GmfTriangles ) > 0 )
    number = 2;
  else if ( GmfStatKwd( fid_ , GmfEdges ) > 0 )
    number = 1;
  else
  {
    printf("no simplices found, assuming vertex mesh\n");
    number = 0;
  }

  // first read the vertices
  luma_assert( GmfGotoKwd( fid_ , GmfVertices ) > 0 );
  read_elements<Simplex>( GmfVertices );

  if (number>=3)
  {
    // read the topologies of tetrahedra
    if ( GmfGotoKwd( fid_ , GmfTetrahedra ) > 0 )
      read_elements<Simplex>( GmfTetrahedra );
  }

  if (number>=2)
  {
    // read the topologies of triangles
    if ( GmfGotoKwd( fid_ , GmfTriangles ) > 0 )
      read_elements<Simplex>( GmfTriangles  );
  }
  if (number>=1)
  {
    // read the topologies of edges
    if ( GmfGotoKwd( fid_ , GmfEdges ) > 0 )
      read_elements<Simplex>( GmfEdges );
  }
  if (GmfGotoKwd(fid_,GmfCorners)>0 )
  {
    // read the corner topologies

    // create a single topology to hold the corners
    index_t cornerTopologyIndex = 0; //topology(0).nb_children();
    std::shared_ptr<Topology<Simplex>> corners = std::make_shared<Topology<Simplex>>(points_,0);

    // create the reference index for this topology to be looked up
    ref_index_.insert( std::pair<int,index_t>(0,cornerTopologyIndex) );

    // read the corners
    read_elements<Simplex>( GmfCorners );
  }

  if (model_!=nullptr)
  {
    index_t nb;
    int status;
    Entity* entity;

    // read the VerticesOnGeometricTriangles
    if ( GmfGotoKwd( fid_ , GmfVerticesOnGeometricTriangles ) > 0 )
    {
      nb = index_t( GmfStatKwd( fid_ , GmfVerticesOnGeometricTriangles ) );
      printf("--> reading %lu vertices on Triangles!\n",nb);
      int v,f;
      real_t u[3]; // u[0] = u, u[1] = v, u[2] = distance
      for (index_t k=0;k<nb;k++)
      {
        status = GmfGetLin( fid_ , GmfVerticesOnGeometricTriangles , &v , &f , u , u+1 , u+2 );
        luma_assert(status==1);
        entity = model_->find_entity( index_t(f) , FACE );
        luma_assert(entity->number()==2);
        if (entity==NULL) printf("ERROR entity not found :(\n");
        points_.set_entity(v-1,entity);
        points_.set_param(v-1,u);
      }
    }

    // read the VerticesOnGeometricEdges
    if ( GmfGotoKwd( fid_ , GmfVerticesOnGeometricEdges ) > 0 )
    {
      nb = index_t( GmfStatKwd( fid_ , GmfVerticesOnGeometricEdges ) );
      printf("--> reading %lu vertices on Edges!\n",nb);
      int v,e;
      real_t u[2]; // u[0] = parameter, u[1] = distance
      for (index_t k=0;k<nb;k++)
      {
        status = GmfGetLin( fid_ , GmfVerticesOnGeometricEdges , &v , &e , &u[0] , &u[1] );
        luma_assert(status==1);
        entity = model_->find_entity( index_t(e) , EDGE );
        if (entity==NULL) printf("ERROR entity not found :(\n");
        luma_assert(entity->number()==1);
        points_.set_entity(v-1,entity);
        points_.set_param(v-1,u);
      }
    }

    // read the VerticesOnGeometricVertices
    if ( GmfGotoKwd( fid_ , GmfVerticesOnGeometricVertices ) > 0 )
    {
      nb = index_t( GmfStatKwd( fid_ , GmfVerticesOnGeometricVertices ) );
      printf("--> reading %lu vertices on Nodes!\n",nb);
      for (index_t k=0;k<nb;k++)
      {
        // read the vertex and which EGADS node this corresponds to
        int v,n;
        status = GmfGetLin( fid_ , GmfVerticesOnGeometricVertices , &v , &n );
        luma_assert( status==1 );
        entity = model_->find_entity( index_t(n) , NODE );
        if (entity==NULL) printf("ERROR entity not found :(\n");
        //entity->print();
        luma_assert(entity->number()==0);
        points_.set_entity(v-1,entity);
      }
    }
  }

  GmfCloseMesh(fid_);

}

} // library

} // luma
