#include "library/meshb.h"

namespace luna
{

namespace meshb
{

#if 0
template<class... Args>
void
get_line( int fid , int keyword , Args... args )
{
  int ret = GmfGetLin(fid,keyword,args...);
  luna_assert_msg( ret > 0 , "ret = %d", ret );
}

template<typename type>
void
read_elements( int fid , int version , int GmfType , Topology<type>& topology )
{
  int status=1;
  index_t nb;
  int indices[6] = {-1,-1,-1,-1,-1};
  real dvalues[4] = {0.,0.,0.};
  float fvalues[4] = {0.,0.,0.};
  index_t simplex[6] = {0,0,0,0,0};

  nb = index_t( GmfStatKwd( fid , GmfType ) );
  printf("reading %lu %lu-elements\n",nb,nv(GmfType)-1);
  for (index_t k=0;k<nb;k++)
  {
    // read the element type
    int ref = -1;
    if (GmfType==GmfTetrahedra)
    {
      status = GmfGetLin( fid , GmfTetrahedra , &indices[0] , &indices[1] , &indices[2] , &indices[3] , &indices[4] );
      ref = indices[4];
    }
    else if (GmfType==GmfTriangles)
    {
      status = GmfGetLin( fid , GmfTriangles , &indices[0] , &indices[1] , &indices[2] , &indices[3] );
      ref = indices[3];
    }
    else if (GmfType==GmfEdges)
    {
      status = GmfGetLin( fid , GmfEdges , &indices[0] , &indices[1] , &indices[2] );
      ref = indices[2];
    }
    else if (GmfType==GmfCorners)
    {
      status = GmfGetLin( fid , GmfCorners , &indices[0] );
      ref = 0;
    }
    else if (GmfType==GmfVertices)
    {
      int domain = -1;
      if (topology.points().dim()==2)
      {
        if (meshVersion_==1)
          status = GmfGetLin( fid , GmfVertices , &fvalues[0] , &fvalues[1] , &domain );
        else
          status = GmfGetLin( fid , GmfVertices , &dvalues[0] , &dvalues[1] , &domain  );
      }
      else if (topology.points().dim()>=3)
      {
        if (meshVersion_==1)
          GmfGetLin( fid , GmfVertices , &fvalues[0] , &fvalues[1] , &fvalues[2] , &domain );
        else
          GmfGetLin( fid , GmfVertices , &dvalues[0] , &dvalues[1] , &dvalues[2] , &domain  );
      }

      if (version==1)
      {
        for (coord_t d=0;d<3;d++)
          dvalues[d] = real(fvalues[d]);
      }
      topology.points().create( dvalues );
      continue; // skip the topology stuff below
    }

    luna_assert( status==1 );

    for (index_t j=0;j<nv(GmfType);j++)
      simplex[j] = index_t(indices[j]-1); // 1-indexed

    int topoIndex = refIndex( ref );
    if (topoIndex<0)
    {
      // we need to create a new topology
      smart_ptr(Topology<Simplex>) t =
            smart_new(Topology<Simplex>)(this->vertices_,nv(GmfType)-1);

      std::string label;
      if (t->number()==0) label = "corners";
      else if (t->number()==1) label = "edges";
      else if (t->number()==2) label = "triangles";
      else if (t->number()==3) label = "tetrahedra";
      else label = stringify(t->number())+"topo";
      label += stringify(this->topology(0).nb_children());
      t->setName(label);

      refIndex_.insert(
        std::pair<int,index_t>( ref , this->topology(0).nb_children() ) );
      this->topology(0).addChild(t);

      t->add(simplex, nv(GmfType) );
    }
    else
    {
      this->topology(0).child( topoIndex )->add( simplex , nv(GmfType) );
    }
  }

}

template<typename type>
void
read( const std::string& filename , Topology<type>& topology )
{
  // open the file
  int dim,version;
  int fid = GmfOpenMesh(filename.c_str(),GmfRead,&version,&dim);
  luna_assert_msg( fid , "could not open mesh file %s ",filename.c_str() );

  printf("vertices dimension = %d\n",dim);

  topology.points().set_dim(coord_t(dim));
  topology.points().set_parameter_dim(coord_t(dim-1));

  // determine the simplex number of the mesh
  int number = -1;
  if ( GmfStatKwd( fid , GmfTetrahedra ) > 0 )
    number = 3;
  else if ( GmfStatKwd( fid , GmfTriangles ) > 0 )
    number = 2;
  else if ( GmfStatKwd( fid , GmfEdges ) > 0 )
    number = 1;
  else
  {
    printf("no simplices found, assuming vertex mesh\n");
    number = 0;
  }

  // the incoming topology is a dummy
  topology.set_number(number);
  topology.set_dummy(true);

  // first read the vertices
  luna_assert( GmfGotoKwd( fid , GmfVertices ) > 0 );
  read_elements( fid , GmfVertices );

  std::map<int,index_t> ref_index; // map from reference index to topology index
  if (number_>=3)
  {
    // read the topologies of tetrahedra
    if ( GmfGotoKwd( fid , GmfTetrahedra ) > 0 )
      read_elements( GmfTetrahedra );
  }

  if (number_>=2)
  {
    // read the topologies of triangles
    if ( GmfGotoKwd( fid , GmfTriangles ) > 0 )
      read_elements( GmfTriangles );
  }

  if (number_>=1)
  {
    // read the topologies of edges
    if ( GmfGotoKwd( fid , GmfEdges ) > 0 )
      read_elements( GmfEdges );
  }

  if (GmfGotoKwd(fid,GmfCorners)>0 )
  {
    // read the corner topologies

    // create a single topology to hold the corners
    index_t cornerTopologyIndex = this->topology(0).nb_children();
    std::shared_ptr<Topology<Simplex>> points =
                            smart_new(Topology<Simplex>)(topology.points(),0);
    topology.add_child(p);

    // create the reference index for this topology to be looked up
    ref_index.insert( std::pair<int,index_t>(0,cornerTopologyIndex) );

    // read the corners
    readElements( GmfCorners );
  }

  GmfCloseMesh(fid);

}

#endif

} // meshb

} // luna
