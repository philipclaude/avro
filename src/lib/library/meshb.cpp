//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/egads/object.h"
#include "geometry/egads/model.h"

#include "library/meshb.h"

#include "mesh/boundary.h"

#include <libmeshb/libmeshb7.h>

namespace avro
{

namespace library
{

template<class... Args>
void
get_line( int fid , int keyword , Args... args )
{
  int ret = GmfGetLin(fid,keyword,args...);
  avro_assert_msg( ret > 0 , "ret = %d", ret );
}

meshb::meshb( const std::string& filename , const EGADS::Model* model ) :
  Mesh(0),
  filename_(filename),
  model_(model),
  main_topology_(nullptr)
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
    else
      avro_implement;

    avro_assert( ref >= 0 );
    avro_assert( status==1 );

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
  avro_assert_not_reached;
  return 0;
}

void
meshb::read()
{
  // open the file
  int dim;
  fid_ = GmfOpenMesh(filename_.c_str(),GmfRead,&version_,&dim);
  avro_assert_msg( fid_ , "could not open mesh file %s ",filename_.c_str() );

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

  set_number(number);
  printf("number = %u\n",number_);

  // first read the points
  avro_assert( GmfGotoKwd( fid_ , GmfVertices ) > 0 );
  read_elements<Simplex>( GmfVertices );

  main_topology_ = std::make_shared<Topology<Simplex>>(points_,number);

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
      printf("--> reading %lu points on Triangles!\n",nb);
      int v,f;
      real_t u[3]; // u[0] = u, u[1] = v, u[2] = distance
      for (index_t k=0;k<nb;k++)
      {
        status = GmfGetLin( fid_ , GmfVerticesOnGeometricTriangles , &v , &f , u , u+1 , u+2 );
        avro_assert(status==1);
        entity = model_->find_entity( index_t(f) , FACE );
        avro_assert(entity->number()==2);
        if (entity==NULL) printf("ERROR entity not found :(\n");
        points_.set_entity(v-1,entity);
        points_.set_param(v-1,u);
      }
    }

    // read the VerticesOnGeometricEdges
    if ( GmfGotoKwd( fid_ , GmfVerticesOnGeometricEdges ) > 0 )
    {
      nb = index_t( GmfStatKwd( fid_ , GmfVerticesOnGeometricEdges ) );
      printf("--> reading %lu points on Edges!\n",nb);
      int v,e;
      real_t u[2]; // u[0] = parameter, u[1] = distance
      for (index_t k=0;k<nb;k++)
      {
        status = GmfGetLin( fid_ , GmfVerticesOnGeometricEdges , &v , &e , &u[0] , &u[1] );
        avro_assert(status==1);
        entity = model_->find_entity( index_t(e) , EDGE );
        if (entity==NULL) printf("ERROR entity not found :(\n");
        avro_assert(entity->number()==1);
        points_.set_entity(v-1,entity);
        points_.set_param(v-1,u);
      }
    }

    // read the VerticesOnGeometricVertices
    if ( GmfGotoKwd( fid_ , GmfVerticesOnGeometricVertices ) > 0 )
    {
      nb = index_t( GmfStatKwd( fid_ , GmfVerticesOnGeometricVertices ) );
      printf("--> reading %lu points on Nodes!\n",nb);
      for (index_t k=0;k<nb;k++)
      {
        // read the vertex and which EGADS node this corresponds to
        int v,n;
        status = GmfGetLin( fid_ , GmfVerticesOnGeometricVertices , &v , &n );
        avro_assert( status==1 );
        entity = model_->find_entity( index_t(n) , NODE );
        if (entity==NULL) printf("ERROR entity not found :(\n");
        avro_assert(entity->number()==0);
        points_.set_entity(v-1,entity);
      }
    }
  }

  GmfCloseMesh(fid_);
}

#if 0
template<typename type>
void
meshb<type>::read_sol( const std::string& filename )
{
  // open the file
  int dim,status;
  int nb_sol,numberType,solSize , TypTab[GmfMaxTyp];
  float fvalues[GmfMaxTyp];
  real_t dvalues[GmfMaxTyp];

  solfile_ = GmfOpenMesh(filename.c_str(),GmfRead,&solVersion_,&dim);
  avro_assert_msg( solfile_ ,
    "could not open sol file %s ",filename.c_str() );

  printf("version = %d, dimension = %d\n",solVersion_,dim);

  // create a field whether this is attached at points or cells

  nb_sol = GmfStatKwd( solfile_ , GmfSolAtVertices , &numberType , &solSize , TypTab );
  printf("nb_sol = %d, numberType = %d, solSize = %d\n",nb_sol,numberType,solSize);

  avro_assert( nb_sol == int(this->points_.nb()) );
  avro_assert( solSize == int(dim*(dim+1)/2) );

  // create a field of tensors
  TargetMetric *metric = new TargetMetric(dim);

  avro_assert( GmfGotoKwd( solfile_ , GmfSolAtVertices ) > 0 );
  for (int k=0;k<nb_sol;k++)
  {

    // read the metric
    if (solVersion_==1)
    {
      status = GmfGetLin( solfile_ , GmfSolAtVertices , fvalues );
      for (index_t j=0;j<6;j++)
        dvalues[j] = real(fvalues[j]);
    }
    else
      status = GmfGetLin( solfile_ , GmfSolAtVertices , dvalues );

    avro_assert( status==1 );

    std::vector<real_t> data(dvalues,dvalues+solSize);
    numerics::SPDT<real_t> m(data);
    metric->add(m);

  }

  this->fields_.make("metric",*metric);
  this->setFields();

}
#endif

template<>
void
meshb::write( const Topology<Simplex>& topology , const std::vector<index_t>& refs )
{
  const coord_t num = topology.number();

  if (refs.size()!=1) avro_assert( topology.nb() == refs.size() );

  int gmfType;
  int indices[ GmfMaxTyp ];

  if (num==0) gmfType = GmfCorners;
  else if (num==1) gmfType = GmfEdges;
  else if (num==2) gmfType = GmfTriangles;
  else if (num==3) gmfType = GmfTetrahedra;
  else avro_assert_not_reached;

  index_t nb = 0;
  for (index_t k=0;k<topology.nb();k++)
  {
    if (topology.ghost(k)) continue;
    nb++;
  }

  GmfSetKwd( fid_ , gmfType , nb );

  for (index_t k=0;k<topology.nb();k++)
  {
    if (topology.ghost(k)) continue;

    // fill the buffer
    for (index_t j=0;j<topology.nv(k);j++)
      indices[j] = int( topology(k,j)+1-topology.points().nb_ghost() );

    if (refs.size() != topology.nb())
      indices[num+1] = refs[0];
    else
      indices[num+1] = refs[k];

    if (num==0)
      GmfSetLin( fid_ , gmfType , indices[0] );
    else if (num==1)
      GmfSetLin( fid_ , gmfType , indices[0] , indices[1] , indices[2] );
    else if (num==2)
      GmfSetLin( fid_ , gmfType , indices[0] , indices[1] , indices[2] , indices[3] );
    else if (num==3)
      GmfSetLin( fid_ , gmfType , indices[0] , indices[1] , indices[2] , indices[3] , indices[4] );
    else if (num==4)
      GmfSetLin( fid_ , gmfType , indices[0] , indices[1] , indices[2] , indices[3] , indices[4] , indices[5] );
    else avro_assert_not_reached;

  }

}

void
meshb::open( int dim , const std::string& filename )
{
  fid_ = GmfOpenMesh(filename.c_str(),GmfWrite,GmfDouble,dim);
  avro_assert( fid_ );
}

void
meshb::close()
{
  GmfCloseMesh(fid_);
}

void
meshb::write( Points& points )
{
  // write the points
  const coord_t dim = points.dim();
  GmfSetKwd( fid_ , GmfVertices , points.nb()-points.nb_ghost() );
  for (index_t k=0;k<points.nb();k++)
  {
    if (k<points.nb_ghost())
    {
      continue;
    }
    const real_t* x = points[k];
    if (dim==2)
      GmfSetLin( fid_ , GmfVertices , x[0], x[1] , points.body(k)+1 );
    else if (dim>=3)
      GmfSetLin( fid_ , GmfVertices , x[0], x[1] , x[2] , points.body(k)+1 );
    else
      avro_assert_not_reached;
  }
}

void
meshb::write( Mesh& mesh , const std::string& filename , bool with_bnd )
{

  printf("--> writing mesh to: %s\n",filename.c_str());
  int dim = mesh.points().dim();
  if (dim==4) dim = 3;
  open(dim,filename);

  // write the points
  write(mesh.points());

  // get all the topologies
  std::vector<const TopologyBase*> topologies;
  mesh.retrieve(topologies);

  index_t ref = 0;
  for (index_t k=0;k<topologies.size();k++)
  {
    if (topologies[k]->dummy()) continue;
    if (topologies[k]->number()!=mesh.number()) continue;

    if (topologies[k]->type_name()=="simplex")
      write<Simplex>( *static_cast<const Topology<Simplex>*>(topologies[k]) , {ref++} );
    else
      avro_implement;
  }

  // option to write the boundary discretization
  if (with_bnd)
  {
    // we don't necessarily know where all the elements are stored
    Topology<Simplex> T(mesh.points(),mesh.number());
    mesh.retrieve(T);
    Boundary<Simplex> boundary(T);
    boundary.extract();

    write<Simplex>( boundary.nodes() , boundary.node_entities() );
    write<Simplex>( boundary.edges() , boundary.edge_entities() );
    if (mesh.number()>2)
      write<Simplex>( boundary.triangles() , boundary.triangle_entities() );
  }

  if (!with_bnd)
  {
    GmfCloseMesh( fid_ );
    return;
  }

  // gather all the geometry information
  std::vector<index_t> nodeVertices;
  std::vector<index_t> edgeVertices;
  std::vector<index_t> faceVertices;
  for (index_t k=0;k<mesh.points().nb();k++)
  {
    Entity* e = mesh.points().entity(k);
    if (e==NULL) continue;
    if (e->number()==0) nodeVertices.push_back(k);
    else if (e->number()==1) edgeVertices.push_back(k);
    else if (e->number()==2) faceVertices.push_back(k);
  }

  // write the GmfVerticesOnGeometricVertices
  GmfSetKwd( fid_ , GmfVerticesOnGeometricVertices , nodeVertices.size() );
  for (index_t k=0;k<nodeVertices.size();k++)
  {
    Entity* e = mesh.points().entity(nodeVertices[k]);

    // get the index of this entity in the body
    index_t id = e->identifier();

    // write the information
    GmfSetLin( fid_ , GmfVerticesOnGeometricVertices ,
                int(nodeVertices[k]+1-mesh.points().nb_ghost()) , int(id) );
  }

  // write the GmfVerticesOnGeometricEdges
  GmfSetKwd( fid_ , GmfVerticesOnGeometricEdges , edgeVertices.size() );
  for (index_t k=0;k<edgeVertices.size();k++)
  {
    Entity* e = mesh.points().entity(edgeVertices[k]);

    // get the index of this entity in the body
    index_t id = e->identifier();

    // write the information
    // TODO retrieve the t parameter of this vertex on the Edge
    GmfSetLin( fid_ , GmfVerticesOnGeometricEdges ,
               int(edgeVertices[k]+1-mesh.points().nb_ghost()) , int(id) , 0. );
  }

  // write the GmfVerticesOnGeometricTriangles
  GmfSetKwd( fid_ , GmfVerticesOnGeometricTriangles , faceVertices.size() );
  for (index_t k=0;k<faceVertices.size();k++)
  {
    Entity* e = mesh.points().entity(faceVertices[k]);

    // get the index of this entity in the body
    index_t id = e->identifier();

    // write the information
    // TODO retrieve the parameters of this vertex on the Face
    GmfSetLin( fid_ , GmfVerticesOnGeometricTriangles ,
               int(faceVertices[k]+1-mesh.points().nb_ghost()) , int(id) , 0. , 0. );
  }

  close();
}

} // library

} // avro
