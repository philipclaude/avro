// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "geometry/body.h"
#include "geometry/model.h"
#include "geometry/tessellation.h"

#include "mesh/mesh.h"
#include "mesh/delaunay/mesher.h"
#include "mesh/topology.h"

extern "C"
{
typedef avro::real REAL;
typedef void VOID;
#define TRILIBRARY
#define ANSI_DECLARATORS
#include <triangle/triangle.h>
}

#define TETLIBRARY
#include <tetgen1.5.0/tetgen.h>

namespace avro
{

template<typename type>
Mesher<type>::Mesher( Model* _model , Mesh<type>& _mesh ) :
  model_(_model),
  mesh_(_mesh) ,
  boundaryFacets_(boundaryVertices_,_model->number())
{
 model_->listEntities(entities_);
}

template<typename type>
Mesher<type>::Mesher( Topology<type>& boundary , Mesh<type>& _mesh ) :
  model_(NULL),
  mesh_(_mesh) ,
  boundaryFacets_(boundaryVertices_,boundary.number())
{
  for (index_t k=0;k<boundary.vertices().nb();k++)
    entities_.push_back(boundary.vertices().entity(k));
  uniquify(entities_);
  printf("nb entities = %lu\n",entities_.size());
}

template<typename type>
int
Mesher<type>::label( Entity* entity ) const
{
  for (index_t j=0;j<entities_.size();j++)
  {
    if (entities_[j]==entity)
      return j;
  }
  return -1;
}

Triangle::Triangle( Model* _model , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(_model,_mesh)
{
  setDimension(2);

  // tessellate the model
  ModelTessellation tess(*model_);

  // compute the expected volume as the integral of x.n over the boundary faces
  expectedVolume_ = tess.volume();

  // retrieve the boundary facets
  boundaryFacets_.setSorted( false ); // do not sort the triangles when we get them from the EGADS tessellation
  tess.retrieveElements( 1 , boundaryFacets_ );

  // retrieve the boundary vertices from the model tessellation
  tess.vertices().copy( boundaryVertices_ );

  // retrieve the internal points from the bodies
  for (index_t k=0;k<tess.nb_internal();k++)
    internalVertices_.create( tess.internalVertex(k) );

}

Triangle::Triangle( Topology<Simplex>& boundary , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(boundary,_mesh)
{
  setDimension(boundary.vertices().dim());

  // the internal vertices remain empty, we simply will ask tetgen for a constrained triangulation
  boundaryFacets_.clear();
  boundaryFacets_.setSorted( boundary.sorted() );
  boundary.retrieveElements( 1 , boundaryFacets_ );

  // compute the expected volume as the integral of x.n over the boundary faces
  expectedVolume_ = fabs( boundary.calculateBoundingVolume() );

  // retrieve the boundary vertices and metadata
  boundary.vertices().copy( boundaryVertices_ );

}

Triangle::Triangle( ModelTessellation& tess , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(tess,_mesh)
{
  setDimension(2);

  // compute the expected volume as the integral of x.n over the boundary faces
  expectedVolume_ = tess.volume();

  // retrieve the boundary facets
  boundaryFacets_.setSorted( false ); // do not sort the triangles when we get them from the EGADS tessellation
  tess.retrieveElements( 1 , boundaryFacets_ );

  // retrieve the boundary vertices from the model tessellation
  tess.vertices().copy( boundaryVertices_ );

  // retrieve the internal points from the bodies
  for (index_t k=0;k<tess.nb_internal();k++)
    internalVertices_.create( tess.internalVertex(k) );
}


void triangulateioInit( struct triangulateio& io )
{

  io.pointlist = NULL;
  io.pointattributelist = NULL;
  io.pointmarkerlist = NULL;
  io.numberofpoints = 0;
  io.numberofpointattributes = 0;

  io.trianglelist = NULL;
  io.triangleattributelist = NULL;
  io.trianglearealist = NULL;
  io.neighborlist = NULL;
  io.numberoftriangles = 0;
  io.numberofcorners = 0;
  io.numberoftriangleattributes = 0;

  io.segmentlist = NULL;
  io.segmentmarkerlist = NULL;
  io.numberofsegments = 0;

  io.holelist = NULL;
  io.numberofholes = 0;

  io.regionlist = NULL;
  io.numberofregions = 0;

  io.edgelist = NULL;
  io.edgemarkerlist = NULL;
  io.normlist = NULL;
  io.numberofedges = 0;
}

void triangulateioDestroy( triangulateio& io )
{
  trifree( (void*)io.pointlist );
  trifree( (void*)io.pointattributelist );
  trifree( (void*)io.pointmarkerlist );

  trifree( (void*)io.trianglelist );
  trifree( (void*)io.triangleattributelist );
  trifree( (void*)io.trianglearealist );
  trifree( (void*)io.neighborlist );

  trifree( (void*)io.segmentlist );
  trifree( (void*)io.segmentmarkerlist );

  trifree( (void*)io.holelist );

  trifree( (void*)io.regionlist );

  trifree( (void*)io.edgelist );
  trifree( (void*)io.edgemarkerlist );
  trifree( (void*)io.normlist );

  // Set all pointers to NULL and counts to 0
  triangulateioInit( io );
}

Triangle::Triangle( Vertices& vertices , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(vertices,_mesh)
{
  setDimension(vertices.dim());

  boundaryFacets_.clear();

  expectedVolume_ = 0.;

  vertices.copy( boundaryVertices_ );
}

void
Triangle::call( const std::string& switches )
{

  // -p Triangulates a Planar Straight Line Graph
  // -e Outputs a list of edges of the triangulation.
  // -z Numbers all items starting from zero (rather than one). This switch is useful when calling Triangle from another program.
  // -n Outputs a list of triangles neighboring each triangle.
  // -Y Prohibits the insertion of Steiner points on the mesh boundary.
  // -V Verbose: Gives detailed information about what Triangle is doing. Add more `V's for increasing amount of detail.
  // -Q Quiet: Suppresses all explanation of what Triangle is doing, unless an error occurs.
  // -q Quality mesh generation with no angles smaller than 20 degrees. An alternate minimum angle may be specified after the `q'.

  struct triangulateio input,output;
  triangulateioInit(input);
  triangulateioInit(output);

  //boundaryVertices_.print("vbnd",true);

  // fill the vertices
  input.numberofpoints = boundaryVertices_.nb();
  input.pointlist = new REAL[input.numberofpoints*2];
  input.pointmarkerlist = new int[input.numberofpoints];
  for (index_t k=0;k<boundaryVertices_.nb();k++)
  {
    for (index_t j=0;j<2;j++)
      input.pointlist[2*k+j] = boundaryVertices_[k][j];
    if (entities_.size()==0)
    {
      input.pointmarkerlist[k] = boundaryVertices_.body(k);
    }
    else
    {
      Entity* e = boundaryVertices_.entity(k);
      if (e==NULL)
        input.pointmarkerlist[k] = -1;
      else
        input.pointmarkerlist[k] = label( boundaryVertices_.entity(k) ) +1;
    }
  }

  // fill the constraining facets
  input.numberofsegments = boundaryFacets_.nb();
  input.segmentlist = new int[ 2*input.numberofsegments ];
  for (index_t k=0;k<boundaryFacets_.nb();k++)
  {
    avro_assert(boundaryFacets_.nv(k)==2);

    input.segmentlist[2*k   ] = boundaryFacets_(k,0);
    input.segmentlist[2*k +1] = boundaryFacets_(k,1);
  }

  // fill the holes
  input.numberofholes = internalVertices_.nb();
  input.holelist = new REAL[input.numberofholes*2];
  for (index_t k=0;k<internalVertices_.nb();k++)
    for (index_t j=0;j<2;j++)
      input.holelist[2*k+j] = internalVertices_[k][j];

  printf("calling triangle with switches: %s\n",switches.c_str());
  triangulate( (char*)switches.c_str() , &input , &output , NULL );

  printf("number of tris = %d\n",output.numberoftriangles);
  if ( index_t(output.numberofpoints)>boundaryVertices_.nb() )
    printf("triangle added vertices!\n");

  // read the vertices provided by tetgen
  mesh_.vertices().clear();
  mesh_.vertices().setDimension(dim_);

  for (index_t k=0;k<index_t(output.numberofpoints);k++)
  {
    std::vector<real> coord(dim_,0.);
    for (index_t d=0;d<2;d++)
      coord[d] = output.pointlist[2*k+d];
    index_t id = mesh_.vertices().nb();
    mesh_.vertices().create( coord.data() );
    if (entities_.size()==0)
    {
      mesh_.vertices().setEntity(id,NULL);
      mesh_.vertices().body(id) = output.pointmarkerlist[k];
    }
    else
    {
      int marker = output.pointmarkerlist[k];
      if (marker > 0)
      {
        Entity* entity = entities_[marker-1];
        //entity->hierarchicalPrint();
        mesh_.vertices().setEntity(id,entity);
        mesh_.vertices().body(id) = entity->body()->index();
      }

    }
  }

  // read the tetrahedra, we only create one topology for the whole volume
  mesh_.clearTopologies();
  Topology_ptr t = smart_new(Topology<Simplex>)(mesh_.vertices(),2);
  for (index_t k=0;k<index_t(output.numberoftriangles);k++)
  {
    index_t tri[3];
    for (index_t j=0;j<3;j++)
      tri[j] = output.trianglelist[ k*3 +j ];
    t->add(tri,3);
  }
  mesh_.addTopology(t);

  // get the actual volume
  actualVolume_ = fabs( t->calculateVolume() );
  printf("expected volume = %.12e, actual = %.12e\n",expectedVolume_,actualVolume_);

  // check the facets are the same
  // TODO

  //triangulateioDestroy(input);
  //triangulateioDestroy(output);

}

TetGen::TetGen( Model* _model , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(_model,_mesh)
{
  setDimension(3);

  // tessellate the model
  real param[3] = {1.8,1,30.};
  ModelTessellation tess(*model_,param);

  // compute the expected volume as the integral of x.n over the boundary faces
  expectedVolume_ = tess.volume();

  // retrieve the boundary facets
  // do not sort the triangles when we get them from the EGADS tessellation
  boundaryFacets_.setSorted( false );
  tess.retrieveElements( 2 , boundaryFacets_ );

  // retrieve the boundary vertices from the model tessellation
  tess.vertices().copy( boundaryVertices_ );

  // retrieve the internal points from the bodies
  for (index_t k=0;k<tess.nb_internal();k++)
    internalVertices_.create( tess.internalVertex(k) );
}

TetGen::TetGen( Topology<Simplex>& boundary , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(boundary,_mesh)
{
  setDimension(3);

  // the internal vertices remain empty, we simply will ask tetgen for a constrained triangulation
  boundaryFacets_.clear();
  boundaryFacets_.setSorted( boundary.sorted() );
  boundary.retrieveElements( 2 , boundaryFacets_ );

  // compute the expected volume as the integral of x.n over the boundary faces
  expectedVolume_ = fabs( boundary.calculateBoundingVolume() );

  boundary.vertices().copy( boundaryVertices_ );

}

TetGen::TetGen( Vertices& vertices , Mesh<Simplex>& _mesh ) :
  Mesher<Simplex>(vertices,_mesh)
{
  setDimension(3);

  boundaryFacets_.clear();

  expectedVolume_ = 0.;

  vertices.copy( boundaryVertices_ );
}

void
TetGen::call( const std::string& switches )
{
  // -p Tetrahedralizes a piecewise linear complex (PLC).
  // -Y Preserves the input surface mesh (does not modify it).
  // -V verbose
  // -q mesh quality (maximum radius-edge ratio)/(minimum dihedral angle)
  // -a maximum volume constraint
  // -f provides the interior+boundry triangular faces
  // -nn to get tet neighbors for each triangular face
  // -k dumps to paraview when last argument is NULL
  // -m Applies a mesh sizing function
  // -C Checks the consistency of the final mesh.

  tetgenio input,output;
  tetgenio::facet *f;
  tetgenio::polygon *p;

  // fill the vertices
  input.numberofpoints = boundaryVertices_.nb();
  input.pointlist = new REAL[input.numberofpoints*dim_];
  input.pointmarkerlist = new int[input.numberofpoints];
  for (index_t k=0;k<boundaryVertices_.nb();k++)
  {
    for (index_t j=0;j<dim_;j++)
      input.pointlist[dim_*k+j] = boundaryVertices_[k][j];

    if (entities_.size()==0)
    {
      input.pointmarkerlist[k] = boundaryVertices_.body(k);
    }
    else
    {
      Entity* e = boundaryVertices_.entity(k);
      if (e==NULL)
        input.pointmarkerlist[k] = -1;
      else
        input.pointmarkerlist[k] = label( boundaryVertices_.entity(k) ) +1;
    }
  }

  // fill the constraining facets
  input.numberoffacets = boundaryFacets_.nb();
  input.facetlist = new tetgenio::facet[ input.numberoffacets ];
  for (index_t k=0;k<boundaryFacets_.nb();k++)
  {
    avro_assert(boundaryFacets_.nv(k)==3);

    f = &input.facetlist[k];

    f->numberofpolygons = 1;
    f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
    f->numberofholes = 0;
    f->holelist = NULL;

    p = &f->polygonlist[0];
    p->numberofvertices = boundaryFacets_.nv(k);
    p->vertexlist = new int[p->numberofvertices];
    for (index_t j=0;j<boundaryFacets_.nv(k);j++)
      p->vertexlist[j] = boundaryFacets_[k][j];
  }

  // fill the holes
  input.numberofholes = internalVertices_.nb();
  input.holelist = new REAL[input.numberofholes*dim_];
  for (index_t k=0;k<internalVertices_.nb();k++)
    for (index_t j=0;j<dim_;j++)
      input.holelist[dim_*k+j] = internalVertices_[k][j];

  tetrahedralize( (char*)switches.c_str() , &input , &output );

  printf("number of tets = %d\n",output.numberoftetrahedra);
  if ( index_t(output.numberofpoints)>boundaryVertices_.nb() )
    printf("tetgen added vertices!\n");

  // read the vertices provided by tetgen
  mesh_.vertices().clear();
  mesh_.vertices().setDimension(dim_);

  for (index_t k=0;k<index_t(output.numberofpoints);k++)
  {
    index_t id = mesh_.vertices().nb();
    mesh_.vertices().create( &output.pointlist[dim_*k] );

    if (entities_.size()==0)
    {
      mesh_.vertices().setEntity(id,NULL);
      mesh_.vertices().body(id) = output.pointmarkerlist[k];
    }
    else
    {
      int marker = output.pointmarkerlist[k];
      if (marker > 0)
      {
        Entity* entity = entities_[marker-1];
        mesh_.vertices().setEntity(id,entity);
        mesh_.vertices().body(id) = entity->body()->index();
      }
    }
  }

  // read the tetrahedra, we only create one topology for the whole volume
  mesh_.clearTopologies();
  Topology_ptr t = smart_new(Topology<Simplex>)(mesh_.vertices(),3);
  for (index_t k=0;k<index_t(output.numberoftetrahedra);k++)
  {
    index_t tet[4];
    for (index_t j=0;j<4;j++)
      tet[j] = output.tetrahedronlist[ k*4 +j ];
    t->add(tet,4);
  }
  mesh_.addTopology(t);

  // get the actual volume
  actualVolume_ = fabs( t->calculateVolume() );
  printf("expected volume = %.12e, actual = %.12e\n",expectedVolume_,actualVolume_);
}

} // avro
