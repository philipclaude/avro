#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/parameters.h"
#include "adaptation/properties.h"

#include "library/meshb.h"

#include "mesh/boundary.h"
#include "mesh/mesh.h"
#include "mesh/topology.h"

#include "numerics/predicates.h"

typedef luna::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <unistd.h>
#include <fstream>

namespace luna
{

template<typename type>
AdaptThread<type>::AdaptThread( Topology<type>& topology , MetricField<type>& metric , AdaptationParameters& params ) :
  topology_(topology),
  metric_(metric),
  params_(params),
  collapser_(topology),
  inserter_(topology),
  smoother_(topology),
  edge_swapper_(topology)
{}

const real_t nb_smooth = 10;
const real_t qt_max = 0.8;
const bool limit_collapse_lmax = false;

template<typename type>
int
call( Topology<type>& topology , Topology<type>& mesh_topology ,
      MetricField<type>& metric , AdaptationParameters& params ,
      Mesh& mesh_out )
{
  const coord_t number = topology.number();

  // retrieve the parameters
  const bool limit_insertion_length = params.limit_insertion_length();
  const bool swapout = false;//params.swapout();
  const real_t lt_min = params.lt_min();
  real_t lt_max = params.lt_max();
  const bool smooth_on = params.use_smoothing();
  const bool fefloa = params.fefloa();

  if (fefloa) lt_max = sqrt(2.0); // and we will do one pass

#if 0
  // project the mesh_topology points to the geometry
  for (index_t k=0;k<mesh_topology.points().nb();k++)
  {
    if (mesh_topology.points().entity(k)==NULL) continue;
    Entity* e = mesh_topology.points().entity(k);
    real_t* x = mesh_topology.points()[k];
    std::vector<real_t> X(x,x+mesh_topology.points().dim());
    e->project(X);
    for (coord_t d=0;d<mesh_topology.points().dim();d++)
      x[d] = X[d];
  }

  // project the metric topology points to the geometry
  for (index_t k=0;k<topology.points().nb();k++)
  {
    if (topology.points().entity(k)==NULL) continue;
    Entity* e = topology.points().entity(k);
    real_t* x = topology.points()[k];
    std::vector<real_t> X(x,x+topology.points().dim());
    e->project(X);
    for (coord_t d=0;d<topology.points().dim();d++)
      x[d] = X[d];
  }
#endif

  // the uv-parameters might not be set for the incoming points
  // for the case of real_t geometries
  if (params.curved())
  {
    // check all parametric coordinates for consistency
    coord_t udim = mesh_topology.points().udim();
    std::vector<real_t> U(udim,mesh_topology.points().INFTY);
    for (index_t k=0;k<mesh_topology.points().nb();k++)
    {
      Entity* e = mesh_topology.points().entity(k);
      if (e==NULL) continue;

      // project the mesh_topology points to the geometry if not provided
      if(!params.has_uv())
      {
        real_t* x = mesh_topology.points()[k];
        std::vector<real_t> X(x,x+mesh_topology.points().dim());
        //e->project(X,U);
        e->inverse(X,U);
        for (coord_t d=0;d<mesh_topology.points().dim();d++)
          x[d] = X[d];

        mesh_topology.points().set_param( k , U );
      }
      else
      {
        for (coord_t d=0;d<udim;d++)
          U[d] = mesh_topology.points().u( k )[d];
      }

      // evaluate the coordinates for the parameters we found
      std::vector<real_t> x_eval(3);
      e->evaluate(U, x_eval);

      // check the coordinates are close
      const real_t *x0 = mesh_topology.points()[k];
      real_t d = numerics::distance2(x0,x_eval.data(),mesh_topology.points().dim() );
      real_t tol=0;
      EGADS::Object* e0 = (EGADS::Object*)e;
      EGADS_ENSURE_SUCCESS( EG_getTolerance( *e0->object(), &tol ) );

      if (d>tol)
      {
        real_t u[2];
        u[0] = U[0];
        u[1] = udim == 2 ? U[1] : mesh_topology.points().INFTY;
        printf("Vertex %lu is not on geometry!\n",k);
        printf("Provided x0 = (%g,%g,%g), x_eval = (%g,%g,%g) with param coordinates (%g,%g)\n",
                  x0[0],x0[1],x0[2],x_eval[0],x_eval[1],x_eval[2],u[0],u[1]);
        printf("For entity:\n");
        e->print(false);
      }

      luna_assert_msg( d <= tol , "d = %1.16e, tol = %1.16e" , d , tol );
    }
  }

  // setup the adaptation and do some checks
  AdaptThread<type> adaptation( mesh_topology , metric , params );
  adaptation.check("initial");

  // initialize the metric-conformity properties
  Properties properties( mesh_topology , metric );

  // collapse short edges in the mesh without swapping upon rejection
  properties.print( "pre-collapse" );
  adaptation.collapse_edges(false,swapout); // no lmax limit, swapout
  adaptation.check("collapses");

  // split long edges without swapping on rejection
  properties.compute( mesh_topology , metric );
  properties.print("pre-insertion" );
  adaptation.split_edges(lt_max,limit_insertion_length,swapout); // blind=false, swapout
  adaptation.check("insertions");

  // swap edges to improve quality
  properties.compute( mesh_topology , metric );
  properties.print("pre-swap" );
  adaptation.swap_edges(0.4,10); // no length limit
  adaptation.swap_edges(qt_max,4); // no length limit
  adaptation.check("swaps");

  // perform some vertex smoothing to drive lengths to 1
  if (smooth_on)
  {
    properties.compute( mesh_topology , metric );
    properties.print("pre-smooth" );
    adaptation.smooth_points(nb_smooth);
    adaptation.check("smoothing");

    // swap edges to improve quality
    properties.compute( mesh_topology , metric );
    properties.print("pre-swap" );
    adaptation.swap_edges(0.4,10,true); // limit length
    adaptation.swap_edges(qt_max,4,true); // limit length
    adaptation.check("swaps");
  }

  if (fefloa) goto done;

  // revisit collapses and splits and swap upon rejection
  properties.compute( mesh_topology , metric );
  properties.print("pre-collapse/split" );
  adaptation.collapse_edges(limit_collapse_lmax,swapout); // limt length, swapout
  adaptation.check("collapses");
  adaptation.split_edges(lt_max,limit_insertion_length,swapout); // swapout
  adaptation.check("insertions");

  // swap facets then edges to improve quality
  properties.compute( mesh_topology , metric );
  properties.print("pre-swap" );
  adaptation.swap_edges(0.4,10); // no length limit
  adaptation.swap_edges(qt_max,2); // no length limit
  adaptation.check("swaps");

  // perform some vertex smoothing to drive lengths to 1
  if (smooth_on)
  {
    properties.compute( mesh_topology , metric );
    properties.print("pre-smooth" );
    adaptation.smooth_points(nb_smooth);
    adaptation.check("smoothing");
  }

  // revisit collapses and splits and swap upon rejection
  properties.compute( mesh_topology , metric );
  properties.print("pre-collapse/split" );
  adaptation.collapse_edges(limit_collapse_lmax,swapout); // lmax limit, swapout
  adaptation.check("collapses");
  adaptation.split_edges(lt_min,limit_insertion_length,swapout); // swapout
  adaptation.check("insertions");

  // swap edges to improve quality
  properties.compute( mesh_topology , metric );
  properties.print("pre-swap" );
  adaptation.swap_edges(0.4,10,true); // limit length
  adaptation.swap_edges(qt_max,4,true); // limit length
  adaptation.check("swaps");

  // perform some vertex smoothing to drive lengths to 1
  if (smooth_on)
  {
    properties.compute( mesh_topology , metric );
    properties.print("pre-smooth" );
    adaptation.smooth_points(nb_smooth);
    adaptation.check("smoothing");
  }

  // revisit collapses and splits and swap upon rejection
  properties.compute( mesh_topology , metric );
  properties.print("pre-collapse/split" );
  adaptation.collapse_edges(limit_collapse_lmax,swapout); // lmax limit, swapout
  adaptation.check("collapses");
  adaptation.split_edges(lt_min,limit_insertion_length,swapout); // blind=false and swapout
  adaptation.check("insertions");

  // perform some vertex smoothing to drive lengths to 1
  if (smooth_on)
  {
    properties.compute( mesh_topology , metric );
    properties.print("pre-smooth" );
    adaptation.smooth_points(nb_smooth);
    adaptation.check("smoothing");
  }

  // swap edges to improve quality
  properties.compute( mesh_topology , metric );
  properties.print("pre-swap" );
  adaptation.swap_edges(0.4,10,true); // limit length
  adaptation.swap_edges(qt_max,4,true); // limit length
  adaptation.check("swaps");

done:

  // compute the final properties and dump to file
  properties.compute( mesh_topology , metric );
  properties.print("final metric conformity" );
  if (params.write_conformity())
    properties.dump( params.directory()+
                     "/properties_"+stringify(params.adapt_iter())+".json");

 // copy back into the original mesh
 topology.TopologyBase::copy( mesh_topology );
 mesh_topology.points().copy( topology.points() );

 // create the output mesh
 if (mesh_out.nb_topologies()==0)
 {
   std::shared_ptr<Topology<type>> ptopology_out =
                std::make_shared<Topology<type>>( mesh_out.points() ,
                                                number );
   mesh_out.add(ptopology_out);
 }

 Topology<type>& topology_out = mesh_out.retrieve<type>(0);
 mesh_topology.points().copy( mesh_out.points() );
 topology_out.TopologyBase::copy( mesh_topology );
 //mesh_topology.neighbours().copy( topology_out.neighbours() );
 //topology_out.inverse().copy( mesh_topology.inverse() );

 // identify ghost elements
 std::vector<index_t> ghost;
 for (index_t k=0;k<topology.nb();k++)
 {
   if (topology.ghost(k))
     ghost.push_back(k);
 }

 // remove ghost elements
 for (index_t k=0;k<ghost.size();k++)
 {
   topology.remove( ghost[k]-k );
   topology_out.remove( ghost[k]-k );
 }

 // remove ghost points
 index_t nb_ghost = mesh_topology.points().nb_ghost();
 for (index_t k=0;k<nb_ghost;k++)
 {
   topology.remove_point( k );
   topology_out.remove_point(k);
   metric.attachment().remove(k,false); // don't check
 }

 topology.set_closed(false);
 topology_out.set_closed(false);

 luna_assert( topology_out.points().nb_ghost()==0 );

 return 0;
}

template<typename type>
int
adapt( AdaptationProblem& problem )
{
  //Gamma<Simplex> gamma; // used to write meshes if necessary

  // standardize the parameters
  AdaptationParameters& params = problem.params;

  const std::string& output_redirect = params.output_redirect();

  int redirected_fd = 0;
  FILE *redirected_fid = 0;
  if (!output_redirect.empty())
  {
    fflush(stdout);
    redirected_fd = dup(fileno(stdout));
    redirected_fid = freopen(output_redirect.c_str(),"w",stdout);
  }

  params.standard();
  params.print();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();

  // retrieve the background mesh and metric field
  Mesh& mesh = problem.mesh_in;
  //VertexField<numerics::SPDT<real_t>>& fld = problem.fld;
  //MetricAttachment& fld = problem.fld;
  std::vector<numerics::SymMatrixD<real_t>>& fld = problem.fld;

  // extract the background topology
  luna_assert_msg( mesh.nb_topologies()==1 ,
                  "nb_topologies = %lu",mesh.nb_topologies() );
  Topology<type>& topology = mesh.retrieve<type>(0);
  if (topology.nb()==0)
  {
    printf("there are no elements?");
    return 0;
  }

  // retrieve some info about the mesh
  const coord_t number = mesh.number();
  luna_assert( number == topology.number() );

  const coord_t dim = mesh.points().dim();
  luna_assert( number == dim );

  // create a mesh topology by copying the input one
  Points points( dim );
  mesh.points().copy( points );
  Topology<type> mesh_topology( points , number );
  mesh_topology.TopologyBase::copy( topology );

  if (!params.prepared())
  {
    // check there are no ghosts
    luna_assert_msg( mesh.points().nb_ghost()==0 , "there are %lu ghosts!" , mesh.points().nb_ghost() );

    // check the number of tensors equals the number of points
    luna_assert_msg( mesh.points().nb() == fld.size() ,
                      "nb_points = %lu, nb_fld = %lu",
                      mesh.points().nb() , fld.size() );

    // close the mesh topology, compute the neighours and mesh inverse
    mesh_topology.close();
    mesh_topology.orient();

    mesh_topology.neighbours().fromscratch() = true; // speed up
    mesh_topology.neighbours().compute();
    mesh_topology.neighbours().fromscratch() = false;
    mesh_topology.inverse().build();

    // copy the data into the background topology used by the discrete metric
    points.copy( topology.points() );
    topology.TopologyBase::copy( mesh_topology );
    mesh_topology.neighbours().copy( topology.neighbours() );

    topology.inverse().build();

  }
  else
  {
    luna_assert_not_reached;

    // use the stored mesh
    // check the number of tensors equals the number of points
    luna_assert_msg( mesh.points().nb()-mesh.points().nb_ghost() == fld.size() ,
      "nb_points = %lu, nb_ghost = %lu, fld.nb = %lu" ,
      mesh.points().nb(), mesh.points().nb_ghost(),fld.size() );

    // both the topology and mesh_topology should be properly oriented
    // we just need to copy in the neighbours into the mesh_topology
    topology.neighbours().copy( mesh_topology.neighbours() );
    mesh_topology.inverse().copy( topology.inverse() );
  }

  // extract the boundaries to check the vertex/geometry association
  #if 0
  Boundary<Simplex> bnd( mesh_topology );
  bnd.extract( params.has_interior_boundaries() );

  Topology<Simplex> bnd_topology( bnd.points() , mesh_topology.number()-1 );
  bnd.retrieve( bnd_topology );

  // check all points on geometries are accounted for in the boundary topology
  Points& verts = mesh_topology.points();
  std::vector<bool> accounted( verts.nb() , false );
  for (index_t k=0;k<accounted.size();k++)
  {
    if (verts.entity(k)==NULL)
      accounted[k] = true;
  }

  for (index_t k=0;k<bnd_topology.nb();k++)
  {
    for (index_t j=0;j<bnd_topology.nv(k);j++)
      accounted[ bnd_topology(k,j) ] = true;
  }

  // there should be no points leftover that were not accounted for
  index_t nerr = 0;
  for (index_t k=0;k<accounted.size();k++)
  {
    if (verts.entity(k)!=NULL && !accounted[k])
    {
      printf("error! vertex %lu is tagged on geometry!\n",k);
      verts.print(k,true);
      nerr++;
    }
  }
  luna_assert_msg(nerr==0,
    "there are interior points tagged on the geometry! if you have interior boundary groups, set this flag to true!");
  #endif

  // check how many ghost points are present and adjust the field
  if (fld.size() < mesh_topology.points().nb() )
  {
    luna_assert( fld.size() == mesh_topology.points().nb() - mesh_topology.points().nb_ghost() );
    index_t nb_ghost = mesh_topology.points().nb_ghost();
    for (index_t k=0;k<nb_ghost;k++)
    {
      Metric tensor(number);
      for (coord_t d=0;d<number;d++)
        tensor(d,d) = 1.;
      fld.push_back( tensor );
    }
  }

  // set the initial cells (for searching) in the field
  // the field holds the adapted set of points but those points point
  // to the background topology
  MetricAttachment field( mesh_topology.points() , fld );
  field.set_cells( topology );
  luna_assert( field.check () );

  // create a discrete metric with the input topology
  MetricField<type> metric( topology , field );

  // initial volume to compare with later
  real_t v0 = mesh_topology.volume();

  // call the adaptation!
  int result = call( topology , mesh_topology , metric , params , problem.mesh_out );

  std::string mesh_file = params.directory()+params.prefix()+"_"+stringify(params.adapt_iter())+".mesh";

  if (params.write_mesh())
  {
    #if 0
    if (topology.number()<=3 && !params.has_interior_boundaries() && params.write_meshb())
    {
      // write the full mesh
      gamma.writeMesh( problem.mesh_out , mesh_file );
      printf("wrote mesh %s\n",mesh_file.c_str());
    }
    else if (topology.number()==4 && params.write_meshb())
    {
      // get the boundary of the mesh
      Boundary<Simplex> tesseract_boundary( mesh_topology );
      tesseract_boundary.extract();

      for (index_t k=0;k<tesseract_boundary.nb_children();k++)
      {
        if (tesseract_boundary.child(k)->number()==3)
        {
          // write out the mesh
          library::Plottable<Simplex> plot( *tesseract_boundary.child(k) );
          gamma.writeMesh( plot , params.directory()+"/"+params.boundarySubdirectory()
                                  +"/"+params.prefix()+"_bnd"+stringify(k)+"_a"+stringify(params.adapt_iter())+".mesh" );
        }
      }
    }
    #else
    //luna_implement;
    #endif

    if (topology.number()==3 && params.has_interior_boundaries())
    {
      #if 0
      Boundary<Simplex> bnd2( mesh_topology );
      bnd2.extract(true);

      Topology<Simplex> bnd_topo( mesh_topology.points() , 2 );
      bnd2.retrieve( bnd_topo );

      library::Plottable<Simplex> plot(bnd_topo);
      gamma.withGeometry() = false;
      std::string boundary_file = params.directory()+params.prefix()+"_boundary_"+stringify(params.adapt_iter())+".mesh";

      gamma.writeMesh( plot , boundary_file , false );
      printf("wrote boundary of mesh %s\n",boundary_file.c_str());

      gamma.writeMesh( problem.mesh_out , mesh_file , false );
      printf("wrote mesh %s\n",mesh_file.c_str());
      #else
      //luna_implement;
      #endif
    }

    if (!params.has_interior_boundaries() && params.write_json())
    {
      #if 0
      const std::string json_file = params.directory()+"mesh_"+stringify(params.adapt_iter())+".json";
      io::writeMesh<Simplex>( json_file , topology , &metric.field() );
      #else
      //luna_implement;
      #endif
    }

  } // option to write the mesh

  // save the metric in case the caller wants information about it
  luna_assert_msg( field.nb() == topology.points().nb() , "|field| = %lu, nv = %lu", field.nb() , topology.points().nb() );
  fld.clear();
  for (index_t k=0;k<field.nb();k++)
    fld.push_back( field[k] );

  // print some information about the volume (for straight-sided geometries)
  real_t v1 = mesh_topology.volume();
  printf("v0 = %1.16e, v1 = %1.16e\n",v0,v1);

  if (!output_redirect.empty())
  {
    fflush(stdout);
    dup2(redirected_fd,fileno(stdout));
    close(redirected_fd);
    fclose(redirected_fid);
  }

  return result;
}

template class AdaptThread<Simplex>;
template int adapt<Simplex>( AdaptationProblem& );

} // luna
