//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/adapt.h"
#include "adaptation/metric.h"
#include "adaptation/properties.h"

#include "geometry/entity.h"
#include "geometry/egads/object.h"

#include "library/meshb.h"
#include "library/spacetime.h"

#include "mesh/boundary.h"
#include "mesh/mesh.h"
#include "mesh/topology.h"

#include "numerics/predicates.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <json/json.hpp>

#include <unistd.h>
#include <fstream>

namespace avro
{

AdaptationParameters::AdaptationParameters() {
  set_defaults();
}

AdaptationParameters::AdaptationParameters( const ParameterSet& params ) :
  ParameterSet(params) {
  set_defaults();
}

void
AdaptationParameters::set_defaults() {
  register_parameter( "swapout" , false , "whether to swap out of restrictive configurations when insertions/collapses are rejected" );
  register_parameter( "insertion volume factor" , -1.0 , "whether to limit insertions by checking complexity after insertion" );
  register_parameter( "use smoothing" , true , "whether to use smoothing during the adaptation" );
  register_parameter( "allow serial" , false , "whether to allow serial adaptation when running in parallel" );
  register_parameter( "has uv" , false , "whether parameter space coordinates are specified for geometry points in the mesh" );
  register_parameter( "smoothing exponent" , index_t(1) , "smoothing exponent" );
  register_parameter( "elems per processor" , index_t(10000) , "target number of elements for each processor" );
  register_parameter( "adapt iter" , index_t(1) , "adaptation iteration" );
  register_parameter( "has interior boundaries" , false , "whether there are geometry entities which are embedded interior to the mesh" );
  register_parameter( "export boundary" , false , "should the boundary be included in the exported mesh" );
  register_parameter( "write conformity" , false , "whether metric conformity information should be written when adaptation finishes" );
  register_parameter( "partitioned" , false , "whether the incoming mesh is already partitioned" );
  register_parameter( "force partition count" , index_t(0) , "whether the number of partitions should be forced during every pass" );
}

template<typename type>
AdaptThread<type>::AdaptThread( Topology<type>& topology , MetricField<type>& metric , AdaptationParameters& params ) :
  topology_(topology),
  metric_(metric),
  params_(params),
  collapser_(topology),
  inserter_(topology),
  smoother_(topology),
  edge_swapper_(topology)
{
  bool curved = params["curved"];
  collapser_.curved() = curved;
  inserter_.curved() = curved;
  smoother_.curved() = curved;
  edge_swapper_.curved() = curved;

  index_t smoothing_exponent = params["smoothing exponent"];
  smoother_.exponent() = smoothing_exponent;
}

const real_t nb_smooth = 10;
const real_t qt_max = 0.8;
const bool limit_collapse_lmax = false;
const index_t nb_swap_pass_qt_min = 5;
const index_t nb_swap_pass_qt_max = 2;

template<typename type>
int
call( Topology<type>& topology , Topology<type>& mesh_topology ,
      MetricField<type>& metric , AdaptationParameters& params ,
      Mesh& mesh_out )
{
  const coord_t number = topology.number();

  // retrieve the parameters
  const bool limit_insertion_length = true;
  const bool swapout = params["swapout"];
  const real_t lt_min = sqrt(2.0);
  real_t lt_max = 2.0;
  const bool smooth_on = params["use smoothing"];
  const bool fefloa = false;

  if (fefloa) lt_max = sqrt(2.0); // and we will do one pass

  // the uv-parameters might not be set for the incoming points
  // for the case of real_t geometries
  bool curved = params["curved"];
  if (curved)
  {
    // check all parametric coordinates for consistency
    coord_t udim = mesh_topology.points().udim();
    std::vector<real_t> U(udim,mesh_topology.points().INFTY);
    for (index_t k=0;k<mesh_topology.points().nb();k++)
    {
      Entity* e = mesh_topology.points().entity(k);
      if (e==NULL) continue;

      // project the mesh_topology points to the geometry if not provided
      bool has_uv = params["has uv"];
      if(!has_uv)
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

      // skip this check if we are operating purely in parameter space
      if (!mesh_topology.element().parameter())
      {
        // evaluate the coordinates for the parameters we found
        std::vector<real_t> x_eval(3);
        e->evaluate(U, x_eval);

        // check the coordinates are close
        const real_t *x0 = mesh_topology.points()[k];
        real_t d = numerics::distance2(x0,x_eval.data(),mesh_topology.points().dim() );
        real_t tol=0;
        EGADS::Object* e0 = (EGADS::Object*)e;
        EGADS_ENSURE_SUCCESS( EG_getTolerance( e0->object(), &tol ) );

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

        avro_assert_msg( d <= tol , "d = %1.16e, tol = %1.16e" , d , tol );
      }
    }
  }

  // setup the adaptation and do some checks
  AdaptThread<type> adaptation( mesh_topology , metric , params );
  adaptation.check("initial");

  // initialize the metric-conformity properties
  Properties properties( mesh_topology , metric );

  if (topology.nb() == 0) goto done;

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
  adaptation.swap_edges(0.4,nb_swap_pass_qt_min); // no length limit
  adaptation.swap_edges(qt_max,nb_swap_pass_qt_max); // no length limit
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
    adaptation.swap_edges(0.4,nb_swap_pass_qt_min,true); // limit length
    adaptation.swap_edges(qt_max,nb_swap_pass_qt_max,true); // limit length
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
  adaptation.swap_edges(0.4,nb_swap_pass_qt_min); // no length limit
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
  adaptation.swap_edges(0.4,nb_swap_pass_qt_min,true); // limit length
  adaptation.swap_edges(qt_max,nb_swap_pass_qt_max,true); // limit length
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
  adaptation.swap_edges(0.4,nb_swap_pass_qt_min,true); // limit length
  adaptation.swap_edges(qt_max,nb_swap_pass_qt_max,true); // limit length
  adaptation.check("swaps");

done:

  // compute the final properties and dump to file
  if ( mesh_topology.nb() > 0 )
  {
    properties.compute( mesh_topology , metric );
    properties.print("final metric conformity" );
    std::string directory = params["directory"];
    std::string prefix = params["prefix"];
    index_t adapt_iter = params["adapt iter"];
    if (params["write conformity"])
      properties.dump( prefix + "_properties_" + stringify(adapt_iter) + ".json");
  }
  else avro_implement;

 // copy back into the original mesh
 topology.TopologyBase::copy( mesh_topology );
 mesh_topology.points().copy( topology.points() );

 // create the output mesh
 if (mesh_out.nb_topologies()==0)
 {
   std::shared_ptr<Topology<type>> ptopology_out =
                std::make_shared<Topology<type>>( mesh_out.points() , number );
   mesh_out.add(ptopology_out);
 }

 Topology<type>& topology_out = mesh_out.retrieve<type>(0);
 mesh_topology.points().copy( mesh_out.points() );
 topology_out.TopologyBase::copy( mesh_topology );

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

 avro_assert( topology_out.points().nb_ghost()==0 );

 return 0;
}

template<typename type>
int
adapt( AdaptationProblem& problem )
{
  // standardize the parameters
  AdaptationParameters& params = problem.params;

  std::string output_redirect = params["output redirect"];

  fpos_t pos;
  int redirected_fd = 0;
  FILE *redirected_fid;
  if (!output_redirect.empty())
  {
    fgetpos(stdout,&pos);
    fflush(stdout);
    redirected_fd = dup(fileno(stdout));
    redirected_fid = freopen(output_redirect.c_str(),"w",stdout);
  }
  UNUSED( redirected_fid );

  //params.standard();
  params.print();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();

  // retrieve the background mesh and metric field
  Mesh& mesh = problem.mesh_in;
  std::vector<symd<real_t>>& fld = problem.fld;

  #if 0
  for (index_t k = 0; k < fld.size(); k++) {
    const symd<real_t>& T = fld[k];
    std::pair< vecd<real_t> , matd<real_t> > decomp = numerics::eig(T);
    real_t hmin = 1e-8;
    for (index_t d = 0; d < T.n(); d++) {
      real_t h = 1./sqrt(decomp.first(d));
      avro_assert( h > hmin );
    }
  }
  #endif

  // extract the background topology
  avro_assert_msg( mesh.nb_topologies()==1 ,
                  "nb_topologies = %lu",mesh.nb_topologies() );
  Topology<type>& topology = mesh.retrieve<type>(0);
  if (topology.nb()==0)
  {
    printf("there are no elements?\n");
    std::shared_ptr<Topology<type>> ptopology_out =
                 std::make_shared<Topology<type>>( problem.mesh_out.points() ,
                                                   topology.number() );
    problem.mesh_out.add(ptopology_out);
    Topology<type>& topology_out = problem.mesh_out.retrieve<type>(0);
    topology.points().copy( problem.mesh_out.points() );
    topology_out.TopologyBase::copy( topology );
    return 0;
  }

  for (index_t k=0;k<topology.nb();k++)
  {
    for (index_t j=0;j<topology.nv(k);j++)
      avro_assert( topology(k,j) < topology.points().nb() );
  }

  // retrieve some info about the mesh
  const coord_t number = mesh.number();
  avro_assert( number == topology.number() );
  const coord_t dim = mesh.points().dim();

  // create a mesh topology by copying the input one
  Points points( dim );
  mesh.points().copy( points );
  Topology<type> mesh_topology( points , number );
  mesh_topology.TopologyBase::copy( topology );
  mesh_topology.element().set_parameter( topology.element().parameter() );

  // check there are no ghosts
  avro_assert_msg( mesh.points().nb_ghost()==0 , "there are %lu ghosts!" , mesh.points().nb_ghost() );

  // check the number of tensors equals the number of points
  avro_assert_msg( mesh.points().nb() == fld.size() ,
                    "nb_points = %lu, nb_fld = %lu",
                    mesh.points().nb() , fld.size() );

  // close the mesh topology, compute the neighours and mesh inverse
  mesh_topology.close();

  for (index_t k=0;k<mesh_topology.nb();k++)
  {
    for (index_t j=0;j<mesh_topology.nv(k);j++)
      avro_assert_msg( mesh_topology(k,j) < mesh_topology.points().nb() , "topology(%lu,%lu) = %lu, but |points| = %lu" , k,j,mesh_topology(k,j),mesh_topology.points().nb() );

  }

  if (dim==number)
    mesh_topology.orient();

  mesh_topology.neighbours().fromscratch() = true; // speed up
  mesh_topology.neighbours().compute();
  mesh_topology.neighbours().fromscratch() = false;
  mesh_topology.inverse().build();

  // copy the data into the background topology used by the discrete metric
  points.copy( topology.points() );
  topology.TopologyBase::copy( mesh_topology );
  topology.element().set_parameter( mesh_topology.element().parameter() );
  mesh_topology.neighbours().copy( topology.neighbours() );

  topology.inverse().build();

  // check how many ghost points are present and adjust the field
  if (fld.size() < mesh_topology.points().nb() )
  {
    avro_assert( fld.size() == mesh_topology.points().nb() - mesh_topology.points().nb_ghost() );
    index_t nb_ghost = mesh_topology.points().nb_ghost();
    for (index_t k=0;k<nb_ghost;k++)
    {
      Metric tensor(number);
      for (coord_t d=0;d<number;d++)
        tensor(d,d) = 1.;
      fld.insert( fld.begin() , tensor );
    }
  }

  // set the initial cells (for searching) in the field
  // the field holds the adapted set of points but those points point
  // to the background topology
  MetricAttachment field( mesh_topology.points() , fld );
  field.set_cells( topology );
  avro_assert( field.check () );

  bool limit_metric = params["limit metric"];
  if (limit_metric)
  {
    // option to limit the metric field from the current mesh-implied metric
    // any problematic vertices will be fixed
    field.limit(mesh_topology,2.0);
  }

  // create a discrete metric with the input topology
  // and pass in the interpolator (null if not provided)
  MetricField<type> metric( topology , field );

  std::shared_ptr< FieldInterpolation<type,Metric> > interpolation = nullptr;
  if (problem.interpolation==nullptr)
  {
		interpolation = std::make_shared< FieldInterpolation<type,Metric> >(&metric);
    metric.set_interpolation(interpolation.get());
  }
	else
		metric.set_interpolation( problem.interpolation );

  topology.element().set_basis( BasisFunctionCategory_Lagrange );
  mesh_topology.element().set_basis( BasisFunctionCategory_Lagrange );

  // call the adaptation!
  int result = call( topology , mesh_topology , metric , params , problem.mesh_out );

  std::string directory = params["directory"];
  std::string prefix = params["prefix"];
  index_t adapt_iter = params["adapt iter"];
  bool write_mesh = params["write mesh"];
  bool has_interior_boundaries = params["has interior boundaries"];
  bool export_boundary = params["export boundary"];

  // option to output the mesh
  std::string mesh_file = directory + prefix + "_" + stringify(adapt_iter) + ".mesh";
  if (write_mesh)
  {
    if (topology.number()<=3 && !has_interior_boundaries)
    {
      // write the full mesh
      library::meshb meshb;
      meshb.write( problem.mesh_out , mesh_file , export_boundary );
      printf("wrote mesh %s\n",mesh_file.c_str());
    }
    else if (topology.number()==4 && export_boundary)
    {
      // get the boundary of the mesh
      Topology_Spacetime<Simplex> spacetime(mesh_topology);
      spacetime.extract();

      std::string filename = directory + "/" + prefix + "_" + stringify(adapt_iter) + ".mesh";
      spacetime.write( filename );
    }
    else if (topology.number() == 4) {

      nlohmann::json jm;

      const TopologyBase& t = problem.mesh_out.topology(0);

      jm["type"]     = mesh_topology.type_name();
      jm["dim"]      = t.points().dim();
      jm["number"]   = mesh_topology.number();
      jm["elements"] = t.data();
      jm["geometry"] = params["geometry"];
      jm["nb_ghost"] = t.points().nb_ghost();
      jm["vertices"] = t.points().data();


      std::ofstream file("mesh-adapt"+std::to_string(adapt_iter)+".avro");
      file << jm;
    }

    if (topology.number()==3 && has_interior_boundaries)
    {
      #if 0
      Boundary<Simplex> bnd2( mesh_topology );
      bnd2.extract(true);

      Topology<Simplex> bnd_topo( mesh_topology.points() , 2 );
      bnd2.retrieve( bnd_topo );

      library::Plottable<Simplex> plot(bnd_topo);
      gamma.withGeometry() = false;
      std::string boundary_file = directory + prefix + "_boundary_" + stringify(adapt_iter) + ".mesh";

      gamma.writeMesh( plot , boundary_file , false );
      printf("wrote boundary of mesh %s\n",boundary_file.c_str());

      gamma.writeMesh( problem.mesh_out , mesh_file , false );
      printf("wrote mesh %s\n",mesh_file.c_str());
      #endif
    }

    if (!has_interior_boundaries)
    {
      #if 0
      const std::string json_file = directory + "mesh_" + stringify(adapt_iter) + ".json";
      io::writeMesh<Simplex>( json_file , topology , &metric.field() );
      #endif
    }

  } // option to write the mesh

  // save the metric in case the caller wants information about it
  avro_assert_msg( field.nb() == topology.points().nb() , "|field| = %lu, nv = %lu", field.nb() , topology.points().nb() );
  fld.clear();
  for (index_t k=0;k<field.nb();k++)
    fld.push_back( field[k] );

  if (!output_redirect.empty())
  {
    fflush(stdout);
    dup2(redirected_fd,fileno(stdout));
    close(redirected_fd);
    clearerr(stdout);
    fsetpos(stdout,&pos);
  }

  avro_assert( problem.mesh_out.nb_topologies() == 1 );

  return result;
}

template class AdaptThread<Simplex>;
template int adapt<Simplex>( AdaptationProblem& );

} // avro
