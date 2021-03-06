//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "programs.h"

#include "adaptation/adapt.h"
#include "adaptation/metric.h"

#include "common/directory.h"
#include "common/tools.h"

#include "geometry/model.h"

#include "graphics/application.h"

#include "library/factory.h"
#include "library/library.h"
#include "library/tesseract.h"

#include "mesh/mesh.h"
#include "avro.h"

#include <stdio.h>
#include <time.h>

namespace avro
{

namespace programs
{

int
adapt( int nb_input , const char** inputs ) {

  // so far only simplex adaptation is supported
  typedef Simplex type;

  if (nb_input < 4 || nb_input == -1) {
    printf("\t\tadapt [input mesh] [geometry] [metric] [output prefix] [optional]\n");
    printf("\t\t--> optional can be:\n");
    printf("\t\t\tcurved=(bool, whether the geometry is curved)\n");
    printf("\t\t\tlimit=(bool, whether target metric is limited from implied mesh metric)\n");
    printf("\t\t\tnb_iter=(int, # iterations for analytic metric only)\n");
    printf("\t\t\thmin=(real, minimum size used by Tesseract Linear case only)\n");
    printf("\t\t\tscale=(real, scale on the metric field entries)\n");
    return 1;
  }

  // options
  bool found;
  const char **options = inputs +4;
  int  nb_options = nb_input -4;

  bool analytic_metric = true;
  bool curved = true;
  coord_t number = 0;

  // option to compute the implied metric
  bool compute_implied = false;
  if (analytic_metric) compute_implied = true;
  if (nb_input>4)
    found = parse<bool>(lookfor(options,nb_options,"limit"),compute_implied);

  real_t href = 2.0;
  if (nb_input > 4)
    found = parse<real_t>(lookfor(options,nb_options,"href"),href);

  // if the metric is analytic, iterate...
  index_t nb_iter = 10;
  if (nb_input > 4)
    found = parse(lookfor(options,nb_options,"nb_iter"),nb_iter);

  // if the metric is analytic, iterate...
  real_t scale = 1.0;
  if (nb_input > 4)
    found = parse(lookfor(options,nb_options,"scale"),scale);

  // get the original input mesh
  std::string meshname( inputs[0] );
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh0 = library::get_mesh(meshname,ptopology);

  // define the input mesh for the adaptation
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(pmesh0->points().dim(),pmesh0->number());
  pmesh0->points().copy( pmesh->points() );
  number = pmesh->number();

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();
  pmesh->add(ptopology);

  // get the input geometry
  std::string geometryname( inputs[1] );
  std::shared_ptr<Model> pmodel;
  pmodel = library::get_geometry( geometryname , curved );
  Model& model = *pmodel;

  // option to overwrite the curved geometry parameter
  found = parse(lookfor(options,nb_options,"curved"),curved);

  // check the points are on the geometry...
  // option to project them
  pmesh->points().attach(model);

  // get the input metric
  std::string metricname( inputs[2] );
  std::shared_ptr<MetricAttachment> pfld;

  // check if we should map the points
  if (number == 4) {

    library::Tesseract& body = static_cast<library::Tesseract&>(pmodel->body(0));
    std::string map_name = lookfor(options,nb_options,"map");
    if (!map_name.empty()) {

      if (map_name == "expansion") {
        body.map_to( &library::Tesseract::expansion , &pmesh->points() );
      }
      else if (map_name == "closingwall") {
        body.map_to( &library::Tesseract::closingwall , &pmesh->points() );
      }
      else {
        printf("unknown map\n");
        avro_assert_not_reached;
      }
      geometryname += "-" + map_name;
    }
  }

  std::vector<real_t> metric_params;
  if (metricname=="Uniform")
    metric_params.push_back(number);
  if (nb_input > 4) {
    real_t hmin;
    found = parse(lookfor(options,nb_options,"hmin"),hmin);
    if (found) metric_params.push_back(hmin);

    if (metricname=="Uniform") avro_assert_msg(found,"specify hmin for uniform metric");
  }
  if (metricname=="Uniform")
    avro_assert_msg( metric_params.size()==2 , "specify hmin=[real] for uniform metric" );

  // define the problem and adapt
  AdaptationParameters params;
  params.set("curved",curved);
  params.set("directory", "./" );
  params.set("write conformity", false);
  params.set("swapout" , false);
  params.set("use smoothing" , true);
  params.set("geometry" , geometryname );

  std::string outputfile(inputs[3]);
  std::vector<std::string> s = split(outputfile,".");
  params.set("prefix", s[0]);

  if (number<4)
    params.set("insertion volume factor", -1.0 );

  clock_t TIME0 = clock();
  for (index_t iter = 0; iter < nb_iter; iter++) {

    // get the mesh and the metric field
    Mesh& mesh = *pmesh;
    pfld = library::get_metric(metricname,mesh.points(),analytic_metric,metric_params);

    // option to do the adaptation from the implied metric
    if (compute_implied)
      pfld->limit( mesh.retrieve<type>(0) , href );

    // adjust the adaptation parameters
    params.set("adapt iter" , iter);

    // create the mesh we will write to
    std::shared_ptr<Mesh> pmesh_out = std::make_shared<Mesh>(mesh.points().dim(),mesh.number());

    // define the problem and adapt
    std::vector<symd<real_t>> metrics( pfld->nb() );
    for (index_t k=0;k<metrics.size();k++)
      metrics[k] = (*pfld)[k]*scale;
    AdaptationProblem problem = {*pmesh,metrics,params,*pmesh_out};
    ::avro::adapt<type>( problem );

    // only one iteration for discrete metrics
    if (!analytic_metric) break;

    // create the mesh for the next adaptation
    pmesh = std::make_shared<Mesh>(number,number);
    pmesh_out->points().copy( pmesh->points() );
    avro_assert( pmesh->points().nb_ghost() == 0 );

    std::shared_ptr<Topology<Simplex>> ptopology = std::make_shared<Topology<Simplex>>(pmesh->points(),number);
    Topology<Simplex>& ptopology_out = pmesh_out->retrieve<type>(0);
    for (index_t k=0;k<ptopology_out.nb();k++)
      ptopology->add(ptopology_out(k),pmesh_out->topology(0).nv(k));
    pmesh->add(ptopology);

    if (params["curved"]) params.set("has uv", true);
  }
  clock_t TIME1 = clock();
  printf("--> time = %g seconds\n",real_t(TIME1-TIME0)/real_t(CLOCKS_PER_SEC));

  Library* lib = Library::get();
  lib->add_mesh_ptr(pmesh);
  char mesh_label[128];
  sprintf(mesh_label,"%p",(void*)pmesh.get());
  lib->add_mesh(mesh_label);

  return 0;
}

} // programs

} // avro
