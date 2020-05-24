#include "avro.h"

#include "programs.h"

#include "adaptation/metric.h"
#include "adaptation/properties.h"

#include "library/factory.h"
#include "library/meshb.h"
#include "library/metric.h"

#include "mesh/mesh.h"

#include "numerics/predicates.h"

typedef avro::real_t REAL;
#include <tetgen1.5.0/predicates.h>

#include <triangle/predicates.h>

#include <cmath>

namespace avro
{

namespace programs
{

int
conformity(int nb_input, const char** inputs)
{
  typedef Simplex type;

  int nb_expected = 3;
  if (nb_input<nb_expected || nb_input==-1)
  {
    printf("\t\tconformity [mesh] [metric] [output] [optional]\n");
    printf("\t\toptions:\n");
    printf("\t\t--> if metric is Uniform, must provide hmin=[real]\n");
    return 0;
  }

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();

  // initialize the predicates
  exactinit(1,0,0,10,10,10);
  exactinit();
  GEO::PCK::initialize();

  // options
  bool found;
  const char **options = inputs +nb_expected;
  int  nb_options = nb_input -nb_expected;
  UNUSED(found);
  UNUSED(options);
  UNUSED(nb_options);

  int nb_expected_elem = -1;
  found = parse(lookfor(options,nb_options,"nb_expected"),nb_expected_elem);

  // get the original input mesh
  std::string meshname( inputs[0] );
  std::shared_ptr<TopologyBase> ptopology = nullptr;
  std::shared_ptr<Mesh> pmesh0 = library::get_mesh(meshname,ptopology);

  // define the input mesh for the adaptation
  std::shared_ptr<Mesh> pmesh = std::make_shared<Mesh>(pmesh0->points().dim(),pmesh0->number());
  pmesh0->points().copy( pmesh->points() );

  // get the topology and add it to the input mesh
  Topology<type>& topology = *static_cast<Topology<type>*>(ptopology.get());
  topology.orient();
  pmesh->add(ptopology);

  // get the input metric
  std::string metricname( inputs[1] );
  std::shared_ptr<MetricAttachment> pfld;

  std::vector<real_t> metric_params;
  if (metricname=="Uniform")
    metric_params.push_back(topology.number());
  if (nb_input>4)
  {
    real_t hmin;
    found = parse(lookfor(options,nb_options,"hmin"),hmin);
    if (found) metric_params.push_back(hmin);

    if (metricname=="Uniform") avro_assert_msg(found,"specify hmin for uniform metric");
  }
  if (metricname=="Uniform")
    avro_assert_msg( metric_params.size()==2 , "specify hmin=[real] for uniform metric" );

  bool analytic_metric;
  pfld = library::get_metric(metricname,pmesh->points(),analytic_metric,metric_params);

  // create the discrete metric
  MetricField<Simplex> metric(topology,*pfld.get());

  // get the metric conformity statistics
  std::string outfile(inputs[2]);
  Properties properties(topology,metric);
  properties.dump(outfile);


  real_t lunit,qunit;
  index_t nb_elem;
  properties.conformity(lunit,qunit,nb_elem);
  printf("lunit = %g %%, qunit = %g %%, nb_elem = %lu (expected %d)\n",lunit,qunit,nb_elem,nb_expected_elem);
  if (nb_expected_elem < 0) nb_expected_elem = nb_elem; // skip check below

  real_t elem_diff = real_t(nb_elem) - real_t(nb_expected_elem);
  if (elem_diff<0) elem_diff *= -1.0;

  if (lunit > 90 && elem_diff/nb_expected_elem < 0.05 )
    return 0;
  return 1;
}

} // programs

} // avro
