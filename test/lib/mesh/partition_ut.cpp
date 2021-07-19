#include "unit_tester.hpp"

#include "common/mpi.hpp"
#include "common/process.h"

#include "geometry/egads/context.h"

#include "graphics/application.h"

#include "library/ckf.h"
#include "library/egads.h"

#include "mesh/field.hpp"
#include "mesh/partition.h"

using namespace avro;

UT_TEST_SUITE( mesh_partition_test_suite )

#if AVRO_MPI

class PartitionNumber : public Field<Simplex,real_t>
{
public:
  PartitionNumber( const Topology<Simplex>& topology , const std::vector<index_t>& parts ) :
    Field<Simplex,real_t>(topology,1,DISCONTINUOUS)
  {
    this->element().set_basis( BasisFunctionCategory_Lagrange );
    this->build();
    for (index_t k=0;k<topology.nb();k++)
    for (index_t j=0;j<topology.nv(k);j++)
      (*this)(k,j) = real_t(parts[k]);
  }

  index_t nb_rank() const override { return 1; }

  std::string get_name( index_t j ) const override
  {
    return std::to_string(j);
  }
};

UT_TEST_CASE( test1 )
{
  coord_t number = 2;

  int rank = mpi::rank();
  if (rank!=0) return;

  EGADS::Context context;
  std::vector<real_t> lens(number,1.);
  EGADS::Cube geometry(&context,lens);

  printf("creating topology\n");
  std::vector<index_t> dims(number,5);

  Points points(number);
  Topology<Simplex> topology( points , number );

  index_t nb_piece = 1;
  for (index_t i=0;i<nb_piece;i++)
  {
    index_t nb_points = points.nb();

    CKF_Triangulation topology1( dims );
    topology1.points().attach(geometry);

    for (index_t k=0;k<topology1.points().nb();k++)
    {
      topology1.points()[k][1] += 2*i; // offset
      points.create( topology1.points()[k] );
      points.set_entity( points.nb()-1 , topology1.points().entity(k) );
    }

    for (index_t k=0;k<topology1.nb();k++)
    {
      std::vector<index_t> tk(topology1(k),topology1(k)+topology1.nv(k));
      for (coord_t j=0;j<tk.size();j++)
        tk[j] += nb_points;
      topology.add( tk.data() , tk.size() );
    }
  }

  printf("computing neighbours\n");
  topology.neighbours().fromscratch() = true;
  topology.neighbours().compute();
  topology.element().set_basis( BasisFunctionCategory_Lagrange );

  Partition<Simplex> partition(topology);

  printf("partitioning mesh..\n");
  index_t npart = 2;
  partition.compute(npart);

  std::vector<std::shared_ptr<Topology_Partition<Simplex>>> parts(npart);
  partition.get(parts);

  std::shared_ptr< PartitionNumber > fld;
  fld = std::make_shared<PartitionNumber>(topology,partition.partition());
  topology.fields().make( "partition" , fld );

  graphics::Viewer vis;
  vis.add(topology);

  vis.run(AVRO_FULL_UNIT_TEST);

}
UT_TEST_CASE_END( test1 )

#endif

UT_TEST_SUITE_END( mesh_partition_test_suite )
