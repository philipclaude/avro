#ifndef avro_LIB_ADAPTATION_PARALLEL_H_
#define avro_LIB_ADAPTATION_PARALLEL_H_

#include "common/types.h"

#include "mesh/points.h"

#include "numerics/matrix.h"

#include <map>
#include <vector>

namespace avro
{

template<typename type> class Topology;

typedef numerics::SymMatrixD<real_t> VertexMetric;
class AdaptationParameters;
template<typename type> class Topology_Partition;
class Mesh;
template<typename type> class PartitionField;

template<typename type>
class AdaptationManager
{
public:
  AdaptationManager( const Topology<type>& topology , const std::vector<VertexMetric>& metric , AdaptationParameters& params );

  // extracts the topology into the working topology (partitions if necessary)
  void initialize( const Topology<type>& topology , const std::vector<VertexMetric>& metrics );

  // fixes the boundary of the partition that is not on the geometry
  void fix_boundary();

  // runs a sequence of [optional: balance]  - adapt - migrate
  void adapt();

  // balances all meshes across processors
  void balance(real_t alpha=-1 , real_t beta=-1);

  // migrates the interface mesh
  void migrate();

  // synchronizes global vertex numbers (at the interfaces) across all processors
  void synchronize();

  // analyzes the metric conformity on this processor,
  // and communicates with all other processors to determine if we are done
  bool analyze();

  // retrieves the adapted and load-balanced partition
  void retrieve( Topology<type>& topology );

  // retrieves the total number of points
  index_t get_total_nb_points() const;

private:

  void migrate_parmetis();
  void migrate_native();

  void send_metrics( index_t receiver , const std::vector<index_t>& global_indices , const std::vector<VertexMetric>& metrics , bool global_flag );
  void receive_metrics( index_t sender , bool overwrite=false );

  void append_partition( index_t p , Topology<type>& topology , Topology_Partition<type>& partition ,
                         std::vector<bool>& created , std::vector<real_t>& partition_index ) const;

  AdaptationParameters& params_;
  Topology_Partition<type> topology_;
  std::vector<VertexMetric> metrics_;

  index_t rank_;
  std::shared_ptr<PartitionField<type>> field_;
};

} // avro

#endif
