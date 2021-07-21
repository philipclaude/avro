#ifndef avro_LIB_ADAPTATION_PARALLEL_H_
#define avro_LIB_ADAPTATION_PARALLEL_H_

#include "avro_types.h"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace avro
{

namespace library
{
class MetricField_Analytic;
}

template<typename type> class Topology;
template<typename T> class symd;

typedef symd<real_t> VertexMetric;
class AdaptationParameters;
template<typename type> class Topology_Partition;
class Mesh;
template<typename type> class PartitionField;

template<typename type>
class AdaptationManager
{
public:
  AdaptationManager( const Topology<type>& topology , const std::vector<VertexMetric>& metric , AdaptationParameters& params );

  void set_analytic( library::MetricField_Analytic* analytic );

  // extracts the topology into the working topology (partitions if necessary)
  void initialize( const Topology<type>& topology , const std::vector<VertexMetric>& metrics );

  // fixes the boundary of the partition that is not on the geometry
  void fix_boundary();

  // runs a sequence of [optional: balance]  - adapt - migrate
  void adapt();

  // synchronizes global vertex numbers (at the interfaces) across all processors
  void synchronize();

  // retrieves the adapted and load-balanced partition
  void retrieve( Topology<type>& topology );

  // retrieves the total number of points
  index_t get_total_nb_points() const;

  // exchanges the elements between the partitions
  void exchange( std::vector<index_t>& repartition );

  const Topology_Partition<type>& topology() const { return topology_; }
  Topology_Partition<type>& topology() { return topology_; }
  void reassign_metrics( const std::vector<VertexMetric>& metrics );

  real_t time_process() const { return time_process_; }
  real_t time_synchronize() const { return time_synchronize_; }
  real_t time_adapt() const { return time_adapt_; }
  real_t time_migrate() const { return time_migrate_; }
  real_t time_partition() const { return time_partition_; }
  real_t time_exchange() const { return time_exchange_; }

private:

  void migrate_balance( index_t nb_part );
  void migrate_interface();

  void send_metrics( index_t receiver , const std::vector<index_t>& global_indices , const std::vector<VertexMetric>& metrics , bool global_flag );
  void receive_metrics( index_t sender , bool overwrite=false );

  void append_partition( index_t p , Topology<type>& topology , Topology_Partition<type>& partition ,
                         std::vector<bool>& created , std::vector<real_t>& partition_index ) const;

  AdaptationParameters& params_;
  Topology_Partition<type> topology_;
  std::vector<VertexMetric> metrics_;

  index_t rank_;
  std::shared_ptr<PartitionField<type>> field_;

  std::vector<index_t> element_offset_;
  std::set<std::pair<index_t,index_t>> pairs_;

  std::vector<index_t> crust_;

  library::MetricField_Analytic* analytic_;

  std::vector<bool> active_;

  real_t time_process_;
  real_t time_synchronize_;
  real_t time_adapt_;
  real_t time_migrate_;
  real_t time_partition_;
  real_t time_exchange_;
};

} // avro

#endif
