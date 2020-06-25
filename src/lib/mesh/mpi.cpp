#include "geometry/entity.h"

#include "mesh/mpi_tags.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#ifdef AVRO_MPI

namespace avro
{

#if 0
template<typename type>
void
Topology<type>::send( mpi::communicator& comm , index_t receiver ) const
{
  mpi::send( mpi::blocking{} , this->data_ ,  receiver , TAG_CELL_INDEX );
  mpi::send( mpi::blocking{} , this->first_ , receiver , TAG_CELL_FIRST );
  mpi::send( mpi::blocking{} , this->last_ ,  receiver , TAG_CELL_LAST );
}

template<typename type>
void
Topology<type>::receive( mpi::communicator& comm , index_t sender )
{
  printf("receiving topology from sender %lu\n",sender);
  this->data_  = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_INDEX);
  this->first_ = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_FIRST);
  this->last_  = mpi::receive<std::vector<index_t>>(sender,TAG_CELL_LAST);

  avro_assert_not_reached;
  Table<index_t>::print();
}

template<typename type>
void
Topology<type>::send_points( mpi::communicator& comm , index_t receiver ) const
{

}
#endif

template class Topology<Simplex>;
template class Topology<Polytope>;

} // avro

#endif // AVRO_MPI
