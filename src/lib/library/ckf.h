#ifndef LUNA_LIB_LIBRARY_CKF_H_
#define LUNA_LIB_LIBRARY_CKF_H_

#include "common/types.h"

#include "mesh/points.h"
#include "mesh/topology.h"

namespace luna
{

class CKF_Triangulation : public Topology<Simplex>
{
public:
  CKF_Triangulation( const std::vector<index_t>& dims );

  void generate();

private:
  void ndgrid( const std::vector<real_t>& dx , coord_t current );
  index_t find_vertex( const real_t* y ) const;

  Points points_;
  std::vector<real_t> p_;
  std::vector<index_t> u_;
  std::vector<index_t> dims_;

  Table<index_t> grid_;

};

} // luna

#endif