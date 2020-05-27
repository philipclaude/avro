#ifndef avro_LIB_LIBRARY_MESHB_H_
#define avro_LIB_LIBRARY_MESHB_H_

#include "mesh/field.h"
#include "mesh/mesh.h"

#include <string>

namespace avro
{

namespace EGADS
{
class Model;
}

namespace library
{

class meshb : public Mesh
{
public:
  meshb( const std::string& filename , const EGADS::Model* model=nullptr );
  meshb() :
    Mesh(0,0)
  {}

  void read();
  void write( Mesh& mesh , const std::string& filename , bool with_bnd );
  template<typename type> void write( const Topology<type>& topology , const std::vector<index_t>& refs );

  index_t nv( const int GmfType ) const;

  template<typename type> void read_elements( int GmfType );

private:
  std::string filename_;
  int64_t fid_;
  int version_;

  const EGADS::Model* model_;

  std::shared_ptr<TopologyBase> main_topology_;

  std::map<int,index_t> ref_index_; // map from reference index to topology index
};

} // library

} // avro

#endif
