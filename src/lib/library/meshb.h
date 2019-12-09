#ifndef LUNA_LIB_LIBRARY_MESHB_H_
#define LUNA_LIB_LIBRARY_MESHB_H_

#include "mesh/field.h"
#include "mesh/mesh.h"

#include <string>

namespace luna
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

  void read();

  index_t nv( const int GmfType ) const;

  template<typename type> void read_elements( int GmfType );

private:
  std::string filename_;
  int64_t fid_;
  int version_;

  const EGADS::Model* model_;

  std::map<int,index_t> ref_index_; // map from reference index to topology index
};

} // library

} // luna

#endif
