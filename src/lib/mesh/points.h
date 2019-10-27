#ifndef URSA_MESH_points_H_
#define URSA_MESH_points_H_

#include "common/error.h"
#include "common/json.h"
#include "common/types.h"

#include <memory>
#include <vector>

namespace ursa
{

class Body;
class Model;

namespace geometrics
{
  class Primitive;
}

template<typename type> class Data;
template<typename type> class Topology;

class Points
{
public:
  Points() : dim_(0), udim_(0), ghost_(0) {}
  Points( const coord_t dim , const coord_t udim );
  Points( const coord_t dim );
  ~Points();

  void copy( Points& v , const bool erase=true , const bool ghosts=true ) const;

  coord_t dim() const { return dim_; }
  coord_t udim() const { return udim_; }
  index_t nb() const { return x_.size() == 0 ? 0 : x_.size()/dim_; }
  void setDimension( const coord_t _dim ) { dim_ = _dim; }
  void setParameterDimension( const coord_t _udim ) { udim_ = _udim; }

  // vertex creation
  void create( const std::vector<real_t>& x );
  void create( const real_t* x );

  // coordinate access
  const real_t* operator[] (const index_t k) const { return x_.data()+k*dim_; }
  real_t* operator[] (const index_t k) { return x_.data()+k*dim_; }
  //const real_t* get( const index_t k ) const { return operator[](k); } // TODO is this used anymore?
  const real_t& operator() (const index_t k, const index_t i ) const
    { return x_[k*dim_+i]; }

  // parameter coordinate access
  const real_t* u( const index_t k ) const { return u_.data()+k*udim_; }
  real_t* u( const index_t k ) { return u_.data()+k*udim_; }
  real_t& u( const index_t k , const index_t i )
    { return u_[k*udim_+i]; }

  void print( std::string pre="\0" , bool info=false ) const;
  void print( const index_t k , bool info=false ) const;

	void clear();

  // read/write body function
  int& body( const index_t k );
  void setPrimitive( const index_t k , geometrics::Primitive* e );
  void setParam( const index_t k , const std::vector<real_t>& u );
  void setParam( const index_t k , const real_t* u );
  void setFixed( const index_t k , bool f ) { fixed_[k] = f; }

  geometrics::Primitive* primitive( const index_t k ) const { return primitive_[k]; }
  bool fixed( const index_t k ) const { return fixed_[k]; }

  // boundary vertex query function
  bool boundary( const index_t k ) const;

  void createGhost();
  bool containsGhost() const { return ghost_>0; }
  bool ghost( const index_t k ) const { return k<ghost_; }
  index_t nb_ghost() const { return ghost_; }

  // manipulation functions
  void remove( const index_t k );

  // duplicates detection (useful when the mesh may have come from an OBJ file)
  // or given a vertex-facet incidence matrix F so that merging is done symbolically
  void duplicates( std::vector<index_t>& idx ,real_t tol=1e-16 ) const;
  void duplicates( std::vector<index_t>& idx , const Data<int>& F ) const;

  void dump( const std::string& filename ) const;

  void findGeometry( const Body& body , index_t ibody=1 ,real_t tol=1e-12 );
  void findGeometry( const Model& model ,real_t tol=1e-12 );
  void projectToGeometry( Body& body );

  // compute the partitioning
  void computePartition( Data<index_t>& data ,
                      std::vector<index_t>& cell_partition , index_t nparts );

  // retrieve the partition of where the vertex lives
  index_t partition(  const index_t k ) const
  {
    ursa_assert( k < partition_.size() );
    return partition_[k];
  }

  const std::vector<real_t>& coordinates() const { return x_; }

  void toJSON( json& J ) const;
  void fromJSON( const json& J , const Model* model=NULL );

 real_t INFTY = 1.0e20;

protected:

  coord_t dim_;  // dimension of ambient space the vertices live in
  coord_t udim_; // dimension of the parameter space of the geometry

  std::vector<real_t> x_; // vertex coordinates
  std::vector<real_t> u_; // geometry parameter coordinates
  std::vector<geometrics::Primitive*> primitive_; // which geometric primitive this vertex is on
  std::vector<int> body_; // which body this vertex lies on (<0 means it is a partition boundary, 0 = volume point, >0 means it is on a geometry body)
  std::vector<bool> fixed_; // flag if the vertex is fixed

  std::vector<index_t> periodic_; // which vertex coordinate this one is
                                  // periodic "with"
                                  // 0 = not periodic,
                                  // otherwise periodic_[k]-1 = vertex index

  index_t ghost_; // how many ghost vertices in the set

  std::vector<index_t> partition_; // which partition this vertex is on
};

} // ursa

#endif
