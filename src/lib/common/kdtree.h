#ifndef avro_LIB_COMMON_KDTREE_H_
#define avro_LIB_COMMON_KDTREE_H_

#include "common/error.h"
#include "common/types.h"

#include <memory>
#include <vector>
#include <math.h>

#include <nanoflann/nanoflann.hpp>

namespace avro
{

class Points;

class PointCloud
{
public:
  PointCloud( Points& _vertices );

  // these functions must be defined for nanoflann
	size_t kdtree_get_point_count() const;
	real_t kdtree_get_pt(const size_t idx, int d) const;
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
  coord_t dim() const;

private:
  Points& points_;

};

class KdTreeNd
{
public:
  KdTreeNd( PointCloud& _cloud ) :
    cloud_(_cloud)
  {}

  // it's ok to have "build" as a virtual function because it gets called once
  // but to have" getNearestNeighbours" as a virtual function is
  // probably a bad idea for performance, but honestly i'm tired and need
  // something that works
  // update: timing results show that it doesn't really matter
  // n.b. that good compilers these days will probably optimize the vtable
  // lookup of the virtual method, meaning it won't be too expensive
  virtual void build() = 0;
  virtual void getNearestNeighbours( const real_t* q, std::vector<index_t>& idx ,
                                     std::vector<real_t>& dist2 ,
                                     index_t& nu ) = 0;

  virtual ~KdTreeNd() {}

protected:
  PointCloud& cloud_;
};

template<coord_t dim>
class KdTree : public KdTreeNd
{
public:
  KdTree( PointCloud& _cloud );
  void build();
  void getNearestNeighbours( const real_t* q, std::vector<index_t>& idx ,
                             std::vector<real_t>& dist2 , index_t& nu );

private:
  std::shared_ptr<nanoflann::KDTreeSingleIndexAdaptor
      <nanoflann::L2_Simple_Adaptor<real_t,PointCloud>,PointCloud,dim>> tree_;
};

inline std::shared_ptr<KdTreeNd>
initializeKdTree( PointCloud& cloud )
{
  // instantiate this for many dims
  // remember the embedded RVD calculation might need high dimensional kdtree
  const coord_t dim = cloud.dim();
  switch (dim)
  {
    case 0:
      avro_assert_not_reached;
    case 1:
      return std::make_shared<KdTree<1>>(cloud);
    case 2:
      return std::make_shared<KdTree<2>>(cloud);
    case 3:
      return std::make_shared<KdTree<3>>(cloud);
    case 4:
      return std::make_shared<KdTree<4>>(cloud);
    case 5:
      return std::make_shared<KdTree<5>>(cloud);
    case 6:
      return std::make_shared<KdTree<6>>(cloud);
    case 7:
      return std::make_shared<KdTree<7>>(cloud);
    case 8:
      return std::make_shared<KdTree<8>>(cloud);
    case 9:
      return std::make_shared<KdTree<9>>(cloud);
    case 10:
      return std::make_shared<KdTree<10>>(cloud);
    case 11:
      return std::make_shared<KdTree<11>>(cloud);
    case 12:
      return std::make_shared<KdTree<12>>(cloud);
    default:
      avro_implement;
  }
  avro_assert_not_reached;
  return std::make_shared<KdTree<3>>(cloud);
}

} // avro

#endif
