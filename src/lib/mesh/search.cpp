#include "common/kdtree.h"

#include "mesh/neighbours.h"
#include "mesh/search.h"
#include "mesh/topology.h"
#include "mesh/points.h"

#include "numerics/geometry.h"

#include <memory>

// tolerance used in the element search
// we may need to adjust this to deal with real geometries
#define TOLER 0.0

namespace luma
{

BoundarySearch::BoundarySearch( const Points& points ) :
  points_(points),
  searcher_(points_)
{}

index_t
BoundarySearch::nearest( real_t* x ) const
{
  std::vector<index_t> pts(1);
  searcher_.closestPoints(x,pts);
  return points_.index( pts[0] );
}

template<typename type>
ElementSearch<type>::ElementSearch( const Topology<type>& _topology ) :
      topology_(_topology), neighbours_(topology_.neighbours()),
      boundary_(topology_.points()),
      nb_brute_force_(0),
      nb_search_(0),
      time_(0),
      nb_steps_(0),
      brute_(true) // enable brute force search by default
{}

template<typename type>
int
ElementSearch<type>::find( real_t* x , const index_t guess )
{
  nb_search_++;

  // reset the tracker on which elements have been visited
  visited_.resize( topology_.nb() );
  std::fill( visited_.begin() , visited_.end() , false );

  // clear the path
  path_.clear();

  // initialize the depth and make sure we start with a non-ghost
  depth_ = 0;
  luma_assert( !topology_.ghost(guess) );

  clock_t t0 = clock();

  // step into the initial guess
  int result = step(x,guess);

  clock_t t1 = clock();

  time_ += (t1 -t0);

  return result;
}

template<typename type>
int
ElementSearch<type>::step( real_t* x , const index_t start )
{
  depth_++;
  nb_steps_++;

  luma_assert( !topology_.ghost(start) );

  // mark this element as visited
  visited_[start] = true;
  path_.push_back( start );

  // this is currently hacked for simplices
  index_t nf = neighbours_.nfacets(); // nfacets(start) in general meshes

  bool found = true;

  // loop through the facets
  for (index_t j=0;j<nf;j++)
  {
    // set the coordinate with this vertex
    std::vector<const real_t*> xj(nf);
    for (index_t i=0;i<nf;i++)
    {
      if (i==j) xj[i] = x;
      else xj[i] = topology_.points()[ topology_(start,i) ];
    }

    // check the barycentric volume
    real_t vol = numerics::simplex_volume( xj , topology_.points().dim() );

    if (vol<TOLER)
    {
      found = false;

      // get the j'th neighbour of this element
      int n = neighbours_(start,j);

      if (n<0)
      {
        // do not step into boundaries
        continue;
      }

      if (visited_[n])
      {
        // do not step into visited elements
        continue;
      }

      if (topology_.ghost(index_t(n)))
      {
        // do not step into ghosts
        visited_[n] = true;
        continue;
      }

      // step into the neighbour
      int k = step( x , index_t(n) );
      if (k>=0)
      {
        // we found the element, so return
        depth_--;
        return int(k);
      }
      else
      {
        // the element was not found in the recursive steps, so return the
        // same information back into the caller
        return -1;
      }
    } // vol < 0
    else
    {
      // volume is > 0 so continue checking the other facets
    }
  } // loop over facets

  // option to brute force the computation which searches the entire topology
  if (!found && brute_)
  {
    return brute(x);
  }

  // we were not able to find the point, it might be outside the topology
  // because we are searching for a projected geometry point
  if (!found) return -1;

  // we found the point in the starting element!
  return int(start);
}

template<typename type>
int
ElementSearch<type>::brute( real_t* x )
{
  nb_brute_force_++;

  const index_t nf = topology_.number()+1;
  std::vector<const real_t*> xk(nf,0);
  real_t vol;
  const coord_t dim = topology_.points().dim();

  for (index_t k=0;k<topology_.nb();k++)
  {
    if (topology_.ghost(k)) continue;

    nb_steps_++;

    bool found = true;
    for (index_t j=0;j<nf;j++)
    {
      for (index_t i=0;i<nf;i++)
      {
        if (i==j) xk[i] = x;
        else xk[i] = topology_.points()[ topology_(k,i) ];
      }
      // check the volume
      vol = numerics::simplex_volume(xk,dim);
      if (vol<TOLER)
      {
        found = false;
        break;
      }
    }
    // all the volumes were positive so the point was found
    if (found)
    {
      return int(k);
    }
  }

  // the point was not found, it is likely outside the topology
  // due to a geometry projection
  return -1;
}

template<typename type>
index_t
ElementSearch<type>::closest( real_t* x , std::vector<real_t>& alpha ) const
{
  // find the closest vertex on the boundary
  index_t q = boundary_.nearest( x );
  real_t dmin = numerics::distance2(topology_.points()[q] , x , topology_.points().dim() );

  #if 0
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    real_t d = numerics::distance2( topology_.points()[k] , x , topology_.points().dim() );
    if (d<dmin)
    {
      dmin = d;
      q    = k;
    }
  }
  dmin = numerics::distance2(topology_.points()[q] , x , topology_.points().dim() );
  #endif

  luma_assert( q >= topology_.points().nb_ghost() );

  // get the ball of the vertex
  std::vector<index_t> B;
  luma_assert( topology_.inverse().created() );
  topology_.inverse().ball( q , B );

  // look for the element with the closest barycenter
  index_t ielem = topology_.nb();

  std::vector<real_t> y( topology_.points().dim() , 0. );
  for (index_t k=0;k<B.size();k++)
  {
    // skip ghosts
    if (topology_.ghost(B[k])) continue;

    // determine the coordinates y closest to x
    real_t distance2 = topology_.master().closest( topology_.points() , topology_(B[k]) , topology_.nv(B[k]) , x , y );

    // compute the barycentric coordinates of y in this simplex
    std::vector<const real_t*> X(topology_.nv(B[k]));
    std::vector<real_t> alphak(topology_.nv(B[k]));
    topology_.get_elem( B[k] , X );
    numerics::barycentric_signed( y.data() , X , topology_.points().dim() , alphak );

    // make sure the point is indeed inside the simplex
    #if 0 // option to debug
    real_t alpha_tot = 0.0;
    for (index_t j=0;j<alphak.size();j++)
    {
      luma_assert_msg( alphak[j] >= -1e-12 , "alpha[%lu] = %g" , j,alphak[j] );
      alpha_tot += alphak[j];
    }
    luma_assert_msg( fabs(alpha_tot-1)<1e-12 , "sum(alpha)-1 = %g" , alpha_tot-1 );
    #endif

    if (distance2<dmin)
    {
      // save the closest element, distance and barycentric coordinates
      dmin  = distance2;
      ielem = B[k];
      alpha = alphak;
    }
  }

  // if no element minimizes the distance, then the original
  // vertex could be the minimizer, pick any element in the ball
  if (ielem==topology_.nb())
  {
    index_t k = 0;
    ielem = k;
    std::vector<const real_t*> X(topology_.nv(B[k]));
    topology_.get_elem( B[k] , X );
    numerics::barycentric_signed( topology_.points()[q] , X , topology_.points().dim() , alpha );
  }

  return ielem;
}

template<typename type>
void
ElementSearch<type>::print() const
{
  printf("nb searches: %lu\n",nb_search_);
  printf("total time:  %4.2g sec\n",timing());
  printf("\tnb brute force searches: %lu\n",nb_brute_force_);
  printf("\tnb neighbour stepping:   %lu\n",nb_search_-nb_brute_force_);
  printf("\ttotal steps:             %lu\n",nb_steps_);
}

PointSearch::PointSearch( Points& _points ) :
    points_(_points)
{
  cloud_ = std::make_shared<PointCloud>(points_);
  kdtree_ = initializeKdTree(*cloud_);
  kdtree_->build();
}

void
PointSearch::closestPoints( real_t* x , std::vector<index_t>& pts ) const
{
  index_t N = pts.size();
  std::vector<real_t> dist2(N);
  index_t nu;
  kdtree_->getNearestNeighbours(x,pts,dist2,nu);
}

void
PointSearch::insert( const real_t* x )
{
  luma_implement;
}

template class ElementSearch<Simplex>;


} // avro
