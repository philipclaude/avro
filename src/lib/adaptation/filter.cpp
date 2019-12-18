#include "adaptation/filter.h"

#include "common/kdtree.h"

#include <algorithm>

namespace luna
{

Filter::Filter( const coord_t dim ) :
  Points(dim),
  lmin_(0.),
  lmax_(0.)
{
}

void
Filter::createPermanent( const real_t* x )
{
  node0_.push_back(-1);
  node1_.push_back(-1);
  s_.push_back( -1. );
  idx_.push_back( nb() );
  Points::create(x);
}

void
Filter::createCandidate( const index_t n0 , const index_t n1 ,
                         const real_t s , const real_t* x,
                         Entity * e, const real_t* params )
{
  luna_assert( n0 < n1 );
  node0_.push_back(n0);
  node1_.push_back(n1);
  s_.push_back(s);
  candidates_.push_back(nb());
  idx_.push_back(-1);
  Points::create(x);
  Points::set_entity(nb()-1,e);
  Points::set_param(nb()-1,params);
}

void
Filter::buildTree()
{
  // setup the kdtree
  cloud_ = std::make_shared<PointCloud>(*this);
  tree_ = initializeKdTree(*cloud_);
  tree_->build();
}

void
Filter::accept( const index_t k , const index_t idx )
{
  luna_assert_msg( node0_[k]>=0 && node1_[k]>=0 ,
      "node %lu is already active",k);

  // save the original values
  int n0 = node0_[k];
  int n1 = node1_[k];

  // look for any other candidate which touches this edge
  for (index_t i=0;i<nb();i++)
  {
    // permanent points do not get modified
    if (permanent(i)) continue;

    if (node0_[i]!=n0) continue;
    if (node1_[i]!=n1) continue;

    // check the s value for which node gets modified
    if (s_[i]<s_[k])
    {
      // node1 gets modified such that the edges is now [node0,idx]
      node1_[i] = idx;
    }
    else
    {
      // node1 gets modified such that the edge is now [idx,node1]
      node0_[i] = idx;
    }
  }

  // turn off the edges, which informs us that the node is active
  node0_[k] = -1;
  node1_[k] = -1;
  s_[k] = -1;
  idx_[k] = idx;

}

void
Filter::clearCandidates()
{
  for (int k=nb()-1;k>=0;k--)
  {
    if (permanent(k)) continue;

    node0_.erase( node0_.begin() +k );
    node1_.erase( node1_.begin() +k );
    s_.erase( s_.begin() +k );
    Points::remove( k );
  }
  candidates_.clear();
}

void
Filter::print() const
{
  for (index_t k=0;k<nb_candidates();k++)
  {
    index_t idx = candidate(k);
    printf("candidate[%lu] = %lu, node0 = %d, node1 = %d, s = %g\n",
                              k,idx,node0_[idx],node1_[idx],s_[idx]);
  }
}

} // luna
