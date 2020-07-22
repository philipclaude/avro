//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "adaptation/cavity.h"

#include "common/tools.h"

#include "mesh/inverse.h"
#include "mesh/neighbours.h"
#include "mesh/points.h"
#include "mesh/topology.h"

#include <unistd.h>

namespace avro
{

template<typename type>
InverseTopology<type>::InverseTopology( const Topology<type>& topology ) :
  topology_(topology)
{}

template<typename type>
index_t
InverseTopology<type>::nb() const
{
  return topology_.points().nb();
}

template<typename type>
bool
InverseTopology<type>::created() const
{
  return elem_.size()==topology_.points().nb();
}

template<typename type>
void
InverseTopology<type>::build()
{
  elem_.resize( topology_.points().nb() , topology_.nb()+1 );

  for (index_t k=0;k<topology_.nb();k++)
  {
    for (index_t j=0;j<topology_.nv(k);j++)
    {
      index_t p = topology_(k,j);

      // skip ghosts
      if (topology_.ghost(k) && p>=topology_.points().nb_ghost()) continue;

      //elem_[p] = std::min(k,elem_[p]);
      elem_[p] = k;
    }
  }

  // make sure all points are set
  for (index_t k=0;k<nb();k++)
  {
    if (k<topology_.points().nb_ghost()) continue;
    if (elem_[k]>= topology_.nb()+1) topology_.points().print(true);
    avro_assert_msg(elem_[k] < topology_.nb()+1 , "inverse not found for vertex %lu" , k );
    avro_assert_msg( topology_.has( elem_[k] , k ) , "element %lu does not have vertex %lu!" , elem_[k] , k );
  }
}

template<typename type>
index_t
InverseTopology<type>::find( index_t p ) const
{
  for (index_t k=0;k<topology_.nb();k++)
  {
    // skip ghosts
    if (topology_.ghost(k) && p>=topology_.points().nb_ghost()) continue;
    for (index_t j=0;j<topology_.nv(k);j++)
    {
      if (topology_(k,j)==p)
        return k;
    }
  }
  avro_assert_not_reached;
  return topology_.nb();
}

template<typename type>
void
InverseTopology<type>::ball( index_t p , std::vector<index_t>& B ) const
{
  std::unordered_set<index_t> b;
  getball( p , elem_[p] , b );

  B.resize( b.size() );
  std::unordered_set<index_t>::const_iterator it;
  index_t count = 0;
  for (it=b.begin();it!=b.end();it++)
    B[count++] = *it;
}

template<typename type>
void
InverseTopology<type>::getball( index_t p , index_t k0 , std::unordered_set<index_t>& b ) const
{
  const index_t nf = topology_.neighbours().nfacets();
  std::unordered_set<index_t>::const_iterator it;
  int k1;

  b.insert( k0 );

  // loop through the neighbours of k, skipping neighbours already in b
  for (index_t j=0;j<nf;j++)
  {
    k1 = topology_.neighbours()(k0,j);
    if (k1<0) continue;

    // skip if this neighbours exists in b
    if (b.find(k1)!=b.end()) { /*printf("exists!\n");*/ continue; }

    // skip if this neighbour doesn't contain p
    if (!topology_.has(k1,p)) { /*printf("elem %d doesn't have %lu!\n",k1,p);*/ continue; }

    // step into k1
    getball( p , k1 , b );
  }
}

template<typename type>
void
InverseTopology<type>::shell( index_t p , index_t q , std::vector<index_t>& S ) const
{
  std::vector<index_t> b;
  ball( p , b );
  int start = -1;
  for (index_t j=0;j<b.size();j++)
  {
    if (topology_.has(b[j],p) && topology_.has(b[j],q))
    {
      start = b[j];
      break;
    }
  }
  if (start<0)
  {
    try
    {
      // this is kind of a disastrous hack around potential non-manifold meshes
      // when really bad partitions are produced in parallel computations
      // todo: check for non-manifold points as a preprocessing step!
      //topology_.all_with( {p,q} , S );
      S.clear();return;
      shell(q,p,S);
      return;
    }
    catch(...)
    {
      avro_implement;
      topology_.all_with( {p,q} , S );
      return;
    }
    std::vector<index_t> bq;
    ball(q,bq);
    std::sort( b.begin() , b.end() );
    std::sort( bq.begin() , bq.end() );
    print_inline( b  , "ball of vertex "+stringify(p)+": ");
    print_inline( bq , "ball of vertex "+stringify(q)+": ");
    printf("elem[%lu] = %lu, elem[%lu] = %lu\n",p,elem_[p],q,elem_[q]);
    printf("nb points = %lu\n",topology_.points().nb());
    topology_.points().print(p,true);
    topology_.points().print(q,true);

    // search the mesh for a containing element
    std::vector<index_t> Bp,Bq;
    for (index_t k=0;k<topology_.nb();k++)
    {
      if (topology_.has(k,p)) Bp.push_back(k);
      if (topology_.has(k,q)) Bq.push_back(k);
      if (topology_.has(k,p) && topology_.has(k,q))
      {
        //start = k;
        printf("element %lu has (%lu,%lu)!\n",k,p,q);
        topology_.neighbours().print(k);
        printf("must be a bug!\n");
      }
    }
    std::sort( Bp.begin() , Bp.end() );
    std::sort( Bq.begin() , Bq.end() );
    print_inline( Bp , "recomputed ball of vertex " +stringify(p) +":" );
    print_inline( Bq , "recomputed ball of vertex " +stringify(q) +":" );
    for (index_t j=0;j<Bp.size();j++)
      topology_.neighbours().print(Bp[j]);
    for (index_t j=0;j<Bq.size();j++)
      topology_.neighbours().print(Bq[j]);
    topology_.neighbours().print();
    fflush(stdout);
  }
  avro_assert_msg( start>=0 ,
      "could not find element touching edge (%lu,%lu)",p,q);

  std::unordered_set<index_t> s;
  getshell( p , q , index_t(start) , s );

  S.resize( s.size() );
  std::unordered_set<index_t>::const_iterator it;
  index_t count = 0;
  for (it=s.begin();it!=s.end();it++)
    S[count++] = *it;
}

template<typename type>
void
InverseTopology<type>::shell( index_t t0 , index_t t1 , index_t t2 , std::vector<index_t>& S ) const
{
  std::vector<index_t> b;
  ball( t0 , b );
  int start = -1;
  for (index_t j=0;j<b.size();j++)
  {
    if (topology_.has(b[j],t0) && topology_.has(b[j],t1) && topology_.has(b[j],t2))
    {
      start = b[j];
      break;
    }
  }
  avro_assert_msg( start>=0 ,
      "could not find element touching triangle (%lu,%lu,%lu)",t0,t1,t2);

  std::unordered_set<index_t> s;
  getshell( t0 , t1 , t2 , index_t(start) , s );

  S.resize( s.size() );
  std::unordered_set<index_t>::const_iterator it;
  index_t count = 0;
  for (it=s.begin();it!=s.end();it++)
    S[count++] = *it;
}

template<typename type>
void
InverseTopology<type>::getshell( index_t p , index_t q , index_t k0 , std::unordered_set<index_t>& s ) const
{
  const index_t nf = topology_.neighbours().nfacets();
  std::unordered_set<index_t>::const_iterator it;
  int k1;

  s.insert( k0 );

  // loop through the neighbours of k, skipping neighbours already in b
  for (index_t j=0;j<nf;j++)
  {
    k1 = topology_.neighbours()(k0,j);
    if (k1<0) continue;

    // skip if this neighbours exists in b
    if (s.find(k1)!=s.end()) continue;

    // skip if this neighbour doesn't contain p or q
    if (!topology_.has(k1,p)) continue;
    if (!topology_.has(k1,q)) continue;

    // step into k1
    getshell( p , q , k1 , s );
  }
}

template<typename type>
void
InverseTopology<type>::getshell( index_t t0 , index_t t1 , index_t t2 , index_t k0 , std::unordered_set<index_t>& s ) const
{
  const index_t nf = topology_.neighbours().nfacets();
  std::unordered_set<index_t>::const_iterator it;
  int k1;

  s.insert( k0 );

  // loop through the neighbours of k, skipping neighbours already in b
  for (index_t j=0;j<nf;j++)
  {
    k1 = topology_.neighbours()(k0,j);
    if (k1<0) continue;

    // skip if this neighbour exists in b
    if (s.find(k1)!=s.end()) continue;

    // skip if this neighbour doesn't contain p or q
    if (!topology_.has(k1,t0)) continue;
    if (!topology_.has(k1,t1)) continue;
    if (!topology_.has(k1,t2)) continue;

    // step into k1
    getshell( t0 , t1 , t2 , k1 , s );
  }
}


template<typename type>
void
InverseTopology<type>::decrement( index_t bar )
{
  //#pragma omp parallel for
  for (index_t k=0;k<nb();k++)
  {
    if (elem_[k]>bar)
    {
      //printf("decrementing element on vertex %lu\n",k);
      elem_[k]--;
    }
  }
}

#if 1
template<typename  type>
void
InverseTopology<type>::update( Cavity<type>& cavity , bool delay )
{

  // adjust the element indices
  if (cavity.nb_insert()>=cavity.nb_cavity())
  {
    // if there are more elements inserted than removed,
    // then the cell numbers of existing cells do not change
    // nothing to do!
  }
  else
  {
    // some cells were removed, so the elements here should be decremented
    // determine the cell numbers of the removed elements
    std::vector<index_t> C( cavity.nb_cavity()-cavity.nb_insert() );
    //printf("number of elements removed = %lu\n",C.size());
    for (index_t k=cavity.nb_insert();k<cavity.nb_cavity();k++)
      C[k-cavity.nb_insert()] = cavity.cavity(k);

    // anything here higher than the elements in C will need to be decremented
    std::sort( C.begin() , C.end() );
    std::reverse( C.begin() , C.end() );
    for (index_t k=0;k<C.size();k++)
    {
      //printf("decrementing anything higher than %lu\n",C[k]);
      decrement( C[k] );
    }

  }

  // first check if any points are deleted
/*  if (!delay)
  {
    avro_assert_not_reached;
    for (index_t k=0;k<cavity.nb_removedNodes();k++)
      this->remove( cavity.removedNode(k) );
  }*/

  // get all the points in the insertion
  index_t elem,m;
  std::vector<index_t> N;
  std::map<index_t,index_t> N2elem;
  for (index_t k=0;k<cavity.nb_insert();k++)
  {
    elem = cavity.inserted()[k]; // index of the inserted element
    //printf("elem = %lu\n",elem);
    for (index_t j=0;j<topology_.nv(elem);j++)
    {
      m = topology_(elem,j);
      //if (m<topology_.points().nb_ghost()) continue;
      N.push_back( m );

      //if (!topology_.ghost(elem) && N2elem.find(m)==N2elem.end())
      if (N2elem.find(m)==N2elem.end())
      {
        N2elem.insert( std::pair<index_t,index_t>(m,elem) );
      }

    }
  }
  uniquify(N);

  if (N.size()!=N2elem.size())
  {
    //cavity.print();
    //print_inline( N , "points to update: " );
    //std::map<index_t,index_t>::const_iterator it;
    //for (it=N2elem.begin();it!=N2elem.end();it++)
    //  printf("N2elem = (%lu,%lu)\n",it->first,it->second);
  }

  /*std::map<index_t,index_t>::const_iterator it;
  for (it=N2elem.begin();it!=N2elem.end();it++)
    printf("N2elem = (%lu,%lu)\n",it->first,it->second);*/

  avro_assert( N.size()==N2elem.size() );
  for (index_t k=0;k<N.size();k++)
    avro_assert_msg( N2elem.find(N[k])!=N2elem.end() ,
                    "vertex %lu should have a non-ghost element" , N[k] );

  // update existing points
  for (index_t k=0;k<N.size();k++)
  {
    elem_[N[k]] = N2elem[N[k]];
  }

}
#endif

template<typename type>
void
InverseTopology<type>::create( index_t nb_new )
{
  for (index_t k=0;k<nb_new;k++)
    elem_.push_back(0); // empty
}


template<typename type>
void
InverseTopology<type>::remove( index_t k )
{
  elem_.erase( elem_.begin() + k );
}

template<typename type>
void
InverseTopology<type>::copy( const InverseTopology<type>& inverse )
{
  // check the points have the same size
  avro_assert( nb()==inverse.nb() );
  elem_.resize( nb() );
  for (index_t k=0;k<inverse.nb();k++)
    elem_[k] = inverse.elem(k);
}

template<typename type>
bool
InverseTopology<type>::check() const
{
  if (elem_.size()!=topology_.points().nb())
  {
    printf("bad sizes: |elem| = %lu , |V| = %lu\n",
            elem_.size(),topology_.points().nb());
    return false;
  }

  index_t noffenders = 0;
  std::vector<index_t> B0,B1;
  for (index_t k=0;k<topology_.points().nb();k++)
  {
    B0.clear();
    B1.clear();
    ball( k , B0 );
    topology_.all_with( {k} , B1 );

    std::sort( B0.begin() , B0.end() );
    std::sort( B1.begin() , B1.end() );

    if (B0.size()!=B1.size())
    {
      printf("vertex %lu:\n",k);
      printf("\t");print_inline(B0,"B0: ");
      printf("\t");print_inline(B1,"B1: ");
      noffenders++;
      continue;
    }

    bool ok = true;
    for (index_t j=0;j<B0.size();j++)
    {
      if (B0[j]!=B1[j])
      {
        ok = false;
        continue;
      }
    }

    if (!ok)
    {
      printf("vertex %lu:\n",k);
      printf("\t");print_inline(B0,"B0: ");
      printf("\t");print_inline(B1,"B1: ");
      noffenders++;
      continue;
    }
  }
  if (noffenders>0) return false;
  return true;
}

template<typename type>
void
InverseTopology<type>::print( bool balls) const
{
  std::vector<index_t> B;
  for (index_t k=0;k<nb();k++)
  {
    printf("vertex[%lu] -> elem %lu",k,elem_[k]);
    if (balls)
    {
      ball( k , B );
      print_inline( B , ", ball = ");
    }
    else printf("\n");
  }
}

template class InverseTopology<Simplex>;
template class InverseTopology<Polytope>;

} // avro
