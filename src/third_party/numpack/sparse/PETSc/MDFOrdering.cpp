// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "MDFOrdering.h"

#include <algorithm>

namespace numpack 
{
namespace SLA
{

struct MDFOrderingNode
{
  MDFOrderingNode() : weight(0), position(-1) {}
  MDFOrderingNode(const PetscScalar& w, const int& pos) : weight(w), position(pos) {}

  PetscScalar weight;
  int position;
};

bool operator >(const MDFOrderingNode& x, const MDFOrderingNode& y)
{
  if (x.weight == y.weight)
    return (x.position > y.position); //tie-breaker
  else
    return (x.weight > y.weight);
}

std::vector<PetscInt>
computeOrdering_MDF(const SparseMatrix_CRS<PetscScalar>& C)
{
  //-----------------------------------------------
  //Compute fill weights

  const int nrow = C.m();
  SANS_ASSERT(C.m() == C.n());

  //min-heap data structure
  std::vector<MDFOrderingNode> weights(nrow);
//  std::priority_queue<MDFOrderingNode, std::vector<MDFOrderingNode>, std::greater<MDFOrderingNode> > weights;

  for (int k = 0; k < nrow; k++)
  {
    PetscScalar weight = 0.0;

    const int row_k_nz = C.rowNonZero(k);
    for (int ki = 0; ki < row_k_nz; ki++)
    {
      int i = C.columnIndex(k,ki);
      if (k == i) continue;  //looking for neighbors of k, so skip diagonal
      //std::cout << "ki: " << ki << ", i: " << i << std::endl;

      if (!C.isNonZero(i,k)) continue; //skip if i,k is zero
      PetscScalar Cik = C.slowAccess(i,k);

      for (int kj = 0; kj < row_k_nz; kj++)
      {
        int j = C.columnIndex(k,kj);
        if (k == j) continue;  //looking for neighbors of k, so skip diagonal
        if (i == j) continue;  //need two unique neighbors of k, so skip if i = j
        //std::cout << "  kj: " << kj << ", j: " << j << std::endl;

        if (C.isNonZero(j,i))
          continue; //skip if i and j are connected

        //std::cout << "k: " << k << ", i: " << i << ", j: " << j << std::endl;
        PetscScalar dC = Cik * C.sparseRow(k,kj); //dC = C(i,k)*C(k,j)
        weight += dC*dC; //discarded fill weight
      }
    }
    weight = sqrt(weight);

    //add current position and weight
    weights[k].position = k;
    weights[k].weight = weight;
    //std::cout << "weight[" << k << "] = " << weight << std::endl;
  }

  //Make a min-heap of the weights list
  std::greater<MDFOrderingNode> comparator;
  std::make_heap(weights.begin(), weights.end(), comparator);

  //for (int m = 0; m < nrow; m++)
  //  std::cout << "weights[" << m << "]: " << weights[m].position << ", " << weights[m].weight << std::endl;

  //-----------------------------------------------
  //Compute MDF ordering

  std::vector<PetscInt> ordering(nrow, -1);
  std::vector<bool> isOrdered(nrow, false);

  for (int m = 0; m < nrow; m++)
  {
    const MDFOrderingNode min_node = weights.front(); //get the node with the minimum weight
    const int p = min_node.position;                  //Pivot with the smallest fill-in
    ordering[m] = p;                                  //save off ordering
    isOrdered[p] = true;                              //keep track of which original indices have already been ordered

    //std::cout << "min-node[" << m << "]: " << min_node.position << ", " << min_node.weight << std::endl;

    //remove min node, this is equivalent to setting weights[p] = infinity
    std::pop_heap(weights.begin(), weights.end(), comparator); //remove min node from heap and move to end of vector
    weights.pop_back(); //remove last element

    //for (std::size_t i = 0; i < weights.size(); i++)
    //  std::cout << "weight[" << i << "] = " << weights[i].position << ", " << weights[i].weight << std::endl;

    const int row_p_nz = C.rowNonZero(p);
    for (int pk = 0; pk < row_p_nz; pk++)
    {
      int k = C.columnIndex(p,pk); //k is a neighbor of p
      if (isOrdered[k]) continue;  //only looking for neighbors of p that are not already ordered

      //Recompute weights for neighbors of p that are not already ordered
      PetscScalar weight = 0.0;

//      std::cout << "k = " << k << std::endl;

      const int row_k_nz = C.rowNonZero(k);
      for (int ki = 0; ki < row_k_nz; ki++)
      {
        int i = C.columnIndex(k,ki);
        if (k == i) continue;       //looking for neighbors of k, so skip diagonal
        if (isOrdered[i]) continue; //skip if neighbor i is already ordered
        SANS_ASSERT(i < nrow);

        if (!C.isNonZero(i,k)) continue; //skip if i,k is zero
        PetscScalar Cik = C.slowAccess(i,k);
        for (int kj = 0; kj < row_k_nz; kj++)
        {
          int j = C.columnIndex(k,kj);
          if (k == j) continue;  //looking for neighbors of k, so skip diagonal
          if (i == j) continue;  //need two unique neighbors of k, so skip if i = j
          if (isOrdered[j]) continue; //skip if neighbor j is already ordered

          if (C.isNonZero(j,i))
            continue; //skip if i and j are connected

          PetscScalar dC = Cik * C.sparseRow(k,kj); //dC = C(i,k)*C(k,j)
          weight += dC*dC; //discarded fill weight
        }
      }
      weight = sqrt(weight);

      bool updated = false;
      for (MDFOrderingNode& node : weights)
      {
        if (node.position == k)
        {
          node.weight = weight;
          updated = true;
          break;
        }
      }
      SANS_ASSERT(updated);
      std::make_heap(weights.begin(), weights.end(), comparator); //update heap

//      std::cout << "min-node2[" << m << "]: " << weights.front().position << ", " << weights.front().weight << std::endl;
//      for (std::size_t i = 0; i < weights.size(); i++)
//        std::cout << "weight2[" << i << "] = " << weights[i].position << ", " << weights[i].weight << std::endl;

    } //k-neighbor loop
  }

  //for (std::size_t i = 0; i < ordering.size(); i++)
  //  std::cout << "ordering[" << i << "] = " << ordering[i] << std::endl;

  return ordering;
}

} //namespace SLA
} //namespace numpack 
