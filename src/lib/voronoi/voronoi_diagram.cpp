//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "common/parallel_for.h"
#include "common/tools.h"

#include "geometry/entity.h"

#include "element/quadrature.h"

#include "mesh/decomposition.h"

#include "numerics/geometry.h"
#include "numerics/integration.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

#if defined(_OPENMP)
#include <omp.h>
#endif
#include <thread>

bool __check_capacity__ = true;

namespace avro
{

namespace delaunay
{

RVDFacets::RVDFacets( const Topology<Simplex>& topology ) :
  topology_(topology)
{
  create();
}

std::string
RVDFacets::generate( const std::vector<index_t>& f ) const
{
  std::string s;
  for (index_t j=0;j<f.size();j++)
  {
    s += stringify(f[j]);
    if (j<f.size()-1) s += "|";
  }
  return s;
}

int
RVDFacets::facet( const std::vector<index_t>& f ) const
{
  std::string s = generate(f);
  int id = 0;
  if (lookup(s,id)<=0)
    avro_assert_msg( lookup(s,id)>=0 , "facet does not exist: %s" , s.c_str() );
  return id;
}

int
RVDFacets::lookup( const std::string& s , int& id ) const
{
  std::map<std::string,int>::const_iterator it = store_.find(s);
  if (it==store_.end()) return -1;
  id = it->second;
  return 1;
}

void
RVDFacets::create()
{
  int id = 0;
  std::string s;

  // mesh facets
  coord_t nf = topology_.number()+1;
  std::vector<index_t> t,f;
  for (index_t k=0;k<topology_.nb();k++)
  {
    t = topology_.get(k);
    for (index_t j=0;j<nf;j++)
    {
      f = t;
      f.erase(f.begin()+j);
      std::sort(f.begin(),f.end());
      s = generate(f);
      if (lookup(s,id)<0)
      {
        id = -store_.size() -1;
        store_.insert( std::pair<std::string,int>(s,id) );
      }
    }
  }
}

void
RVDFacets::print() const
{
  std::map<std::string,int>::const_iterator it;
  printf("nb facets = %lu\n",store_.size());
  for (it=store_.begin();it!=store_.end();it++)
  {
    printf("facet[%d] = %s\n",it->second,it->first.c_str());
  }
}

} // delaunay

} // avro
