
#include "element/simplex.h"
#include "mesh/topology.h"

#include "voronoi/new/facets.h"

namespace avro
{

namespace voronoi
{

bool
operator==( const Bisector& bx , const Bisector& by ) {
  if (bx.p0 != by.p0) return false;
  if (bx.p1 != by.p1) return false;
  return true;
}

// needed to create a map of elements
bool
operator<( const Bisector& f , const Bisector& g ) {
  if (f.p0 < g.p0) return true;
  if (f.p0 == g.p0) {
    if (f.p1 < g.p1) return true;
  }
  return false;
}


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

} // voronoi

} // avro
