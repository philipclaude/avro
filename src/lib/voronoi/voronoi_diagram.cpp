#include "common/parallel_for.h"
#include "common/tools.h"

#include "geometry/entity.h"

#include "mesh/decomposition.h"

#include "numerics/geometry.h"

#include "voronoi/delaunay.h"
#include "voronoi/voronoi.h"

#if defined(_OPENMP)
#include <omp.h>
#endif
#include <thread>

namespace avro
{

namespace delaunay
{

RVDFacets::RVDFacets( const Topology<Simplex>& topology ) :
  topology_(topology)
{
  create();
}

RVDFacets::RVDFacets( const Topology<Simplex>& topology , const std::vector<index_t>& S ) :
  topology_(topology)
{
  create(S);
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
RVDFacets::create( const std::vector<index_t>& S )
{
  int id = 0;
  std::string s;

  // mesh facets
  coord_t nf = topology_.number()+1;
  std::vector<index_t> t,f;
  for (index_t k=0;k<S.size();k++)
  {
    t = topology_.get(S[k]);
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

RestrictedVoronoiDiagram::RestrictedVoronoiDiagram( const Topology<Simplex>& _mesh , Delaunay& _delaunay ) :
  Topology<Polytope>( vertices_ , _mesh.number() ),
  vertices_( _delaunay.dim() ),
  mesh_(_mesh) , delaunay_(_delaunay),
  neighbours_( delaunay_ ),
  parallel_(false), gpu_(false),
  outdir_(".")
{
  neighbours_.compute();
}

void
RestrictedVoronoiDiagram::compute( const bool exact )
{
  RVDFacets facets( mesh_ );
  simplices_.clear();

  // create all the restricted simplices
  simplices_.resize( mesh_.nb() );
  for (index_t k=0;k<mesh_.nb();k++)
  {
    simplices_[k] = std::make_shared<RestrictedVoronoiSimplex>( k , mesh_ ,
                        facets  , delaunay_ , neighbours_ , exact );
  }

  // dispatch the computation of all the simplices
  if (!parallel_)
  {
    for (index_t k=0;k<nb_simplices();k++)
      clip(k);
  }
  else
  {
    if (!gpu_)
    {
      ProcessCPU::parallel_for (
        parallel_for_member_callback( this , &thisclass::clip ),
        0,nb_simplices()
      );
    }
    else
    {
      // need to define the kernel for this
      avro_implement;
    }
  }

  printf("accumulating\n");

  // accumulate the result
  accumulate();
}

void
RestrictedVoronoiDiagram::compute( const std::vector<index_t>& S , const bool exact )
{
  RVDFacets facets( mesh_ , S );
  simplices_.clear();

  // create all the restricted simplices
  simplices_.resize( S.size() );
  for (index_t k=0;k<S.size();k++)
  {
    simplices_[k] = std::make_shared<RestrictedVoronoiSimplex>( S[k] , mesh_ ,
                        facets  , delaunay_ , neighbours_ , exact );
  }

  // dispatch the computation of all the simplices
  if (!parallel_)
  {
    for (index_t k=0;k<nb_simplices();k++)
      clip(k);
  }
  else
  {
    if (!gpu_)
    {
      ProcessCPU::parallel_for (
        parallel_for_member_callback( this , &thisclass::clip ),
        0,nb_simplices()
      );
    }
    else
    {
      // need to define the kernel for this
      avro_implement;
    }
  }

  // accumulate the result
  accumulate();

}

typedef struct
{
  std::vector<int> indices;
} SymbolicVertex;

// needed to create a set/map of elements
bool
operator==( const SymbolicVertex& fx , const SymbolicVertex& fy )
{
  // assumes fx and fy have the same topological dimension
  // and that the indices are sorted
  avro_assert( fx.indices.size()==fy.indices.size() );
  for (index_t j=0;j<fx.indices.size();j++)
    if (fx.indices[j]!=fy.indices[j])
      return false;
  return true;
}

// needed to create a map of elements
bool
operator<( const SymbolicVertex& f , const SymbolicVertex& g )
{
  // lexicographically compare the indices
  return std::lexicographical_compare(f.indices.begin(), f.indices.end(),
                                      g.indices.begin(), g.indices.end());
}

void
RestrictedVoronoiDiagram::accumulate()
{
  // clear the site information
  fields_.remove("sites");

  // clear the points
  points_.clear();
  this->child_.clear();
  this->Topology<Polytope>::clear();

  std::vector<index_t> sites;
  for (index_t k=0;k<simplices_.size();k++)
  {
    // add the cells
    for (index_t j=0;j<simplices_[k]->topology().nb();j++)
    {
      //std::vector<index_t> cell = simplices_[k]->topology()[j];
      std::vector<index_t> cell = simplices_[k]->topology().get(j);
      for (index_t i=0;i<cell.size();i++)
        cell[i] += points_.nb();
      add( cell.data() , cell.size() );
      sites.push_back( simplices_[k]->seed(j) );
    }

    // add the points
    for (index_t j=0;j<simplices_[k]->points().nb();j++)
      points_.create( simplices_[k]->points()[j] );

    // add the vertex-to-facet information
    const Table<int>& vf = simplices_[k]->topology().master().incidence();
    for (index_t j=0;j<vf.nb();j++)
    {
      std::vector<int> f = vf.get(j);
      points_.incidence().add( f.data() , f.size() );
    }
  }

  sites_ = std::make_shared<VoronoiSites>(*this);
  sites_->build();

  for (index_t k=0;k<sites.size();k++)
    sites_->value(k) = sites[k];

  fields_.make("sites",sites_);

  // merge the vertices
  #if 0
  std::map<SymbolicVertex,index_t> symbolic;

  std::vector<index_t> merge( points_.nb() );
  for (index_t k=0;k<points_.nb();k++)
  {
    merge[k] = k;
    SymbolicVertex s;
    s.indices = points_.incidence().get(k);
    std::sort( s.indices.begin() , s.indices.end() );
    if (symbolic.find(s)==symbolic.end())
      symbolic.insert( {s,k} );
    else
      merge[k] = symbolic[s];
  }

  for (index_t k=0;k<this->nb();k++)
  {
    for (index_t j=0;j<this->nv(k);j++)
      (*this)(k,j) = merge[ (*this)(k,j) ];
  }
  #endif
}

void
RestrictedVoronoiDiagram::computeCentroids( Points& centroids )
{

  // retrieve the site information
  //Topology<Polytope>& rvd = topology(0);
  Topology<Polytope>& rvd = (*this);

  // copy the points so the master can triangulate
  Points points0( rvd.points().dim() );
  rvd.points().copy( points0 );

  // initialize size of centroids
  centroids.clear();
  std::vector<real_t> x0( delaunay_.dim() , 0. );
  for (index_t k=0;k<delaunay_.nb();k++)
    centroids.create( x0.data() );

  std::vector<real_t> V( delaunay_.nb() , 0. );

  //rvd.Table<index_t>::print();
  //rvd.points().print();

  coord_t number = rvd.number();
  SimplicialDecomposition<Polytope> decomposition(rvd);
  decomposition.extract();

  std::vector<index_t> simplices;
  std::vector<index_t> parents;
  decomposition.get_simplices( rvd.number() , simplices , parents );

  index_t nb_simplices = simplices.size()/(number+1);

  const Simplex& master = decomposition.master();
  for (index_t k=0;k<nb_simplices;k++)
  {
    real_t vk = master.volume( decomposition.points() , &simplices[k*(number+1)] , number+1 );

    index_t s = sites_->value( parents[k] );
    V[s] += vk;

    // compute the centroid of this piece
    std::vector<real_t> c( delaunay_.dim() );
    numerics::centroid( &simplices[k*(number+1)] , number+1 , decomposition.points() , c );

    // add the contribution to x*V and V
    for (coord_t d=0;d<delaunay_.dim();d++)
      centroids[s][d] += c[d]*vk;
  }

  /*

  for (index_t k=0;k<rvd.nb();k++)
  {
    // get the site this piece is associate with
    index_t s = sites_->value(k);//sites_[k] -1;

    // compute the volume of this piece
    real_t vk = rvd.master().volume( points0 , rvd(k) , rvd.nv(k) );

    // compute the centroid of this piece
    std::vector<real_t> c( delaunay_.dim() );
    numerics::centroid( rvd(k) , rvd.nv(k) , rvd.points() , c );

    // add the contribution to x*V and V
    for (coord_t d=0;d<delaunay_.dim();d++)
      centroids[s][d] += c[d]*vk;

    V[s] += vk;
  }
  */

  // go back and divide by the volume
  for (index_t k=0;k<centroids.nb();k++)
  for (coord_t d=0;d<centroids.dim();d++)
  {
    // set the centroid to the original seed if the cell volume is zero
    if (V[k]<1e-12)
      centroids[k][d] = delaunay_[k][d];
    else
      centroids[k][d] /= V[k];
  }
}

real_t
cvtEnergy( const std::vector<real_t>& x , void* data )
{
  real_t* z = (real_t*) data;

  real_t e = 0.;
  for (coord_t d=0;d<x.size();d++)
    e += ( x[d] -z[d] )*( x[d] -z[d] );
  return e;
}

real_t
volume_integrand( const std::vector<real_t>& x , void* data )
{
  return 1.;
}

#if 0
real_t
RestrictedVoronoiDiagram::energy()
{
  // retrieve the topology where the rvd is stored
  Topology<Polytope>& rvd = topology(0);

  // copy the points so the master can triangulate
  Points points0( rvd.points().dim() );
  rvd.points().copy( points0 );

  // compute the CVT energy
  real_t E = 0.;
  for (index_t k=0;k<rvd.nb();k++)
  {
    // get the site this piece is associated with
    index_t s = 0;avro_implement;//sites_[k] -1;
    real_t* z = delaunay_[s];

    // integrate the deviation between the site and the centroid
    real_t dE = rvd.master().integrate( points0 , rvd(k) , rvd.nv(k) ,  &cvtEnergy , (void*) z );

    E += dE;
  }

  return E;
}

real_t
RestrictedVoronoiDiagram::energy_nd()
{
  // retrieve the topology where the rvd is stored
  Topology<Polytope>& rvd = topology(0);

  numerics::SimplexQuadrature quadrature( rvd.number() , rvd.points().dim() , 2 );

  // copy the points so the master can triangulate
  Points points0( rvd.points().dim() );
  rvd.points().copy( points0 );

  // compute the CVT energy
  real_t E = 0.;
  real_t V = 0.;
  for (index_t k=0;k<rvd.nb();k++)
  {
    // get the site this piece is associate with
    index_t s = 0;avro_implement;//sites_[k] -1;
    real_t* z = delaunay_[s];

    // integrate the deviation between the site and the centroid
    real_t dE = rvd.master().integrate( points0 , rvd(k) , rvd.nv(k) ,  &cvtEnergy , (void*) z , &quadrature );
    real_t dV = rvd.master().integrate( points0 , rvd(k) , rvd.nv(k) ,  &volume_integrand , NULL , &quadrature );

    E += dE;
    V += dV;
  }

  return E;
}
#endif

void
RestrictedVoronoiDiagram::optimise( const index_t nb_iter , bool exact , FILE* fid )
{
  // perform some iterations of Lloyd relaxation
  for (index_t iter=0;iter<nb_iter;iter++)
  {
    // compute the rvd
    compute(exact);

    // compute the centroids
    Points centroids( delaunay_.dim() );
    computeCentroids( centroids );

    real_t omega = 0.2;

    // move the points to the centroids
    real_t dx = 0;
    for (index_t k=0;k<centroids.nb();k++)
    {
      if (delaunay_.entity(k)) continue;
      for (coord_t d=0;d<centroids.dim();d++)
      {
        dx += pow( delaunay_[k][d] -centroids[k][d] , 2. );

        delaunay_[k][d] = delaunay_[k][d] +omega*(centroids[k][d] -delaunay_[k][d]);
      }
    }
    dx = std::sqrt(dx/delaunay_.nb());


    // recompute the nearest neighbour information
    neighbours_.compute();

    energy_ = 0;//avro_implement;//energy_nd();
    printf("iter[%lu]: e = %.6e, dx = %.6e\n",iter,energy_,dx);
    if (fid!=NULL)
      fprintf(fid,"%.12e\n",energy_);
  }
}

void
RestrictedVoronoiDiagram::extract( Topology<Simplex>& triangulation ) const
{
  // first make a unique set of points based on the vertex facet matrix
  const Table<int>& vfm = child(0).master().incidence();

  // find unique entries of the vfm
  std::unordered_set<std::string> labels;
  std::unordered_set<std::string> simplices;
  for (index_t k=0;k<vfm.nb();k++)
  {
    std::vector<int> bisectors = vfm.get(k);
    avro_assert(bisectors.size()==this->number_);

    bool meshfacet = false;
    for (index_t j=0;j<bisectors.size();j++)
    {
      if (bisectors[j]<0) meshfacet = true;
    }
    if (meshfacet) continue;

    std::sort(bisectors.begin(),bisectors.end());

    // create a string representation of these bisectors
    std::string s;
    for (index_t j=0;j<bisectors.size();j++)
    {
      s += stringify(bisectors[j]);
      if (j<bisectors.size()-1)
        s += "|";
    }

    if (labels.find(s)==labels.end())
    {
      labels.insert(s);

      // get the seeds of each bisector
      std::vector<index_t> simplex;
      for (index_t j=0;j<bisectors.size();j++)
      {
        index_t zi,zj;
        delaunay_.seeds( bisectors[j] , zi , zj );
        simplex.push_back(zi);
        simplex.push_back(zj);
      }
      uniquify(simplex);
      avro_assert(simplex.size()==index_t(this->number_+1));

      // create a label for this simplex
      std::string S;
      std::sort( simplex.begin() , simplex.end() );
      for (index_t j=0;j<simplex.size();j++)
      {
        S += stringify(simplex[j]);
        if (j<simplex.size()-1)
          S += "|";
      }
      if (simplices.find(S)==simplices.end())
      {
        simplices.insert(S);
        triangulation.add( simplex.data() , simplex.size() );
      }
    }
  }
  printf("RDT contains %lu simplices\n",triangulation.nb());
}

} // delaunay

} // avro
