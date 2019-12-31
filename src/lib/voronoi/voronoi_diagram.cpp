// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include "common/parallel_for.h"
#include "common/stringify.h"

#include "geometry/entity.h"

#include "graph/neighbours.h"

#include "mesh/delaunay/delaunay.h"
#include "mesh/delaunay/voronoi.h"

#include "mesh/geometrics.h"

#include "numerics/integrand.h"
#include "numerics/quadrature.h"
#include "numerics/tools.h"

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

RestrictedVoronoiDiagram::RestrictedVoronoiDiagram( const Topology<Simplex>& _mesh ,
    Delaunay& _delaunay ) :
  Mesh<ConvexPolytope>( _delaunay.dim() , _mesh.number() ),
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
    simplices_[k] = smart_new(RestrictedVoronoiSimplex)( k , mesh_ ,
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

void
RestrictedVoronoiDiagram::compute( const std::vector<index_t>& S , const bool exact )
{
  RVDFacets facets( mesh_ , S );
  simplices_.clear();

  // create all the restricted simplices
  simplices_.resize( S.size() );
  for (index_t k=0;k<S.size();k++)
  {
    simplices_[k] = smart_new(RestrictedVoronoiSimplex)( S[k] , mesh_ ,
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

void
RestrictedVoronoiDiagram::accumulate()
{
  // create the topology where the rvd will be stored
  smart_ptr(Topology<ConvexPolytope>) t = smart_new(Topology<ConvexPolytope>)
                                                      (vertices_,number_);

  // clear the site information
  fields_.remove("sites");
  sites_.clear();

  // clear the vertices
  vertices_.clear();
  topology_.clear();

  for (index_t k=0;k<simplices_.size();k++)
  {

    // add the cells
    for (index_t j=0;j<simplices_[k]->topology().nb();j++)
    {
      std::vector<index_t> cell = simplices_[k]->topology()[j];
      for (index_t i=0;i<cell.size();i++)
        cell[i] += vertices_.nb();
      t->add( cell.data() , cell.size() );
      sites_.addCell( real(simplices_[k]->seed(j)) );
    }

    // add the vertices
    for (index_t j=0;j<simplices_[k]->vertices().nb();j++)
      vertices_.create( simplices_[k]->vertices()[j] );

    // add the vertex facet information
    const Data<int>& vf = simplices_[k]->topology().master().vertexFacetMatrix();
    for (index_t j=0;j<vf.nb();j++)
    {
      std::vector<int> f = vf.get(j);
      t->master().vertexFacetMatrix().add( f.data() , f.size() );
    }

  }
  addTopology(t);
  fields_.make("sites",sites_);
  setFields();
}

void
RestrictedVoronoiDiagram::computeCentroids( Vertices& centroids )
{

  // retrieve the site information
  Topology<ConvexPolytope>& rvd = topology(0);

  // copy the vertices so the master can triangulate
  Vertices vertices0( rvd.vertices().dim() );
  rvd.vertices().copy( vertices0 );

  std::vector<real> V( delaunay_.nb() , 0. );

  // initialize size of centroids
  centroids.clear();
  std::vector<real> x0( delaunay_.dim() , 0. );
  for (index_t k=0;k<delaunay_.nb();k++)
    centroids.create( x0.data() );

  for (index_t k=0;k<rvd.nb();k++)
  {
    // get the site this piece is associate with
    index_t s = sites_[k] -1;

    // compute the volume of this piece
    real vk = rvd.master().volume( vertices0 , rvd(k) , rvd.nv(k) );

    // compute the centroid of this piece
    std::vector<real> c( delaunay_.dim() );
    geometrics::centroid( rvd(k) , rvd.nv(k) , rvd.vertices() , c );

    // add the contribution to x*V and V
    for (coord_t d=0;d<delaunay_.dim();d++)
      centroids[s][d] += c[d]*vk;

    V[s] += vk;
  }

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

real
cvtEnergy( const std::vector<real>& x , void* data )
{
  real* z = (real*) data;

  real e = 0.;
  for (coord_t d=0;d<x.size();d++)
    e += ( x[d] -z[d] )*( x[d] -z[d] );
  return e;
}

real
volume_integrand( const std::vector<real>& x , void* data )
{
  return 1.;
}

real
RestrictedVoronoiDiagram::energy()
{
  // retrieve the topology where the rvd is stored
  Topology<ConvexPolytope>& rvd = topology(0);

  // copy the vertices so the master can triangulate
  Vertices vertices0( rvd.vertices().dim() );
  rvd.vertices().copy( vertices0 );

  // compute the CVT energy
  real E = 0.;
  for (index_t k=0;k<rvd.nb();k++)
  {
    // get the site this piece is associated with
    index_t s = sites_[k] -1;
    real* z = delaunay_[s];

    // integrate the deviation between the site and the centroid
    real dE = rvd.master().integrate( vertices0 , rvd(k) , rvd.nv(k) ,  &cvtEnergy , (void*) z );

    E += dE;
  }

  return E;
}

real
RestrictedVoronoiDiagram::energy_nd()
{
  // retrieve the topology where the rvd is stored
  Topology<ConvexPolytope>& rvd = topology(0);

  numerics::SimplexQuadrature quadrature( rvd.number() , rvd.vertices().dim() , 2 );

  // copy the vertices so the master can triangulate
  Vertices vertices0( rvd.vertices().dim() );
  rvd.vertices().copy( vertices0 );

  // compute the CVT energy
  real E = 0.;
  real V = 0.;
  for (index_t k=0;k<rvd.nb();k++)
  {
    // get the site this piece is associate with
    index_t s = sites_[k] -1;
    real* z = delaunay_[s];

    // integrate the deviation between the site and the centroid
    real dE = rvd.master().integrate( vertices0 , rvd(k) , rvd.nv(k) ,  &cvtEnergy , (void*) z , &quadrature );
    real dV = rvd.master().integrate( vertices0 , rvd(k) , rvd.nv(k) ,  &volume_integrand , NULL , &quadrature );

    E += dE;
    V += dV;
  }

  return E;
}

void
RestrictedVoronoiDiagram::optimise( const index_t nb_iter , bool exact , FILE* fid )
{
  // perform some iterations of Lloyd relaxation
  for (index_t iter=0;iter<nb_iter;iter++)
  {

    // compute the rvd
    compute(exact);

    // compute the centroids
    Vertices centroids( delaunay_.dim() );
    computeCentroids( centroids );

    real omega = 0.2;

    // move the vertices to the centroids
    real dx = 0;
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

    energy_ = energy_nd();
    printf("iter[%lu]: e = %.6e, dx = %.6e\n",iter,energy_,dx);
    if (fid!=NULL)
      fprintf(fid,"%.12e\n",energy_);
  }
}

void
RestrictedVoronoiDiagram::extract( Topology<Simplex>& triangulation ) const
{
  // first make a unique set of vertices based on the vertex facet matrix
  Data<int>& vfm = topology_[0]->master().vertexFacetMatrix();

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
