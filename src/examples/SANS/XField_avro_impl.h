// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// No include blocker here. This file should only be included once in a cpp file for explicit instantiation

#if !defined(XFIELD_AVRO_INSTANTIATE) && !defined(SANS_HEADERCOMPILE_CHECK)
#error "This file should only be included in a cpp file for explicit instantiations"
#endif

#include <map>
#include <sstream>

#include "XField_avro.h"

#include "tools/linspace.h"
#include "tools/minmax.h"
#include "tools/safe_at.h"

#include "Field/Partition/XField_Lagrange.h"
#include "Field/tools/for_each_CellGroup.h"
#include "Field/tools/for_each_BoundaryTraceGroup.h"

// This should always be in a cpp file or impl file
#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"

#ifdef SANS_MPI
#include "MPI/continuousElementMap.h"
#include "MPI/serialize_DenseLinAlg_MatrixS.h"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/all_gather.hpp>
#include <boost/mpi/collectives/scatter.hpp>
#include <boost/serialization/map.hpp>
#endif

namespace SANS
{

//---------------------------------------------------------------------------//
template< class BaseType >
XFieldBase_avro<BaseType>::XFieldBase_avro( mpi::communicator& comm , int num , int dim , int udim ) :
  BaseType(comm),
  context_(num,dim,udim),
  udim_(udim)
{
  context_.get_geometry_ids(geometry_);
}

//---------------------------------------------------------------------------//
template< class BaseType >
XFieldBase_avro<BaseType>::XFieldBase_avro( const std::shared_ptr<mpi::communicator> comm , int num , int dim , int udim ) :
  BaseType(comm),
  context_(num,dim,udim),
  udim_(udim)
{
  context_.get_geometry_ids(geometry_);
}

//---------------------------------------------------------------------------//
template< class BaseType >
XFieldBase_avro<BaseType>::XFieldBase_avro( mpi::communicator& comm , const avro::Context& ctx ) :
  BaseType(comm),
  context_(ctx),
  udim_(ctx.udim())
{
  context_.get_geometry_ids(geometry_);
}

//---------------------------------------------------------------------------//
template< class BaseType >
XFieldBase_avro<BaseType>::XFieldBase_avro( const std::shared_ptr<mpi::communicator> comm , const avro::Context& ctx ) :
  BaseType(comm),
  context_(ctx),
  udim_(ctx.udim())
{
  context_.get_geometry_ids(geometry_);
}


//---------------------------------------------------------------------------//
template< class BaseType >
XFieldBase_avro<BaseType>::XFieldBase_avro( avro::Context& context ) :
   context_(context),
   udim_(context_.udim())
{
  context_.get_geometry_ids(geometry_);
}

//---------------------------------------------------------------------------//
template< class BaseType >
XFieldBase_avro<BaseType>::XFieldBase_avro( std::shared_ptr<mpi::communicator> comm , avro::Context& context ) :
  BaseType(comm),
  context_(context),
  udim_(context_.udim())
{
  context_.get_geometry_ids(geometry_);
}

//----------------------------------------------------------------------------//
//  Maps geometry objects to XField boundary groups
template<class PhysDim>
class MapBoundaryGroups :
    public GroupFunctorBoundaryTraceType< MapBoundaryGroups<PhysDim> >
{
public:

  // Save off the boundary trace integrand and the residual vectors
  MapBoundaryGroups( const avro::Context& context ,
                    const std::vector<int>& geometry ,
                    const std::vector<int>& traceGroups,
                    std::map<int,int>& bndGroupMap ) :
   context_(context),
   geometry_(geometry),
   traceGroups_(traceGroups),
   bndGroupMap_(bndGroupMap) {}

  std::size_t nBoundaryTraceGroups() const          { return traceGroups_.size(); }
  std::size_t boundaryTraceGroup(const int n) const { return traceGroups_[n];     }

//----------------------------------------------------------------------------//
  // Integration function that integrates each element in the trace group
  template <class TopologyTrace>
  void
  apply( const typename XField<PhysDim, typename TopologyTrace::CellTopoDim>::template FieldTraceGroupType<TopologyTrace>& xfldTrace,
         const int traceGroupGlobal )
  {
    int nodeDOFmap[TopologyTrace::NNode];
    avro::index_t simplex[TopologyTrace::NNode];

    // loop just one element at the most
    const int nelem = MIN(xfldTrace.nElem(),1);
    for (int elem = 0; elem < nelem; elem++)
    {
      // get the global node element indices
      xfldTrace.associativity(elem).getNodeGlobalMapping(nodeDOFmap, TopologyTrace::NNode );

      // copy to avro data structure
      for (int i = 0; i < TopologyTrace::NNode; i++)
        simplex[i] = nodeDOFmap[i];

      // there must be an entity on the boundary...
      int gid = context_.facet_geometry( simplex , TopologyTrace::NNode );
      SANS_ASSERT( gid >= 0 );

      // find the geometry id that represents the boundary with this element
      bool found = false;
      for (std::size_t i = 0; i < geometry_.size(); i++)
      {
        if (geometry_[i] == gid)
        {
          bndGroupMap_[i] = traceGroupGlobal;
          found = true;
        }
      }
      SANS_ASSERT_MSG(found, "did not find a boundary trace group?");
    }
  }

protected:
  const avro::Context& context_;
  const std::vector<int>& geometry_;
  const std::vector<int>& traceGroups_;
  std::map<int,int>& bndGroupMap_;
};

//---------------------------------------------------------------------------//
template< class BaseType >
void
XFieldBase_avro<BaseType>::mapBoundaryGroups()
{
  // geometry must be set before generating the map
  SANS_ASSERT(geometry_.size() > 0);

  std::vector<int> boundaryTraceGroups;
  if (this->nBoundaryTraceGroups() > 0)
    boundaryTraceGroups = linspace(0,this->nBoundaryTraceGroups()-1);

  // assign the value to the trace
  for_each_BoundaryTraceGroup<typename BaseType::TopoDim>::apply(
      MapBoundaryGroups<typename BaseType::PhysDim>(context_, geometry_, boundaryTraceGroups, bndGroupMap_),
      *this );

#ifdef SANS_MPI
  // synchronize the map across all processors
  std::vector<std::map<int,int>> bndGroupMaps(this->comm_->size());
  boost::mpi::all_gather(*this->comm_, bndGroupMap_, bndGroupMaps );

  for (const std::map<int,int>& bnd : bndGroupMaps)
    bndGroupMap_.insert(bnd.begin(), bnd.end());
#endif
}

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
XField_avro<PhysDim, TopoDim>::XField_avro( mpi::communicator& comm, const avro::Context& ctx ) :
  BaseType(comm,ctx)
{}

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
XField_avro<PhysDim, TopoDim>::XField_avro( const XField_avro& xfld, const XFieldEmpty ) :
  BaseType(xfld.comm(),xfld.context())
{
  // set the boundary group map
  bndGroupMap_ = xfld.bndGroupMap_;
}

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
XField_avro<PhysDim, TopoDim>::XField_avro( const XField<PhysDim,TopoDim>& xfld, const avro::Context& ctx ) :
  BaseType(xfld.comm(),ctx)
{
  //amodel_->template initGeometry<PhysDim>();
  this->cloneFrom(xfld);
  attachToGeometry();
}

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
void
XField_avro<PhysDim, TopoDim>::attachToGeometry()
{
  using avro::coord_t;
  using avro::index_t;
  using avro::real_t;

  const coord_t dim = BaseType::D;

  // load the coordinates into the avro context
  index_t i = 0;
  std::vector<real_t> coordinates( this->nDOF()*dim );
  for (int k = 0; k< this->nDOF(); k++)
  {
    for (coord_t d = 0;d < dim; d++)
      coordinates[i++] = this->DOF(k)[d];
  }
  context_.load_coordinates(coordinates);

  std::vector<real_t> parameters;
  context_.attach_geometry();
  context_.retrieve_geometry( vertexOnGeometry_ , parameters );
  SANS_ASSERT_MSG( vertexOnGeometry_.size() == this->nDOF() , "|geometry| = %lu, nDOF = %d" , vertexOnGeometry_.size() , this->nDOF() );
  SANS_ASSERT( parameters.size() == udim_ * this->nDOF() );

  // save the parameters
  paramOnGeometry_.assign( parameters.begin() , parameters.end() );

  // map the boundary groups
  this->mapBoundaryGroups();
}

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
void
XField_avro<PhysDim, TopoDim>::import( const avro::Context& context )
{
  using avro::index_t;
  using avro::coord_t;
  using avro::real_t;

  SANS_ASSERT( PhysDim::D == context.dim() );

  XField_Lagrange<PhysDim> xfld(*this->comm());
  typename XField_Lagrange<PhysDim>::VectorX X;

  std::vector<index_t> simplices;
  std::vector<real_t> coordinates;
  context.retrieve_mesh( coordinates , simplices );

  index_t nVertex = coordinates.size() / context.dim();
  xfld.sizeDOF( nVertex );

  if (this->comm_->rank() == 0)
  {
    index_t i = 0;
    for (index_t k = 0; k < nVertex; k++)
    {
      for (coord_t d = 0; d < context.dim(); d++)
        X[d] = coordinates[i++];
      xfld.addDOF( X );
    }
  }

  index_t nv = context_.number()+1;
  index_t nf = context_.number();
  index_t nb_simplices = simplices.size() / nv;
  xfld.sizeCells( nb_simplices );

  // read the simplices
  std::vector<int> simplex(nv);

  int order = 1;
  int cellgroup = 0;

  if (this->comm_->rank() == 0)
  {
    for (index_t k=0;k<nb_simplices;k++)
    {
      for (index_t j=0;j<nv;j++)
        simplex[j] = simplices[k*nv+j];

      if (TopoDim::D == 2)
        xfld.addCell(cellgroup, eTriangle, order, simplex );
      else if (TopoDim::D == 3)
        xfld.addCell(cellgroup, eTet , order , simplex );
      else if (TopoDim::D == 4)
        xfld.addCell( cellgroup , ePentatope , order , simplex );
      else
        SANS_DEVELOPER_EXCEPTION( "unsupported simplex dimension" );
    }
  }

  // get the boundary groups
  std::vector<std::vector<index_t>> boundary;
  std::vector<int> geometry;
  context.retrieve_boundary(boundary,geometry);
  SANS_ASSERT( boundary.size() == geometry.size() );

  // count how many d-1 boundaries there are
  int nbnd = 0;
  for (index_t k=0;k<boundary.size();k++)
    nbnd += boundary[k].size() / nf;
  xfld.sizeBoundaryTrace(nbnd);

  if (this->comm_->rank() == 0)
  {
    if (bndGroupMap_.size() == 0)
    {
      #if 0
      // construct a boundary map based on the stacking bodies and using the bodyIndex of each
      // body within the stack
      std::vector<index_t> bodyOffset(context.nb_bodies(),0);
      for (index_t k = 0; k < boundary.size(); k++)
      {
        // determine the entity this boundary is on
        index_t id = geometry[k];

        std::size_t ibody = 0;
        for (ibody = 0; ibody < bodyOffset.size(); ibody++)
        {
          if (entity->body()->object() == amodel_->body(ibody)->object())
            break;
        }
        if (ibody == bodyOffset.size()-1) continue;

        // find the maximum body index
        bodyOffset[ibody+1] = MAX(bodyOffset[ibody+1], entity->bodyIndex()-1);
      }

      for (index_t k = 0; k < bnd.nb_children(); k++)
      {
        if (bnd.child(k)->number()!=topology.number()-1) continue;

        // determine the entity this boundary is on
        avro::Entity* entity = bnd.entity(k);

        std::size_t ibody = 0;
        for (ibody = 0; ibody < bodyOffset.size(); ibody++)
        {
          if (entity->body()->object() == amodel_->body(ibody)->object())
            break;
        }

        // stack the groups consistent with the facet numbering of the body offset by each body
        bndGroupMap_[ this->label(entity) ] = bodyOffset[ibody] + entity->bodyIndex()-1;
      }
      #else
      SANS_DEVELOPER_EXCEPTION("implement");
      #endif
    }

    std::vector<int> facet(nf);
    int group = 0;
    for (index_t k=0;k<boundary.size();k++)
    {
      // retrieve the boundary group of the associated geometry id
      group = safe_at( bndGroupMap_ , geometry[k] );

      index_t nb_facets = boundary[k].size() / nf;
      for (index_t j=0;j<nb_facets;j++)
      {
        // retrieve the facet indices
        for (index_t i=0;i<nf;i++)
          facet[i] = boundary[k][j*nf+i];

        if (TopoDim::D == 4)
          xfld.addBoundaryTrace(group, eTet , facet);
        else if (TopoDim::D == 3)
          xfld.addBoundaryTrace(group, eTriangle, facet);
        else if (TopoDim::D == 2)
          xfld.addBoundaryTrace(group, eLine , facet );
        else
          SANS_DEVELOPER_EXCEPTION( "unsupported simplex dimension" );
      }

      printf("-> imported boundary group %d\n",group);
    }
  }

  // construct the grid which will generate the local2nativeDOFmap
  this->buildFrom(xfld);
  this->importVertexOnGeometry(xfld,context);
}

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
void
XField_avro<PhysDim, TopoDim>::importVertexOnGeometry(const XField_Lagrange<PhysDim>& xfld,const avro::Context& context)
{
  // get the paramter dimension
  udim_ = context_.udim();

#ifdef SANS_MPI
  using avro::coord_t;
  using avro::index_t;
  using avro::real_t;

  SANS_ASSERT(udim_ == PhysDim::D-1);
  typedef DLA::VectorS<PhysDim::D-1,Real> VectorP;

  std::map<int,int> vertexOnGeometry_local;
  std::map<int,VectorP> paramOnGeometry_local;

  if (this->comm_->rank() == 0)
  {
    std::vector<int> g;
    std::vector<real_t> u;
    context.retrieve_geometry(g,u);
    const coord_t udim = context.udim();

    // send global2local maps to rank 0
    std::vector<std::map<int,int>> global2local(this->comm_->size());
    boost::mpi::gather(*this->comm_, xfld.global2localDOFmap(), global2local, 0 );

    // fill vertexOnGeometry for each processor
    std::vector<std::map<int,int>>     vertexOnGeometry_buffer(this->comm_->size());
    std::vector<std::map<int,VectorP>> paramOnGeometry_buffer(this->comm_->size());
    for (int rank = 0; rank < this->comm_->size(); rank++)
    {
      for (const std::pair<int,int>& gl : global2local[rank])
      {
        vertexOnGeometry_buffer[rank][gl.second] = g[gl.first];

        VectorP param;
        for (int d = 0; d < VectorP::M; d++)
          param[d] = u[gl.first*udim+d];

        paramOnGeometry_buffer[rank][gl.second] = param;
      }
    }

    // scatter back vertexOnGeometry to each processor
    boost::mpi::scatter(*this->comm_, vertexOnGeometry_buffer, vertexOnGeometry_local, 0 );
    boost::mpi::scatter(*this->comm_, paramOnGeometry_buffer , paramOnGeometry_local , 0 );
  }
  else // send the to and receive from rank 0
  {
    boost::mpi::gather(*this->comm_, xfld.global2localDOFmap(), 0 );

    boost::mpi::scatter(*this->comm_, vertexOnGeometry_local, 0 );
    boost::mpi::scatter(*this->comm_, paramOnGeometry_local , 0 );
  }

  // populate the local array
  vertexOnGeometry_.resize( vertexOnGeometry_local.size() );
  for (const std::pair<int,int>& vOnG : vertexOnGeometry_local)
    vertexOnGeometry_[vOnG.first] = vOnG.second;

  paramOnGeometry_.resize( udim_*paramOnGeometry_local.size() );
  for (const std::pair<int,VectorP>& pOnG : paramOnGeometry_local)
    for (coord_t d = 0; d < udim_; d++)
      paramOnGeometry_[udim_*pOnG.first+d] = pOnG.second[d];

#else

  context.retrieve_geometry( vertexOnGeometry_ , paramOnGeometry_ );

#endif
}

template <class PhysDim, class TopoDim>
class addCellGroupTopology : public GroupFunctorCellType< addCellGroupTopology<PhysDim,TopoDim> >
{
public:
  addCellGroupTopology( const XField_avro<PhysDim,TopoDim>& xfld,
                        avro::Context& context,
                        const int nCellGroups )
    : xfld_(xfld), context_(context), nCellGroups_(nCellGroups) {}

  std::size_t nCellGroups() const          { return nCellGroups_; }
  std::size_t cellGroup(const int n) const { return n;            }

  //----------------------------------------------------------------------------//
  // Function that is applied to each cell group
  template< class Topology >
  void
  apply(const typename XField<PhysDim,TopoDim>::template FieldCellGroupType<Topology>& xfldCell,
        const int cellGroupGlobal)
  {
    SANS_ASSERT( context_.number() == TopoDim::D );
    using avro::index_t;

    int nElem = xfldCell.nElem();
    int nBasis = xfldCell.nBasis();
    SANS_ASSERT( nBasis == TopoDim::D+1 );

    std::vector<int> nodeDOFmap( nBasis );
    std::vector<index_t> simplex( nBasis );

#ifdef SANS_MPI

    int comm_rank = xfld_.comm()->rank();

    // construct a continuous element ID map for the elements in the current group
    const std::vector<int>& cellIDs = xfld_.cellIDs(cellGroupGlobal);
    std::map<int,int> globalElemMap;
    continuousElementMap( *xfld_.comm(), cellIDs, xfldCell, nElem, globalElemMap);

    int comm_size = xfld_.comm()->size();
    int nElemLocal = xfldCell.nElem();

    // maximum element chunk size that rank 0 will write at any given time
    int Elemchunk = nElem / comm_size + nElem % comm_size;

    // send one chunk of elements at a time to rank 0
    for (int rank = 0; rank < comm_size; rank++)
    {
      int elemLow  =  rank   *Elemchunk;
      int elemHigh = (rank+1)*Elemchunk;

      std::map<int,std::vector<int>> buffer;

      for (int elem = 0; elem < nElemLocal; elem++)
      {
        int elemID = globalElemMap[cellIDs[elem]];

        if (elemID >= elemLow && elemID < elemHigh &&
            xfldCell.associativity(elem).rank() == comm_rank)
        {
          xfldCell.associativity(elem).getGlobalMapping( nodeDOFmap.data() , nodeDOFmap.size() );

          // transform back to native DOF indexing
          for (std::size_t i = 0; i < nodeDOFmap.size(); i++)
            nodeDOFmap[i] = xfld_.local2nativeDOFmap(nodeDOFmap[i]);

          buffer[cellIDs[elem]] = nodeDOFmap;
        }
      }

      if (comm_rank == 0)
      {
        std::vector<std::map<int,std::vector<int>>> bufferOnRank;
        boost::mpi::gather(*xfld_.comm(), buffer, bufferOnRank, 0 );

        // collapse down the buffer from all other ranks
        for (std::size_t i = 1; i < bufferOnRank.size(); i++)
          buffer.insert(bufferOnRank[i].begin(), bufferOnRank[i].end());

        // create the simplex and add it to the topology
        std::vector<index_t> simplices( nElem*nBasis );

        index_t count = 0;
        for ( const auto& nodesPair : buffer)
        {
          for (int n = 0; n < nBasis; n++)
            simplices[count++] = nodesPair.second[n];
        }
        context_.load_simplices(simplices);
      }
      else // send the buffer to rank 0
        boost::mpi::gather(*xfld_.comm(), buffer, 0 );

    } // loop over ranks

#else // SANS_MPI

    std::vector<index_t> simplices( nElem*nBasis );

    index_t count = 0;
    for (int elem=0;elem<nElem;elem++)
    {
      // get the global element indices
      xfldCell.associativity(elem).getGlobalMapping(nodeDOFmap.data(), nodeDOFmap.size() );
      for (int i=0;i<nBasis;i++)
        simplices[count++] = nodeDOFmap[i];
    }

    context_.load_simplices(simplices);

#endif // SANS_MPI

  }

private:
  const XField_avro<PhysDim,TopoDim>& xfld_;
  avro::Context& context_;
  int nCellGroups_;
};

//---------------------------------------------------------------------------//
template< class PhysDim, class TopoDim >
void
XField_avro<PhysDim, TopoDim>::convert( avro::Context& context ) const
{
  using avro::index_t;
  using avro::coord_t;
  using avro::real_t;

  const coord_t dim = context.dim();
  const coord_t udim = context.udim();

  std::vector<Real> x( PhysDim::D );
  SANS_ASSERT( dim == PhysDim::D );
  SANS_ASSERT( int(vertexOnGeometry_.size()) == this->nDOF() );

  index_t nDOFnative = this->nDOFnative();

#ifdef SANS_MPI
  typedef DLA::VectorS<PhysDim::D,Real> VectorX;
  typedef DLA::VectorS<PhysDim::D-1,Real> VectorP;

  std::map<int,VectorX> buffer_coord;
  std::map<int,int>     buffer_geom;
  std::map<int,VectorP> buffer_param;

  int comm_rank = this->comm()->rank();

  for (int i = 0; i < this->nDOF(); i++)
  {
    int iDOF_native = this->local2nativeDOFmap(i);
    buffer_coord[iDOF_native] = this->DOF(i);
    buffer_geom[iDOF_native]  = vertexOnGeometry_[i];

    VectorP param;
    for (int d = 0; d < VectorP::M; d++)
      param[d] = this->paramOnGeometry(i)[d];
    buffer_param[iDOF_native] = param;
  }

  if (comm_rank == 0)
  {
    std::vector<int> geometry( nDOFnative , -1 );
    std::vector<real_t> coordinates( nDOFnative*dim , 1e20 );
    std::vector<real_t> parameters( nDOFnative*udim , 1e20 );

    std::vector<std::map<int,VectorX>> bufferCoordOnRank;
    std::vector<std::map<int,int>>     bufferGeomOnRank;
    std::vector<std::map<int,VectorP>> bufferParamOnRank;
    boost::mpi::gather(*this->comm(), buffer_coord, bufferCoordOnRank, 0 );
    boost::mpi::gather(*this->comm(), buffer_geom,  bufferGeomOnRank , 0 );
    boost::mpi::gather(*this->comm(), buffer_param, bufferParamOnRank, 0 );

    // collapse down the buffer from all ranks (this removes duplicates)
    int nbuffer = bufferCoordOnRank.size();
    SANS_ASSERT( int(bufferGeomOnRank.size())==nbuffer );
    for (std::size_t i = 0; i < std::size_t(nbuffer); i++)
    {
      buffer_coord.insert(bufferCoordOnRank[i].begin(), bufferCoordOnRank[i].end());
      buffer_geom.insert( bufferGeomOnRank[i].begin(),  bufferGeomOnRank[i].end() );
      buffer_param.insert(bufferParamOnRank[i].begin(), bufferParamOnRank[i].end());
    }

    // create the vertices
    for ( const auto& DOFpair : buffer_coord)
    {
      // set the coordinates
      for (coord_t d=0;d<dim;d++)
        coordinates[ DOFpair.first * dim + d] = DOFpair.second[d];
    }

    // attach the geometry
    for ( const auto& geom : buffer_geom )
    {
      if (geom.second>=0)
        geometry[ geom.first ] = vertexOnGeometry_[geom.first];
    }

    // set the paramter values
    for ( const auto& param : buffer_param )
    for ( coord_t d=0;d<udim;d++)
      parameters[ param.first * udim + d ] = param.second[d];

    // copy the vertex information into the avro context
    context.load_coordinates(coordinates);
    context.load_geometry( geometry , parameters );

  }
  else // send the buffer to rank 0
  {
    boost::mpi::gather(*this->comm(), buffer_coord, 0 );
    boost::mpi::gather(*this->comm(), buffer_geom,  0 );
    boost::mpi::gather(*this->comm(), buffer_param, 0 );
  }

#else // SANS_MPI

  std::vector<int> geometry( nDOFnative , -1 );
  std::vector<real_t> coordinates( nDOFnative*dim , 1e20 );
  std::vector<real_t> parameters( nDOFnative*udim , 1e20 );

  // create the vertex information
  for (index_t k = 0; k < index_t(this->nDOFnative()); k++)
  {
    // set the coordinates
    for (coord_t d = 0; d < dim; d++)
      coordinates[k*dim+d] = this->DOF(k)[d];

    // set the geometry info
    for (coord_t d = 0; d < udim; d++)
      parameters[k*udim+d] = this->paramOnGeometry(k)[d];
    geometry[k] = vertexOnGeometry_[k];
  }

  // copy the vertex information into the avro context
  context.load_coordinates(coordinates);
  context.load_geometry( geometry , parameters );

#endif // SANS_MPI

  // loop over cell groups and write out the connectivity
  int nCellGroup = this->nCellGroups();
  for_each_CellGroup<TopoDim>::apply( addCellGroupTopology<PhysDim, TopoDim>(*this,context, nCellGroup), *this );

  // put a barrier for processor communication
  this->comm()->barrier();
}

} // SANS
