// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include <iostream>
#include <sstream>
#include <vector>

#include "tools/SANSnumerics.h"     // Real
#include "Topology/Dimension.h"
#include "Topology/ElementTopology.h"

#include "MesherInterface_avro.h"

#include "XField_avro.h"

#include "Meshing/libMeshb/WriteSolution_libMeshb.h"

#include "Field/XFieldLine.h"
#include "Field/XFieldArea.h"
#include "Field/XFieldVolume.h"
#include "Field/XFieldSpacetime.h"

#include "Field/FieldLine_CG_Cell.h"
#include "Field/FieldArea_CG_Cell.h"
#include "Field/FieldVolume_CG_Cell.h"
#include "Field/FieldSpacetime_CG_Cell.h"

#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"
#include "MPI/MPI_sleep.h"

#include <avro.h>

#ifdef SANS_MPI
#include "MPI/serialize_DenseLinAlg_MatrixS.h"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/serialization/map.hpp>
#endif

namespace SANS
{

// cppcheck-suppress passedByValue
void avroParams::checkInputs(PyDict d)
{
  std::vector<const ParameterBase*> allParams;
  allParams.push_back(d.checkInputs(params.FilenameBase));
  allParams.push_back(d.checkInputs(params.Curved));
  allParams.push_back(d.checkInputs(params.DisableCall));
  allParams.push_back(d.checkInputs(params.BoundarySubdirectory));
  allParams.push_back(d.checkInputs(params.InsertionVolumeFactor));
  allParams.push_back(d.checkInputs(params.WriteMesh));
  allParams.push_back(d.checkInputs(params.WriteConformity));
  allParams.push_back(d.checkInputs(params.HasInteriorBoundaries));
  allParams.push_back(d.checkInputs(params.LimitInsertionLength));
  allParams.push_back(d.checkInputs(params.SwapOut));
  allParams.push_back(d.checkInputs(params.MinInsertionLengthTarget));
  allParams.push_back(d.checkInputs(params.MaxInsertionLengthTarget));
  allParams.push_back(d.checkInputs(params.UseSmoothing));
  allParams.push_back(d.checkInputs(params.FefloaStyle));
  d.checkUnknownInputs(allParams);
}
avroParams avroParams::params;

template <class PhysDim, class TopoDim>
MesherInterface<PhysDim, TopoDim, avroMesher>::MesherInterface(
  int adapt_iter, const PyDict& paramsDict ) :
  adapt_iter_(adapt_iter),
  paramsDict_(paramsDict),
  conforms_(false)
{}

template <class PhysDim, class TopoDim>
std::shared_ptr<XField<PhysDim, TopoDim>>
MesherInterface<PhysDim, TopoDim, avroMesher>::
adapt( const Field_CG_Cell<PhysDim,TopoDim,MatrixSym>& metric_request ,
       const XField<PhysDim,TopoDim>& mesh_in )
{
  SANS_ASSERT(mesh_in.derivedTypeID() == typeid(XField_avro<PhysDim,TopoDim>));

  // use the avro types
  using avro::coord_t;
  using avro::index_t;
  using avro::real_t;

  index_t number = TopoDim::D;
  index_t nb_rank = number*(number +1)/2;

  // cast the input mesh and retrieve the context
  const XField_avro<PhysDim,TopoDim>& mesh = static_cast<const XField_avro<PhysDim,TopoDim>&>(mesh_in);

  avro::Context context( mesh.context() );

  // create the output mesh
  std::shared_ptr< XField_avro<PhysDim,TopoDim> > mesh_out;
  mesh_out = std::make_shared< XField_avro<PhysDim,TopoDim>>(mesh, XFieldEmpty());

  // create the list of nodal metrics
#ifdef SANS_MPI
  int nDOF = metric_request.nDOFpossessed();

  nDOF = boost::mpi::all_reduce(*mesh.comm(), nDOF, std::plus<int>());
#else
  int nDOF = metric_request.nDOF();
#endif

  int comm_rank = mesh.comm()->rank();

  // size the metric field for avro
  std::vector<real_t> metrics;
  if (comm_rank == 0)
    metrics.resize( nDOF * nb_rank , 0.0 );

#ifdef SANS_MPI

  std::map<int,MatrixSym> buffer;

  for (int i = 0; i < metric_request.nDOFpossessed(); i++)
  {
    int iDOF_native = metric_request.local2nativeDOFmap(i);
    buffer[iDOF_native] = metric_request.DOF(i);
  }

  if (comm_rank == 0)
  {
    std::vector<std::map<int,MatrixSym>> bufferOnRank;
    boost::mpi::gather(*metric_request.comm(), buffer, bufferOnRank, 0 );

    // collapse down the buffer from all ranks (this remove any possible duplicates)
    for (std::size_t i = 0; i < bufferOnRank.size(); i++)
      buffer.insert(bufferOnRank[i].begin(), bufferOnRank[i].end());

    // Write out the metric values
    for ( const auto& DOFpair : buffer)
    {
      const MatrixSym& m = DOFpair.second;
      index_t count = 0;
      for (int i = 0; i < TopoDim::D; i++)
        for (int j = i; j < TopoDim::D; j++)
          metrics[ DOFpair.first*nb_rank + count++ ] = m(i,j);
          //fld[ DOFpair.first ](i,j) = m(i,j);
    }

  }
  else // send the buffer to rank 0
    boost::mpi::gather(*mesh.comm(), buffer, 0 );

#else

  index_t count = 0;
  for (int k=0;k<metric_request.nDOF();k++)
  {
    const MatrixSym& m = metric_request.DOF(k);
    for (int i=0;i<TopoDim::D;i++)
      for (int j=i;j<TopoDim::D;j++)
        metrics[count++] = m(i,j);
  }
#endif

  // setup some parameters
  avro::AdaptationParameters& params = context.parameters();
  std::string filename_base = paramsDict_.get(avroParams::params.FilenameBase);
  params.directory() = filename_base;
  params.adapt_iter() = adapt_iter_;
  params.curved() = paramsDict_.get(avroParams::params.Curved);
  params.insertion_volume_factor() = paramsDict_.get(avroParams::params.InsertionVolumeFactor);
  params.write_mesh() = true;//paramsDict_.get(avroParams::params.WriteMesh);
  params.write_conformity() = paramsDict_.get(avroParams::params.WriteConformity);
  params.has_interior_boundaries() = paramsDict_.get(avroParams::params.HasInteriorBoundaries);
  params.limit_insertion_length() = paramsDict_.get(avroParams::params.LimitInsertionLength);
  params.swapout() = paramsDict_.get(avroParams::params.SwapOut);
  params.lt_min() = paramsDict_.get(avroParams::params.MinInsertionLengthTarget);
  params.lt_max() = paramsDict_.get(avroParams::params.MaxInsertionLengthTarget);
  params.use_smoothing() = paramsDict_.get(avroParams::params.UseSmoothing);
  params.fefloa() = paramsDict_.get(avroParams::params.FefloaStyle);
  params.limit_metric() = true;
  params.has_uv() = true;

  // convert the XField to a mesh
  mesh.convert( context );

  if (mesh.comm()->rank()==0)
  {
    // adapt the mesh
    int result;
    try
    {
      result = context.adapt( metrics );
    }
    catch (const std::exception& e)
    {
      SANS_DEVELOPER_EXCEPTION(e.what());
    }
    conforms_ = true;

  } // rank = 0

  // wait for rank 0 (1 second at a time) to finish executing avro
  // using barrier causes waiting processor to use 100% CPU
  MPI_sleep(*metric_request.getXField().comm(), 0, 1000);

  // store the mesh from the context into the output XField
  mesh_out->import( context );

  return mesh_out;
}

//Explicit instantiations
template class MesherInterface<PhysD2, TopoD2, avroMesher>;
template class MesherInterface<PhysD3, TopoD3, avroMesher>;
template class MesherInterface<PhysD4, TopoD4, avroMesher>;

}
