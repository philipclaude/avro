// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef XFIELD_AVRO_H
#define XFIELD_AVRO_H

#ifndef SANS_AVRO
#error "XField_avro.h should only be included if SANS_AVRO is defined"
#endif

#include <vector>
#include <memory> //shared_ptr

#include "MPI/communicator_fwd.h"

#include "Field/XField.h"

#include <avro.h>

#include <boost/type_traits/remove_const.hpp>

namespace SANS
{

template< class BaseType >
class XFieldBase_avro : public BaseType
{
public:
  // these constructors will initialize a new avro context
  XFieldBase_avro( mpi::communicator& comm , int num , int dim , int udim );
  XFieldBase_avro( const std::shared_ptr<mpi::communicator> comm , int num , int dim , int udim );
  XFieldBase_avro( mpi::communicator& comm , const avro::Context& ctx );
  XFieldBase_avro( const std::shared_ptr<mpi::communicator> comm , const avro::Context& ctx );

protected:
  // these constructors will initialize the context using an incoming one
  XFieldBase_avro( avro::Context& context );
  XFieldBase_avro( std::shared_ptr<mpi::communicator> comm, avro::Context& context );

public:

  virtual ~XFieldBase_avro() {}

  std::vector<int>& vertexOnGeometry() { return vertexOnGeometry_; }

  // functions to go back and forth between ego/integer labeling
  // not sure these are needed anymore with the revised avro context
  //ego object( const int k ) const;
  //int label( ego object ) const;

  avro::Context& context() { return context_; }
  const avro::Context& context() const { return context_; }

  const std::map<int,int>& getBndGroupMap() const { return bndGroupMap_; }

protected:

  // constructs the ego to boundary group map
  void mapBoundaryGroups();

  mutable avro::Context context_;
  std::vector<int> vertexOnGeometry_;
  avro::coord_t udim_;
  std::vector<avro::real_t> paramOnGeometry_;
  std::vector<int> geometry_;

  const avro::real_t* paramOnGeometry(avro::index_t k) const { return paramOnGeometry_.data() + udim_*k; }

  // maps ego's in geometry_ to XField boundary group
  std::map<int,int> bndGroupMap_;
};

class XFieldEmpty {};

template< class PhysDim, class TopoDim >
class XField_avro : public XFieldBase_avro< XField<PhysDim,TopoDim> >
{
public:
  typedef XFieldBase_avro< XField<PhysDim,TopoDim> > BaseType;

  XField_avro( const XField_avro& xfld, const XFieldEmpty);
  XField_avro( mpi::communicator& comm, const avro::Context& context );
  XField_avro( const XField<PhysDim,TopoDim>& xfld, const avro::Context& context );

  void import( const avro::Context& context );
  void convert( avro::Context& context ) const;

  // use inverse evaluation to attach the XField to the avro::Model
  void attachToGeometry();

  virtual ~XField_avro() {}

  void fill( const std::vector<double>& params = {0.025,0.001,15.} );

  virtual const std::type_info& derivedTypeID() const { return typeid(XField_avro); }

protected:
  void importVertexOnGeometry(const XField_Lagrange<PhysDim>& xfld , const avro::Context& context );

  using BaseType::context_;
  using BaseType::vertexOnGeometry_;
  using BaseType::geometry_;
  using BaseType::udim_;
  using BaseType::paramOnGeometry_;
  using BaseType::bndGroupMap_;
};

} // SANS

#endif // XFIELD_AVRO_H
