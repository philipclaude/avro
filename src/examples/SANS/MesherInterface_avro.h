// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef MESHERINTERFACE_AVRO_H_
#define MESHERINTERFACE_AVRO_H_

//Python must be included first
#include "tools/PyDict.h"
#include "tools/Parameter.h"

#include "Field/XField.h"

#include <memory>
#include <cmath>

namespace SANS
{

template<class PhysDim,class TopoDim> class XField_avro;

template <class PhysDim, class TopoDim, class Mesher>
class MesherInterface;

class avroMesher;

struct avroParams : noncopyable
{
  const ParameterString FilenameBase{"FilenameBase", "tmp/", "Default filepath prefix for generated files"};
  const ParameterBool Curved{"Curved", true, "Mesh for curved geometries."};
  const ParameterBool DisableCall{"DisableCall", false, "Disable avroMesher call for debugging?"};
  const ParameterString BoundarySubdirectory
  { "BoundarySubdirectory", ".",
    "Default path relative to FilenameBase where avro should dump boundary files (only used for tesseract meshes)"
  };
  const ParameterNumeric<Real> InsertionVolumeFactor
  { "InsertionVolumeFactor", std::sqrt(2.0), -1, NO_LIMIT,
    "if factor*number_expected_simplices < number_created_simplices, the insertion is rejected...set to -1 to ignore this check"
  };
  const ParameterBool WriteMesh{"WriteMesh",true,"option to ask avro to write the output mesh"};
  const ParameterBool WriteConformity{"WriteConformity",false,"option to ask avro to write the output mesh-metric conformity"};
  const ParameterBool HasInteriorBoundaries
  {"HasInteriorBoundaries",false,
   "whether avro should extract interior boundary groups (used for wakes)"
  };
  const ParameterBool LimitInsertionLength
  {"LimitInsertionLength",true,
   "Whether splits are limited by checking if short edges are created (HIGHLY RECOMMENDED)"
  };
  const ParameterBool SwapOut
  {"SwapOut",true,
   "Whether swaps are used to get out of topological configurations that cause insertions/collapses to get rejected"
  };
  const ParameterNumeric<Real> MinInsertionLengthTarget
  {"MinInsertionLengthTarget",std::sqrt(2.0),1.4,2.0,
   "edge split target length on the second stage of insertions"
  };
  const ParameterNumeric<Real> MaxInsertionLengthTarget
  {"MaxInsertionLengthTarget",2.0 ,1.4,3.0,
   "edge split target length on the first stage of insertions"
  };
  const ParameterBool UseSmoothing{"UseSmoothing",true,"whether vertex smoothing is used (HIGHLY RECOMMENDED)"};
  const ParameterBool FefloaStyle{"FefloaStyle",false,"whether to (try to) emulate fefloa's algorithm"};

  static void checkInputs(PyDict d);
  static avroParams params;
};

//Class for transferring meshes and metric requests between sola and avro

template <class PhysDim, class TopoDim>
class MesherInterface<PhysDim, TopoDim, avroMesher>
{
public:
  static const int D = PhysDim::D;
  typedef DLA::MatrixSymS<D,Real> MatrixSym;

  MesherInterface(int adapt_iter, const PyDict& paramsDict);

  std::shared_ptr<XField<PhysDim, TopoDim>>
  adapt( const Field_CG_Cell<PhysDim,TopoDim,MatrixSym>& metric_request ,
         const XField<PhysDim,TopoDim>& mesh );

  bool conforms() const { return conforms_; }

protected:
  int adapt_iter_;
  const PyDict paramsDict_;

  bool conforms_;

};

}

#endif /* MESHERINTERFACE_AVRO_H_ */
