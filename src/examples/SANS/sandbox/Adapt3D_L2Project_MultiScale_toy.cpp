// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

// Adapt3D_AD_TripleBoundaryLayer_btest
// Testing of the MOESS framework on the advection-diffusion pde

//#define BOUNDARYOUTPUT

#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp> // to automagically make the directories

#include "SANS_btest.h"

#include "tools/SANSnumerics.h"     // Real

#include <iostream>

#include "pyrite_fstream.h"

#include "LinearAlgebra/DLA/StaticSize/MatrixS_Det.h"

#include "pde/AnalyticFunction/ScalarFunction3D.h"
#include "pde/AdvectionDiffusion/AdvectionDiffusion_Traits.h"

#include "pde/NDConvert/OutputNDConvertSpace3D.h"
#include "pde/NDConvert/SolnNDConvertSpace3D.h"
#include "pde/NDConvert/FunctionNDConvertSpace3D.h"

#include "pde/OutputCell_SolutionSquared.h"
#include "pde/OutputCell_SolutionErrorSquared.h"

#include "Discretization/Galerkin/IntegrandCell_Galerkin_Output.h"
#include "Discretization/Galerkin/AlgebraicEquationSet_Project.h"

#include "Discretization/Galerkin/FunctionalCell_Galerkin.h"
#include "Discretization/Galerkin/IntegrandCell_ProjectFunction.h"
#include "Discretization/Galerkin/IntegrandCell_Project.h"

#include "Adaptation/callMesher.h"
#include "Adaptation/MeshAdapter.h"

#include "Field/FieldVolume_CG_Cell.h"

#include "LinearAlgebra/SLA/SparseLinAlg_LinearSolver.h"
#include "NonLinearSolver/NewtonSolver.h"

#include "Meshing/refine/MultiScale_metric.h"
#include "Meshing/libMeshb/WriteMesh_libMeshb.h"
#include "Meshing/libMeshb/WriteSolution_libMeshb.h"

#include "unit/UnitGrids/XField3D_Box_Tet_X1.h"

#define MPI_COMMUNICATOR_IN_CPP
#include "MPI/communicator.h"

#ifdef SANS_AVRO
#include "Meshing/avro/XField_avro.h"
#include "Meshing/libMeshb/WriteMesh_libMeshb.h"
#endif

#include "Python/PyDict.h"

using namespace std;

//Explicitly instantiate the classes so that coverage information is correct
namespace SANS
{

}

using namespace SANS;

//############################################################################//
BOOST_AUTO_TEST_SUITE( Adapt3D_CG_AD_L2_MultiScale_test_suite )

//----------------------------------------------------------------------------//
BOOST_AUTO_TEST_CASE( Adapt3D_CG_AD_L2_MultiScale_test )
{
#define USE_TRIPLEBL 1
#if USE_TRIPLEBL
  typedef ScalarFunction3D_TripleBL SolutionExact;
  std::string func = "tripleBL";
#else
  typedef ScalarFunction3D_Tanh3 SolutionExact;
  std::string func = "tanh3";
  //typedef ScalarFunction3D_SinATan3 SolutionExact;
  //std::string func = "sinatan3";
  //typedef ScalarFunction3D_SinFun3 SolutionExact;
  //std::string func = "sinfun3";
#endif

  typedef IntegrandCell_ProjectFunction<FunctionNDConvertSpace<PhysD3,Real>,IntegrandCell_ProjFcn_detail::FcnX> IntegrandCellFunctionClass;

  typedef AlgebraicEquationSet_Project< XField<PhysD3, TopoD3>,
      IntegrandCellFunctionClass, TopoD3, AlgEqSetTraits_Sparse > ProjectionEquationSetClass;

  typedef ProjectionEquationSetClass::ArrayQ ArrayQ;

  typedef ProjectionEquationSetClass::SystemMatrix SystemMatrixClass;
  typedef ProjectionEquationSetClass::SystemVector SystemVectorClass;

  mpi::communicator world;

  int powerL = 3, powerH = 8;
  std::string mesher = "avro";
  std::string recon = "kexact";

  int argc = boost::unit_test::framework::master_test_suite().argc;
  char **argv = boost::unit_test::framework::master_test_suite().argv;

  if (argc == 4)
  {
    powerL = powerH = std::stoi(argv[1]);
    mesher = std::string(argv[2]);
    recon = std::string(argv[3]);
  }

  // Create a function
#if USE_TRIPLEBL

  Real a = 1.0;
  Real b = 1.0;
  Real c = 1.0;
  Real nu = 0.01;

  // Create a solution dictionary
  PyDict solnArgs;
  solnArgs[SolutionExact::ParamsType::params.a] = a;
  solnArgs[SolutionExact::ParamsType::params.b] = b;
  solnArgs[SolutionExact::ParamsType::params.c] = c;
  solnArgs[SolutionExact::ParamsType::params.nu] = nu;
  solnArgs[SolutionExact::ParamsType::params.offset] = 1;
  solnArgs[SolutionExact::ParamsType::params.scale] = -1;

  SolutionExact solnExact( solnArgs );
#else
  SolutionExact solnExact;
#endif

  FunctionNDConvertSpace<PhysD3,Real> fcnND(solnExact);

  // Nonlinear solver dicts
  PyDict SolverContinuationDict, NonlinearSolverDict, NewtonSolverDict, LinearSolverDict, LineUpdateDict;

#if defined(SANS_PETSC)
  std::cout << "Linear solver: PETSc" << std::endl;

  PyDict PreconditionerDict;
  PyDict PreconditionerILU;
  PyDict PETScDict;

  PreconditionerILU[SLA::PreconditionerASMParam::params.SubPreconditioner.Name] = SLA::PreconditionerASMParam::params.SubPreconditioner.ILU;
  PreconditionerILU[SLA::PreconditionerILUParam::params.PreconditionerSide] = SLA::PreconditionerILUParam::params.PreconditionerSide.Right;
  PreconditionerILU[SLA::PreconditionerILUParam::params.Ordering] = SLA::PreconditionerILUParam::params.Ordering.QDM;
  PreconditionerILU[SLA::PreconditionerILUParam::params.FillLevel] = 2;

  PreconditionerDict[SLA::PETScSolverParam::params.Preconditioner.Name] = SLA::PETScSolverParam::params.Preconditioner.ASM;
  PreconditionerDict[SLA::PreconditionerASMParam::params.SubPreconditioner] = PreconditionerILU;

  PETScDict[SLA::LinearSolverParam::params.LinearSolver.Solver] = SLA::LinearSolverParam::params.LinearSolver.PETSc;
  PETScDict[SLA::PETScSolverParam::params.RelativeTolerance] = 1e-9;
  PETScDict[SLA::PETScSolverParam::params.AbsoluteTolerance] = 1e-10;
  PETScDict[SLA::PETScSolverParam::params.MaxIterations] = 2000;
  PETScDict[SLA::PETScSolverParam::params.GMRES_Restart] = 300;
  PETScDict[SLA::PETScSolverParam::params.Preconditioner] = PreconditionerDict;
  PETScDict[SLA::PETScSolverParam::params.Verbose] = true;
  PETScDict[SLA::PETScSolverParam::params.computeSingularValues] = false;
  PETScDict[SLA::PETScSolverParam::params.ResidualHistoryFile] = "";
  //PETSCDict[SLA::PETScSolverParam::params.FilenameBase] = filename_base;

  LinearSolverDict[SLA::LinearSolverParam::params.LinearSolver] = PETScDict;

#elif defined(INTEL_MKL)
  std::cout << "Linear solver: MKL_PARDISO" << std::endl;
  PyDict MKL_PARDISODict;
  MKL_PARDISODict[SLA::LinearSolverParam::params.LinearSolver.Solver] = SLA::LinearSolverParam::params.LinearSolver.MKL_PARDISO;
  LinearSolverDict[SLA::LinearSolverParam::params.LinearSolver] = MKL_PARDISODict;
  NewtonSolverDict[NewtonSolverParam::params.LinearSolver] = MKL_PARDISODict;
#else
  std::cout << "Linear solver: UMFPACK" << std::endl;
  PyDict UMFPACKDict;
  UMFPACKDict[SLA::LinearSolverParam::params.LinearSolver.Solver] = SLA::LinearSolverParam::params.LinearSolver.UMFPACK;
  LinearSolverDict[SLA::LinearSolverParam::params.LinearSolver] = UMFPACKDict;
  NewtonSolverDict[NewtonSolverParam::params.LinearSolver] = UMFPACKDict;
#endif


  PyDict MesherDict;
  #if 0
  if (mesher == "refine")
  {
    MesherDict[MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.Name] = MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.refine;
    MesherDict[refineParams::params.DumpRefineDebugFiles] = false;
#ifdef SANS_REFINE
    MesherDict[refineParams::params.CallMethod] = refineParams::params.CallMethod.API;
    //MesherDict[refineParams::params.CallMethod] = refineParams::params.CallMethod.system;
#endif
  }
  else if (mesher == "EPIC")
  {
    MesherDict[MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.Name] = MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.Epic;
    std::vector<int> allBC = {0,1,2,3,4,5};
    MesherDict[EpicParams::params.SymmetricSurf] = allBC;
  }
  else if (mesher == "fefloa")
  {
    MesherDict[MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.Name] = MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.FeFloa;
  }
  else
  #endif
#ifdef SANS_AVRO
  if (mesher == "avro")
  {
    MesherDict[MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.Name] = MeshAdapterParams<PhysD3, TopoD3>::params.Mesher.avro;
    MesherDict[avroParams::params.Curved] = false;
  }
#endif
  else
    BOOST_REQUIRE_MESSAGE(false, "Unknown mesh generator.");

  typedef OutputNDConvertSpace<PhysD3, OutputCell_SolutionErrorSquared<AdvectionDiffusionTraits<PhysD3>, SolutionExact>> L2ErrorClass;
  typedef IntegrandCell_Galerkin_Output<L2ErrorClass> L2ErrorIntegrandClass;

  //Output functional
  L2ErrorClass L2ErrorOutput(solnExact);
  L2ErrorIntegrandClass L2ErrorIntegrand(L2ErrorOutput, {0});

  typedef OutputCell_SolutionSquared<AdvectionDiffusionTraits<PhysD3>> OutputClass;
  typedef OutputNDConvertSpace<PhysD3, OutputClass> NDOutputClass;
  typedef IntegrandCell_Galerkin_Output<NDOutputClass> OutputIntegrandClass;

  //Output functional
  NDOutputClass fcnOutput;
  OutputIntegrandClass outputIntegrand(fcnOutput, {0});

  //--------ADAPTATION LOOP--------

  const int maxIter = 30;

  // Adjoint equation is integrated with 2 times adjoint order for MOESS
  const int order = 1;
  const int quadOrder = -1; //2*(order+1);

  for (int power = powerL; power <= powerH; power++ )
  {
    int nk = pow(2,power);

    const int string_pad = 6;
    std::string int_pad = std::string(string_pad - std::to_string(nk).length(), '0') + std::to_string(nk) + "k";

    std::string filename_base = "tmp/L2Project_" + func + "/" + mesher + "_";

    filename_base += "CG_" + int_pad + "_" + recon + "_MultiScale/";

    boost::filesystem::path base_dir(filename_base);
    if ( not boost::filesystem::exists(base_dir) )
      boost::filesystem::create_directories(filename_base);

    std::shared_ptr<XField<PhysD3, TopoD3>> pxfld;

#if defined(SANS_AVRO)
    std::shared_ptr<avro::Context> context;
    if (mesher == "avro")
    {
      using avro::coord_t;
      using avro::index_t;

      // define the geometry and initial mesh using the avro library
      avro::Context context(3,3,2);
      context.define_geometry("box");

      XField3D_Box_Tet_X1 xfld0( world, 5, 5, 5 );

      // copy the mesh into the domain and attach the geometry
      pxfld = std::make_shared< XField_avro<PhysD3,TopoD3> >(xfld0,context);
    }
    else
#endif
    {
      pxfld = std::make_shared<XField3D_Box_Tet_X1>( world, 5, 5, 5 );
    }


    std::string output_filename = filename_base + "output.dat";
    fstream foutputhist;
    if ( world.rank() == 0 )
    {
      foutputhist.open( output_filename, fstream::out );
      BOOST_REQUIRE_MESSAGE(foutputhist.good(), "Error opening file: " + output_filename);
    }

    MesherDict[avroParams::params.FilenameBase] = filename_base;

    if (world.rank() == 0 )
    {
      // write the header to the output file
      foutputhist << "VARIABLES="
                  << std::setw(5)  << "\"Iter\""
                  << std::setw(10) << "\"DOF\""
                  << std::setw(10) << "\"Elements\""
                  << std::setw(20) << "\"u<sup>2</sup>\""
                  << std::setw(20) << "\"L<sup>2</sup> Error\""
                  << std::endl;

      foutputhist << "ZONE T=\"MultiScale " << mesher << " " << nk << "k\"" << std::endl;
    }

    int iter = 0;
    while (true)
    {
      std::cout<<"-----Adaptation Iteration " << iter << "-----" << std::endl;

      Field_CG_Cell<PhysD3, TopoD3, ArrayQ> qfld(*pxfld, order, BasisFunctionCategory_Hierarchical);
      qfld = 0;

      QuadratureOrder quadratureOrder( *pxfld, quadOrder );
      IntegrandCellFunctionClass fcnCell( fcnND, {0} );
      const std::array<Real,1> tol_projection = {{1e-15}};
      ProjectionEquationSetClass ProjEqSet(*pxfld, qfld, fcnCell, quadratureOrder, tol_projection );

      SystemVectorClass q(ProjEqSet.vectorStateSize());
      SystemVectorClass rhs(ProjEqSet.vectorEqSize());

      ProjEqSet.fillSystemVector(q);

      rhs = 0;
      ProjEqSet.residual(rhs);

      // solver
      SLA::LinearSolver< SystemMatrixClass > solver( LinearSolverDict, ProjEqSet );

      SystemVectorClass dq(ProjEqSet.vectorStateSize());
      dq = 0; // set initial guess

      solver.solve(rhs, dq);

      q -= dq;

      ProjEqSet.setSolutionField(q);

      std::string qfld_filename = filename_base + "qfld_a" + std::to_string(iter) + ".dat";
      if (iter == maxIter)
        output_Tecplot( qfld, qfld_filename );

      Real u2 = 0;
      IntegrateCellGroups<TopoD3>::integrate( FunctionalCell_Galerkin( outputIntegrand, u2 ),
                                              *pxfld, qfld,
                                              quadratureOrder.cellOrders.data(),
                                              quadratureOrder.cellOrders.size() );

      Real L2error = 0;
      IntegrateCellGroups<TopoD3>::integrate( FunctionalCell_Galerkin( L2ErrorIntegrand, L2error ),
                                              *pxfld, qfld,
                                              quadratureOrder.cellOrders.data(),
                                              quadratureOrder.cellOrders.size() );
      L2error = sqrt(L2error);

      int nDOFtotal = qfld.nDOFnative();

#ifdef SANS_MPI
      // count the number of elements possessed by this processor
      int nElem = 0;
      for (int elem = 0; elem < pxfld->nElem(); elem++ )
        if (pxfld->getCellGroupGlobal<Tet>(0).associativity(elem).rank() == world.rank())
          nElem++;

      int nElemtotal = 0;
      boost::mpi::reduce(*pxfld->comm(), nElem, nElemtotal, std::plus<int>(), 0 );
#else
      int nElemtotal = pxfld->nElem();
#endif


      if (world.rank() == 0 )
      {
        foutputhist << std::setw(5) << iter
                    << std::setw(10) << nDOFtotal
                    << std::setw(10) << nElemtotal
                    << std::setw(20) << std::setprecision(10) << std::scientific << u2
                    << std::setw(20) << std::setprecision(10) << std::scientific << L2error
                    << std::endl;
      }
      if ( iter == maxIter ) break;

      try
      {
        refine_recon_reconstructions recons;

        if (recon == "L2")
          recons = refine_recon_L2projection;
        else if (recon == "kexact")
          recons = refine_recon_Kexact;
        else
          SANS_DEVELOPER_EXCEPTION("Unknown reconstruction: %s", recon.c_str());

        int p = 2;
        Real gradation = 5;

        int targetCost = 1000*nk;

        Real nDOFperCell = 4;
        Real equiVol = sqrt(2.0)/12.0;
        Real targetComp = targetCost*equiVol/nDOFperCell;

        typedef DLA::MatrixSymS<PhysD3::D, Real> MatrixSym;

        // compute the multi-scale metric
        Field_CG_Cell<PhysD3, TopoD3, MatrixSym> metric_req(*pxfld, 1, BasisFunctionCategory_Hierarchical);
        MultiScale_metric(recons, p, gradation, targetComp, qfld, metric_req);

        // pass metric along to the mesher
        pxfld = callMesher<PhysD3, TopoD3>::call(*pxfld, metric_req, iter, MesherDict).first;

      }
      catch (...)
      {
        WriteMesh_libMeshb( *pxfld , filename_base+"/refine_debug_"+std::to_string(iter)+".meshb" );
        WriteSolution_libMeshb( qfld , filename_base+"/refine_debug_"+std::to_string(iter)+".solb" );
        throw;
      }
      iter++;
    }
  }

}

//############################################################################//
BOOST_AUTO_TEST_SUITE_END()
