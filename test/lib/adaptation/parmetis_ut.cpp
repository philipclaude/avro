//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "unit_tester.hpp"

#include "common/mpi.hpp"
#include "common/tools.h"
#include "avro_types.h"

#if AVRO_MPI

#include <parmetis.h>
#if PARMETIS_MAJOR_VERSION == 3
#define PARM_INT idxtype
#define PARM_REAL float
#elif PARMETIS_MAJOR_VERSION == 4
#define PARM_INT idx_t
#define PARM_REAL float
#else
#error "unknown version of parmetis"
#endif

#endif

using namespace avro;

UT_TEST_SUITE( parmetis_test_suite )

#if AVRO_MPI

#if TEST_NUM_PROCS == 3

UT_TEST_CASE( test1 )
{
  index_t rank = mpi::rank();

  std::vector<PARM_INT> xadj;
  std::vector<PARM_INT> adjncy;
  std::vector<PARM_INT> egwt;
  std::vector<PARM_INT> vtxdist = {0,5,10,15};
  if (rank == 0)
  {
    xadj   = {0,2,5,8,11,13};
    adjncy = {1,5,0,2,6,1,3,7,2,4,8,3,9};
    egwt   = {1,1,1,1,1,1,1,1,1,1,1,1,1};
  }
  else if (rank == 1)
  {
    xadj   = {0,3,7,11,15,18};
    adjncy = {0,6,10,1,5,7,11,2,6,8,12,3,7,9,13,4,8,14};
    egwt   = {1,1, 1,1,1,1, 1,1,1,1, 1,1,1,1, 1,1,1, 1};
  }
  else if (rank == 2)
  {
    xadj   = {0,2,5,8,11,13};
    adjncy = {5,11,6,10,12,7,11,13,8,12,14,9,13};
    egwt   = {1, 1,1, 1, 1,1, 1, 1,1, 1, 1,1, 1};
  }
  else
    avro_assert_not_reached;

    // setup the parmetis version of the adjacency graph
    PARM_INT *pvwgt = NULL;
    PARM_INT *padjwgt = egwt.data();//NULL;
    PARM_INT wgtflag = 1; // 0 = no weights, 1 = edge weights, 2 = vertex weights, 3 = both
    PARM_INT edgecut = 0;
    PARM_INT bias = 0;
    PARM_INT ncon = 1;
    PARM_INT nparts = mpi::size();

    std::vector<PARM_REAL> tpwgts(ncon*nparts,1./nparts);
    std::vector<PARM_REAL> ubvec(ncon,1.05);
    PARM_INT options[METIS_NOPTIONS];
    /*options[0] = 1;
    options[1] = 3;
    options[2] = 0;
    options[3] = PARMETIS_PSR_COUPLED;*/

    //METIS_SetDefaultOptions[METIS_NOPTIONS];
    //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    //options[METIS_OPTION_DBGLVL] = METIS_DBG_INFO;



    mpi::barrier();

    index_t nb_vert = 5;

    // partition the graph!
    std::vector<PARM_INT> part(nb_vert,rank);
    MPI_Comm comm = MPI_COMM_WORLD;
    #if 1
    int result = ParMETIS_V3_PartKway( vtxdist.data() , xadj.data() , adjncy.data() ,
                            pvwgt, padjwgt, &wgtflag,
                            &bias , &ncon , &nparts,
                            tpwgts.data() , ubvec.data(),
                            options,
                            &edgecut,
                            part.data(),
                            &comm);
    #else
    PARM_REAL itr = 1000.;
    index_t mem_vertex = 1*sizeof(index_t);
    std::vector<PARM_INT> vsize(nb_vert,mem_vertex);
    int result = ParMETIS_V3_AdaptiveRepart( vtxdist.data() , xadj.data() , adjncy.data() ,
                            pvwgt, vsize.data(), padjwgt, &wgtflag,
                            &bias , &ncon , &nparts,
                            tpwgts.data() , ubvec.data(), &itr,
                            options,
                            &edgecut,
                            part.data(),
                            &comm);
    #endif
    printf("edge cut = %d\n",edgecut);
    mpi::barrier();

    print_inline(part);

    UT_ASSERT_EQUALS( result , METIS_OK );

}
UT_TEST_CASE_END(test1)

#endif // TEST_NUM_PROCS

#endif

UT_TEST_SUITE_END( parmetis_test_suite )
