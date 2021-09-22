//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//

// no include guards because it should only be included once

Sign avro_side1_nd_exact_pck(const double* p0, const double* p1, const double* q0 ,unsigned short dim );

Sign avro_side2_nd_exact_pck( const double* p0, const double* p1, const double* p2,
                             const double* q0, const double* q1 ,unsigned short dim );

Sign avro_side3_nd_exact_pck( const double* p0,const double* p1,const double* p2,const double* p3,
                             const double* q0,const double* q1,const double* q2 ,unsigned short dim );

Sign avro_side4_nd_exact_pck( const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,
                             const double* q0,const double* q1,const double* q2,const double* q3 ,unsigned short dim );

Sign avro_side4_3d_exact_pck( const double* p0,const double* p1,const double* p2,const double* p3,const double* p4 , bool sos );

Sign avro_side5_nd_exact_pck( const double* p0,const double* p1,const double* p2,const double* p3,const double* p4,const double* p5,
                             const double* q0,const double* q1,const double* q2,const double* q3,const double* q4 ,unsigned short dim );
