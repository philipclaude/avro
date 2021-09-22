//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#include "geometry/egads/data.h"

#include <egads.h>

namespace avro
{

namespace EGADS
{

std::string
get_error_message( const int ierror )
{
	if (ierror==-32) return "EGADS_READERR";
	if (ierror==-31) return "EGADS_TESSTATE";
	if (ierror==-30) return "EGADS_EXISTS";
	if (ierror==-29) return "EGADS_ATTRERR";
	if (ierror==-28) return "EGADS_TOPOCNT";
	if (ierror==-27) return "EGADS_OCSEGFLT";
	if (ierror==-26) return "EGADS_BADSCALE";
	if (ierror==-25) return "EGADS_NOTORTHO";
	if (ierror==-24) return "EGADS_DEGEN";
	if (ierror==-23) return "EGADS_CONSTERR";
	if (ierror==-22) return "EGADS_TOPOERR";
	if (ierror==-21) return "EGADS_GEOMERR";
	if (ierror==-20) return "EGADS_NOTBODY";
	if (ierror==-19) return "EGADS_WRITERR";
	if (ierror==-18) return "EGADS_NOTMODEL";
	if (ierror==-17) return "EGADS_NOLOAD";
	if (ierror==-16) return "EGADS_RANGERR";
	if (ierror==-15) return "EGADS_NOTGEOM";
	if (ierror==-14) return "EGADS_NOTTESS";
	if (ierror==-13) return "EGADS_EMPTY";
	if (ierror==-12) return "EGADS_NOTTOPO";
	if (ierror==-11) return "EGADS_REFERCE";
	if (ierror==-10) return "EGADS_NOTXFORM";
	if (ierror==-9)  return "EGADS_NOTCNTX";
	if (ierror==-8)  return "EGADS_MIXCNTX";
	if (ierror==-7)  return "EGADS_NODATA";
	if (ierror==-6)  return "EGADS_NONAME";
	if (ierror==-5)  return "EGADS_INDEXERR";
	if (ierror==-4)  return "EGADS_MALLOC";
	if (ierror==-3)  return "EGADS_NOTOBJ";
	if (ierror==-2)  return "EGADS_NULLOBJ";
	if (ierror==-1)  return "EGADS_NOTFOUND";
	if (ierror==0)   return "EGADS_SUCCESS";
	if (ierror==1)   return "EGADS_OUTSIDE";

	// in case the API changes and this function needs to be updated
	return "EGADS_UNKNOWN_ERROR";
}

} // EGADS

} // avro
