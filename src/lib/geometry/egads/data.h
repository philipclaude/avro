//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_LIB_GEOMETRY_EGADS_DATA_H_
#define avro_LIB_GEOMETRY_EGADS_DATA_H_

#include "common/error.h"
#include "avro_types.h"

#include "geometry/egads/egads_fwd.h"

#include <egads.h>

#include <string>

namespace avro
{

namespace EGADS
{

std::string get_error_message( const int ierror );

#define EGADS_CHECK_SUCCESS(X) { int ierror = (X); if (ierror!=0) printf("egads error: returned %d (%s) at %s (line %d)\n",ierror,get_error_message(ierror).c_str(),__FILE__,__LINE__); }
#define EGADS_ENSURE_SUCCESS(X) { int ierror = (X); if (ierror!=0) avro_throw(EGADS::get_error_message(ierror).c_str()); }

typedef struct
{
  int object_class, member_type;
  int nb_children, *senses;
  ego reference, *children;
  ego previous,next;
  double data[4];
} egoData;

namespace utilities
{

inline bool
object_tessellatable(int object_class, int member_type)
{
	// degenerate edges and loops don't directly have tessellations
	if (object_class==EDGE && member_type==DEGENERATE) return false;
	if (object_class==LOOP) return false;
	if (object_class==BODY) return false;
	if (object_class==SHELL) return false;
	if (object_class==EDGE || object_class==FACE)
		return true;
	if (object_class==NODE) return true;
	return false;
}

inline coord_t
topological_number(int object_class,int member_type)
{
  avro_assert_msg( object_class>5 , "object class = %d" , object_class );
  avro_assert( object_class>5 );
  avro_assert( member_type>-2 );

  coord_t number = 0;
  if (object_class==20)
  {
    number = 0;
  }
  else if (object_class==10 || object_class==11 || object_class==21 || object_class==22)
  {
    number = 1;
  }
  else if (object_class==12 || object_class==23 || object_class==24)
  {
    number = 2;
  }
  else if (object_class==25)
  {
    if (member_type==6) number = 1; // WIREBODY
    else if (member_type>=7 && member_type<=9) number = 2;
    else printf("unrecognized mtype %d\n",member_type);
  }
  else
  {
    printf("unrecognized topological number for oclass = %d\n",object_class);
    avro_assert(false);
  }
  return number;
}

inline const std::string
object_class_name(int object_class)
{

	if (object_class==0) return "CONTEXT";
	if (object_class==1) return "TRANSFORM";
	if (object_class==2) return "TESSELLATION";
	if (object_class==3) return "NIL";
	if (object_class==4) return "EMPTY";
	if (object_class==5) return "REFERENCE";

	if (object_class==10) return "PCURVE";
	if (object_class==11) return "CURVE";
	if (object_class==12) return "SURFACE";

	if (object_class==20) return "NODE";
	if (object_class==21) return "EDGE";
	if (object_class==22) return "LOOP";
	if (object_class==23) return "FACE";
	if (object_class==24) return "SHELL";
	if (object_class==25) return "BODY";
	if (object_class==26) return "MODEL";

  if (object_class==27) return "TESSERACT_NODE";
  if (object_class==28) return "TESSERACT_EDGE";
  if (object_class==29) return "TESSERACT_FACE";
  if (object_class==30) return "TESSERACT_CUBE";

	return "placeholder";
}

inline const std::string
member_type_name(int object_class,int member_type)
{

	if ( object_class_name(object_class)=="PCURVE" ||
       object_class_name(object_class)=="CURVE" )
	{
		if (member_type==1) return "LINE";
		if (member_type==2) return "CIRCLE";
		if (member_type==3) return "ELLIPSE";
		if (member_type==4) return "PARABOLA";
		if (member_type==5) return "HYPERBOLA";
		if (member_type==6) return "TRIMMED";
		if (member_type==7) return "BEZIER";
		if (member_type==8) return "BSPLINE";
		if (member_type==9) return "OFFSET";
	}
	else if (object_class_name(object_class)=="SURFACE")
	{
		if (member_type==1) return "PLANE";
		if (member_type==2) return "SPHERICAL";
		if (member_type==3) return "CYLINDRICAL";
		if (member_type==4) return "REVOLUTION";
		if (member_type==5) return "TOROIDAL";
		if (member_type==10) return "CONICAL";
		if (member_type==11) return "EXTRUSION";
	}
	else if ( object_class>=20 && object_class<=26 )
	{
		if (member_type==-1) return "SREVERSE";
		if (member_type==0) return "NOMTYPE";
    if (object_class!=21 && member_type==1) return "SFORWARD";
		if (member_type==1) return "ONENODE";
		if (member_type==2) return "TWONODE";
		if (member_type==3) return "OPEN";
		if (member_type==4) return "CLOSED";
		if (member_type==5) return "DEGENERATE";
		if (member_type==6) return "WIREBODY";
		if (member_type==7) return "FACEBODY";
		if (member_type==8) return "SHEETBODY";
		if (member_type==9) return "SOLIDBODY";
	}
  else if ( object_class>=27 && object_class<=30 )
  {
    return "NOMTYPE";
  }
	return "placeholder";
}

} // utilities

} // EGADS

} // avro

#endif
