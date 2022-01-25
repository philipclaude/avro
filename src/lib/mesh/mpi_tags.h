//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_MESH_MPI_TAGS_H_
#define AVRO_LIB_MESH_MPI_TAGS_H_

// point tags
#define TAG_COORDINATE 0
#define TAG_PARAMETER  1
#define TAG_GEOMETRY   2
#define TAG_GEOMETRY_NUMBERS 3
#define TAG_LOCAL2GLOBAL 4
#define TAG_FIXED 5
#define TAG_AGE 6

// topology tags
#define TAG_CELL_INDEX 11
#define TAG_CELL_FIRST 12
#define TAG_CELL_LAST  13

// tag for passing miscellaneous data
#define TAG_MISC 200

#endif
