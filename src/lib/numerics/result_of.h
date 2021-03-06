//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef AVRO_LIB_NUMERICS_RESULT_OF_H_
#define AVRO_LIB_NUMERICS_RESULT_OF_H_

namespace avro {

template<typename S,typename T> class result_of;

template<int N> class result_of<real_t,SurrealS<N>> { public: typedef SurrealS<N> type; };
template<int N> class result_of<SurrealS<N>,real_t> { public: typedef SurrealS<N> type; };
template<int N> class result_of<SurrealS<N>,SurrealS<N>>  { public: typedef SurrealS<N> type; };
template<> class result_of<real_t,real_t> { public: typedef real_t type; };
template<> class result_of<float,float> { public: typedef float type; };

} // avro

#endif
