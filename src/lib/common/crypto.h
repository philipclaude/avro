//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef WEBGLPP_CRYPTO_H_
#define WEBGLPP_CRYPTO_H_

#include <cstdint>
#include <cstdlib>

int b64_encode_string(const char *in, int in_len, char *out, int out_size);
unsigned char* SHA1(const unsigned char *d, size_t n, unsigned char *md);

#endif
