// luma: Adaptive Voronoi Remesher
// Copyright 2017-2018, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef NLOHMANN_JSON_FWD_HPP
#define NLOHMANN_JSON_FWD_HPP

#include <cstdint> // int64_t, uint64_t
#include <map> // map
#include <memory> // allocator
#include <string> // string
#include <vector> // vector

/*!
@brief namespace for Niels Lohmann
@see https://github.com/nlohmann
@since version 1.0.0
*/
namespace nlohmann
{
/*!
@brief default JSONSerializer template argument
This serializer ignores the template arguments and uses ADL
([argument-dependent lookup](http://en.cppreference.com/w/cpp/language/adl))
for serialization.
*/
template<typename = void, typename = void>
struct adl_serializer;

template<template<typename U, typename V, typename... Args> class ObjectType =
         std::map,
         template<typename U, typename... Args> class ArrayType = std::vector,
         class StringType = std::string, class BooleanType = bool,
         class NumberIntegerType = std::int64_t,
         class NumberUnsignedType = std::uint64_t,
         class NumberFloatType = double,
         template<typename U> class AllocatorType = std::allocator,
         template<typename T, typename SFINAE = void> class JSONSerializer =
         adl_serializer>
class basic_json;

/*!
@brief JSON Pointer
A JSON pointer defines a string syntax for identifying a specific value
within a JSON document. It can be used with functions `at` and
`operator[]`. Furthermore, JSON pointers are the base for JSON patches.
@sa [RFC 6901](https://tools.ietf.org/html/rfc6901)
@since version 2.0.0
*/
template<typename BasicJsonType>
class json_pointer;

/*!
@brief default JSON class
This type is the default specialization of the @ref basic_json class which
uses the standard template types.
@since version 1.0.0
*/
using json = basic_json<>;
}

#endif
