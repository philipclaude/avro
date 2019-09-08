// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#ifndef TYPE2STRING_H
#define TYPE2STRING_H

#include <string>

namespace boost
{
namespace python
{
// Forward declare list here to reduce pre-compiled header size
class list;
namespace api
{
class object;
}
using api::object;
}
}

#include <vector>

namespace SANS
{

struct PyDict;
class DictKeyPair;

template<class T>
struct Type2String;

template<> struct Type2String<int>                   { static std::string str() { return "int";    } };
template<> struct Type2String<float>                 { static std::string str() { return "float";  } };
template<> struct Type2String<double>                { static std::string str() { return "double"; } };
template<> struct Type2String<bool>                  { static std::string str() { return "bool";   } };
template<> struct Type2String<std::string>           { static std::string str() { return "string"; } };
template<> struct Type2String<PyDict>                { static std::string str() { return "dict";   } };
template<> struct Type2String<boost::python::list>   { static std::string str() { return "list";   } };

template<> struct Type2String<std::vector<int>>      { static std::string str() { return "list of int";   } };
template<> struct Type2String<std::vector<double>>   { static std::string str() { return "list of double";   } };
template<> struct Type2String<DictKeyPair>           { static std::string str() { return "dict";   } };

template<> struct Type2String<boost::python::object> { static std::string str() { return "python object"; } };

} //namespace SANS

#endif //TYPE2STRING_H
