//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2020, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_COMMON_ERROR_H_
#define avro_COMMON_ERROR_H_

#include <cstdio>
#include <exception>
#include <string>

namespace avro
{

#ifndef unlikely
#ifdef __GNUC__
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define unlikely(x) x
#endif
#endif

#define avro_throw(args...) throw(avro::Exception(__FILE__,__LINE__,args));

#if (1 || defined(avro_DEBUG))
  //#if AVRO_MPI
  //  #define avro_assert(X) if(unlikely(!(X))) { printf("\nfailed to assert %s in file %s line %d\n",(#X),__FILE__,__LINE__); avro_throw("assertion error");}
  //#else
    #define avro_assert(X) if(unlikely(!(X))) { printf("\nfailed to assert %s in file %s line %d\n",(#X),__FILE__,__LINE__);avro_throw("assertion error");}
  //#endif

  #define avro_assert_msg(X,...) { try{ avro_assert(X); } catch(...) { printf(__VA_ARGS__); avro_throw("assertion error"); } }
  #define avro_assert_else(X,Y) { try{ avro_assert(X); } catch(...) { Y ; avro_throw("assertion error"); } }
#else
  #define avro_assert(X) {}
  #define avro_assert_msg(X,...) {}
#endif

#define avro_implement avro_throw("not implemented");

#define avro_assert_not_reached {printf("\nthis should not have been reached!\n"); avro_assert(false);}

void call_backtrace(const int start=1,const int end=2);

class Exception : public std::exception {

  public:

    explicit Exception(const char *file, int line,const char *fmt, ...);
    virtual ~Exception() throw() {}
    virtual const char* what() const throw()
      { return message.c_str(); }

  protected:
    std::string message;
};

}

#endif
