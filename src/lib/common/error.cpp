// avro: Adaptive Voronoi Remesher
// Copyright 2017-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
#include "common/error.h"

#if !defined(__CYGWIN__) && !defined(__MINGW32__) && !defined(_MSC_VER)

#include <execinfo.h>
#include <cxxabi.h>
#include <dlfcn.h>
#include <sys/types.h>
#include <unistd.h>
#include <cstdlib>
#include <cstdarg>
#include <csignal>

namespace avro
{

Exception::Exception(const char *file, int line,const char *fmt, ...) {
  char buffer[512];
  char info[512];
  va_list argp;

  va_start(argp,fmt);
  vsprintf(buffer,fmt,argp);
  va_end(argp);

  sprintf(info," (file %s, line %d) ===\n",file,line);

  message = "\n=== caught avro::Exception";
  message += std::string(info);
  message += std::string(buffer);
  printf("%s\n\n=== Backtrace ===\n",message.c_str());
  call_backtrace(1);
}

inline void
call_backtrace(const int start,const int end) {
  int 	i;
  enum 	{MAX_DEPTH=10};
  void 	*trace[MAX_DEPTH];
  char 	*demangled;
  int 	trace_size,status=0;
  Dl_info dlinfo;
  const char *symname;

  trace_size = backtrace(trace,MAX_DEPTH);

  for (i=start;i<trace_size-end;i++)
  {
    if (!dladdr(trace[i],&dlinfo)) continue;
    symname = dlinfo.dli_sname;

    demangled = abi::__cxa_demangle(symname, NULL, 0, &status);
    if(status==0 && demangled)
      symname = demangled;

    if (symname)
      printf("%s\n",symname);

    if (demangled) {
      free(demangled);
    }
  }
	printf("\n");
}

} // avro

#endif // __CYGWIN__
