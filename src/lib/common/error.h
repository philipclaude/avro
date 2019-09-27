#ifndef URSA_COMMON_ERROR_H_
#define URSA_COMMON_ERROR_H_

#include <cstdio>
#include <exception>
#include <string>

namespace ursa
{
#ifndef unlikely
#ifdef __GNUC__
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define unlikely(x) x
#endif
#endif

#define ursa_throw(args...) throw(ursa::Exception(__FILE__,__LINE__,args));

#if (1 || defined(URSA_DEBUG))
#define ursa_assert(X) if(unlikely(!(X))) { printf("\nfailed to assert %s in file %s line %d\n",(#X),__FILE__,__LINE__);ursa_throw("assertion error");}
#define ursa_assert_msg(X,...) { try{ ursa_assert(X); } catch(...) { printf(__VA_ARGS__); ursa_throw("assertion error"); } }
#else
#define ursa_assert(X) {}
#define ursa_assert_msg(X,...) {}
#endif

#define ursa_implement ursa_throw("not implemented");

#define ursa_assert_not_reached {printf("\nthis should not have been reached!\n"); ursa_assert(false);}

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
