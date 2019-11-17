#ifndef LUNA_COMMON_ERROR_H_
#define LUNA_COMMON_ERROR_H_

#include <cstdio>
#include <exception>
#include <string>

namespace luna
{
#ifndef unlikely
#ifdef __GNUC__
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define unlikely(x) x
#endif
#endif

#define luna_throw(args...) throw(luna::Exception(__FILE__,__LINE__,args));

#if (1 || defined(LUNA_DEBUG))
#define luna_assert(X) if(unlikely(!(X))) { printf("\nfailed to assert %s in file %s line %d\n",(#X),__FILE__,__LINE__);luna_throw("assertion error");}
#define luna_assert_msg(X,...) { try{ luna_assert(X); } catch(...) { printf(__VA_ARGS__); luna_throw("assertion error"); } }
#else
#define luna_assert(X) {}
#define luna_assert_msg(X,...) {}
#endif

#define luna_implement luna_throw("not implemented");

#define luna_assert_not_reached {printf("\nthis should not have been reached!\n"); luna_assert(false);}

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
