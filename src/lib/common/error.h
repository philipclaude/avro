#ifndef luma_COMMON_ERROR_H_
#define luma_COMMON_ERROR_H_

#include <cstdio>
#include <exception>
#include <string>

namespace luma
{
#ifndef unlikely
#ifdef __GNUC__
#define unlikely(x) __builtin_expect(!!(x), 0)
#else
#define unlikely(x) x
#endif
#endif

#define luma_throw(args...) throw(luma::Exception(__FILE__,__LINE__,args));

#if (1 || defined(luma_DEBUG))
#define luma_assert(X) if(unlikely(!(X))) { printf("\nfailed to assert %s in file %s line %d\n",(#X),__FILE__,__LINE__);luma_throw("assertion error");}
#define luma_assert_msg(X,...) { try{ luma_assert(X); } catch(...) { printf(__VA_ARGS__); luma_throw("assertion error"); } }
#define luma_assert_else(X,Y) { try{ luma_assert(X); } catch(...) { Y ; luma_throw("assertion error"); } }
#else
#define luma_assert(X) {}
#define luma_assert_msg(X,...) {}
#endif

#define luma_implement luma_throw("not implemented");

#define luma_assert_not_reached {printf("\nthis should not have been reached!\n"); luma_assert(false);}

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
