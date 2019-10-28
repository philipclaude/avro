#ifndef URSA_COMMON_TOOLS_H_
#define URSA_COMMON_TOOLS_H_

#include "common/error.h"
#include "common/types.h"

#include <algorithm>
#include <cstdio>
#include <sstream>
#include <stdarg.h>
#include <string>
#include <vector>

namespace ursa
{

#define UNUSED(x) (void)(x);

real_t randomWithin( const real_t lo , const real_t hi );
int randomWithin( const int lo , const int hi );

template <typename T>
inline std::string
stringify( const T& x )
{
  std::ostringstream ss;
  ss << x;
  return ss.str();
}

template <typename type> inline type unstringify(const std::string &txt);

template <>
inline int
unstringify<int>(const std::string &txt)
{
  return atoi(txt.c_str());
}

template <>
inline double
unstringify<double>(const std::string &txt)
{
  return atof(txt.c_str());
}

template <>
inline std::string
unstringify<std::string>(const std::string &txt)
{
  return txt;
}

template<>
inline index_t
unstringify<index_t>( const std::string &txt)
{
  return (index_t) atoi(txt.c_str());
}

template<>
inline coord_t
unstringify<coord_t>( const std::string &txt)
{
  return (coord_t) atoi(txt.c_str());
}

template<>
inline bool
unstringify<bool>( const std::string &txt)
{
  if (txt=="true" || txt=="True") return true;
  return false;
}

inline std::vector<std::string>
split(const std::string &txt,const std::string &delim)
{
  std::vector<std::string> strings;
  std::string::const_iterator start,end;
  if (delim.empty()) {
    strings.push_back(txt);
    return strings;
  }

  start = txt.begin();
  while (true) {
    end = search(start,txt.end(),delim.begin(),delim.end());
    std::string s(start,end);
    if (s.size()>0)
      strings.push_back(s);
    if (end==txt.end())
      break;
    start = end +delim.size();
  }
  return strings;
}

inline std::string
lowercase( const std::string& s )
{
  std::string S = s;
  std::transform( S.begin() , S.end() , S.begin() , ::tolower );
  return S;
}

template<typename type>
index_t
indexof( const type* x , const std::vector<type*>& X )
{
	for (index_t k=0;k<X.size();k++)
	{
		if (X[k]==x) return k;
	}
	ursa_assert_not_reached;
	return 0; // to avoid compiler warnings
}

template <typename T>
void
uniquify( std::vector<T> &x )
{
  // uniquify the list of indices
  std::sort( x.begin() , x.end() );
  typename std::vector<T>::iterator it = std::unique( x.begin() , x.end() );
  x.erase(it,x.end());
}

template<typename T>
std::vector<T>
linspace( T lo , T hi , int n=-1 )
{
  if (n==-1) n = hi -lo +1;
  T dt = T(hi -lo)/T(n -1);
  std::vector<T> y(n);
  for (index_t k=0;k<index_t(n);k++)
    y[k] = lo +k*dt;
  return y;
}

inline
std::vector<index_t>
linspace( index_t n )
{
	std::vector<index_t> y(n,0);
	for (index_t k=1;k<n;k++)
		y[k] = y[k-1] +1;
	return y;
}

template<typename T>
std::string
unique_label( T& x0 , T& x1 )
{
  ursa_assert(x0!=x1);
  if (x0>x1) std::swap(x0,x1);
  std::string s = stringify(x0)+"|"+stringify(x1);
  return s;
}

template<typename T>
std::string
unique_label( std::vector<T>& x )
{
  ursa_assert( x.size()>0 );
  std::sort( x.begin() , x.end() );
  std::string s = stringify(x[0]);
  for (index_t k=1;k<x.size();k++)
    s += "|"+stringify(x[k]);
  return s;
}

template<typename T>
std::string
unique_label_skip( T* x , index_t n , index_t kskip )
{
  std::vector<index_t> f(x,x+n);
  f.erase( f.begin() +kskip );
  return unique_label(f);
}

inline void
printInfo()
{

  printf("\nursa compiled with ");

  if (__cplusplus == 201103L) printf("c++11");
  else if (__cplusplus == 199711L ) printf("c++98");
  else printf("pre-standard c++");

  printf("deleting memory upon termination using");
  #ifdef URSA_SMART_PTR
    printf(" internal smart pointers.\n");
  #else
    printf(" std::shared_ptr (c++11).\n");
  #endif
  printf("\n");
}

void ursaPrintf( const char* fmt , ... );

template<typename type>
static void
printValue( const type& x );

template<>
inline void
printValue( const index_t& x )
{
  printf("%d ",int(x));
}

template<>
inline void
printValue( const int& x )
{
  printf("%d ",int(x));
}


template<>
inline void
printValue( const real_t& x )
{
  printf("%g ",x);
}

template<>
inline void
printValue( const float& x )
{
  printf("%g ",x);
}

template<>
inline void
printValue( const std::vector<real_t>& x )
{
  printf(" [ ");
  for (index_t j=0;j<x.size();j++)
    printf("%g ",x[j]);
  printf("] ");
}

template<typename type>
static void
printInline( const std::vector<type>& s , const std::string& name=std::string() , const int id=-1 , const index_t nt=0 )
{
  for (index_t k=0;k<nt;k++) printf("  ");
  if (name.size()>0)
      printf("%s",name.c_str());
  if (id>=0)
      printf("[%d]: ",id);
  printf("( ");
  for (index_t j=0;j<s.size();j++)
      printValue(s[j]);
  printf(")\n");
}

inline void
tabbedPrint( const index_t nt , const char *fmt , ... )
{
  for (index_t k=0;k<nt;k++) printf("  ");
  va_list args;
  va_start(args,fmt);
  char buffer[256];
  vsprintf(buffer,fmt,args);
  va_end(args);
}

inline void
trim_string( std::string& str )
{
  const char* white_space = " \t\n\r";
  size_t location;
  location = str.find_first_not_of(white_space);
  str.erase(0,location);
  location = str.find_last_not_of(white_space);
  str.erase(location + 1);
}

template<typename type>
struct SortBy
{
  const std::vector<type>& target_;
  explicit SortBy(const std::vector<type>& _target) : target_(_target) {}
  bool operator()(index_t a, index_t b) const
  {
    return target_[a] < target_[b];
  }
};

template void printInline( const std::vector<index_t>& s , const std::string& name , const int id , const index_t nt );
template void printInline( const std::vector<double>& s , const std::string& name , const int id , const index_t nt );
template void printInline( const std::vector<float>& s , const std::string& name , const int id , const index_t nt );

}

#endif
