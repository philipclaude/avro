//
// avro - Adaptive Voronoi Remesher
//
// Copyright 2017-2021, Philip Claude Caplan
// All rights reserved
//
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php
//
#ifndef avro_TESTING_FRAMEWORK_H_
#define avro_TESTING_FRAMEWORK_H_

#include <exception>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <fcntl.h>
#include <cstdio>
#include <math.h>

#ifdef MACHII_LIBRARY_LOCATION
#define BASE_TEST_DIR std::string(MACHII_LIBRARY_LOCATION)
#else
#define BASE_TEST_DIR std::string("library")
#endif

class TestSuite;

class TestResult
{
public:
  TestResult() :
    failed_(0),
    passed_(0),
    asserted_(0),
    exceptions_(0)
  {}

  void summary()
  {
    printf("Summary:\n");
    printf("%lu assertions passed out of %lu with %lu failures and %lu exceptions.\n",
            passed_,asserted_,failed_,exceptions_);
  }

  bool successful() const { return (failed_==0 && exceptions_==0); }
  unsigned long nb_failed() const { return failed_; }
  unsigned long nb_assert() const { return asserted_; }
  unsigned long nb_exceptions() const { return exceptions_; }

  void failed() { failed_++; }
  void passed() { passed_++; }
  void asserted() { asserted_++; }
  void exception() { exceptions_++; }

private:
  unsigned long failed_;
  unsigned long passed_;
  unsigned long asserted_;
  unsigned long exceptions_;
};

class TestCase
{
public:
  TestCase( const char* name , TestSuite* suite ) :
    name_(name), suite_(suite)
  {}

  virtual ~TestCase() {}

  virtual void run( TestResult& __result__ ) = 0;

  const std::string& name() { return name_; }

protected:
  std::string name_;
  TestSuite* suite_;
};

class TestDriver;
class TestSuite
{
public:

  TestSuite( const char* name ) :
    name_(name)
  {}

  TestSuite( const char* name , TestDriver* driver );

  const std::string& name() const { return name_; }

  void add(TestCase* t)
  {
    cases_.push_back(t);
  }

  std::size_t ntests() const { return cases_.size(); }

  int run( TestResult& __result__ )
  {
    unsigned long nb_failed0 = __result__.nb_failed();
    (void)(nb_failed0); // suppresses the warning of unused nb_failed0
    #ifndef STANDALONE

    #ifndef STDOUT_REDIRECT
    #error "must define where stdout is redirected"
    #endif

    // redirection for stdout
    int fd;
    fflush(stdout);
    fd = dup(fileno(stdout));
    FILE* x = freopen(STDOUT_REDIRECT,"a",stdout);
    (void)(x); // suppresses the warning of unused x

    // redirection for cout
    std::stringstream out;
    std::streambuf *coutbuf = std::cout.rdbuf(); // save old buf
    std::cout.rdbuf(out.rdbuf()); // redirect std::cout to out

    printf("stdout output (printf):\n");

    #endif

    printf("\n\n==== Test suite: %s ====\n\n",name_.c_str());
    printf("running %lu tests...\n",cases_.size());

    clock_t t0 = clock();
    for (std::size_t i=0;i<cases_.size();i++)
    {
      try
      {
        printf("\n--> running case %s:\n\n",cases_[i]->name().c_str());
        cases_[i]->run(__result__);
      }
      catch( std::exception& E )
      {
        printf("unexpected exception!!\n");
	std::cout << E.what() << std::endl;
        __result__.exception();
      }
    }
    clock_t t1 = clock();

    printf("\n ==== Total time = %g seconds. ====\n",double(t1-t0)/CLOCKS_PER_SEC);

    #ifndef STANDALONE

    // write the cout re-direction to the file too
    std::flush(std::cout);
    printf("\nstd::cout output:\n%s\n",out.str().c_str());

    fflush(stdout);
    dup2(fd,fileno(stdout));
    close(fd);

    // reset to standard output
    std::cout.rdbuf(coutbuf);

    if (__result__.nb_failed() > nb_failed0)
      printf("suite %s failed with %lu new errors!\n",
              name_.c_str(),__result__.nb_failed()-nb_failed0);

    #endif

    if (__result__.successful()) return 0;
    return -1;
  }

private:
  std::string name_;
  std::vector<TestCase*> cases_;
};

class TestDriver
{
public:
  void add(TestSuite* suite)
  {
    suites_.push_back(suite);
  }

  int run(TestResult& __result__)
  {
    int result = 0;
    for (std::size_t i=0;i<suites_.size();i++)
    {
      int nb_assert0 = __result__.nb_assert();
      printf("running suite: %-30s with %3lu tests ... ",
              suites_[i]->name().c_str(),suites_[i]->ntests());
      clock_t t0 = clock();
      try
      {
        if (suites_[i]->run(__result__) < 0)
          result = -1;
      }
      catch (...)
      {
        result = -1;
      }
      int nb_assert = __result__.nb_assert() -nb_assert0;
      printf(" %4d assertions ... ",nb_assert);
      clock_t t1 = clock();
      double s = double(t1-t0)/CLOCKS_PER_SEC;
      double ms = 1000*s -1000*floor(s);
      printf("done [%3d s : %-3d ms]\n",int(floor(s)),int(floor(ms)));
    }
    return result;
  }
private:
  std::vector<TestSuite*> suites_;
};
extern TestDriver* __driver__;

inline
TestSuite::TestSuite( const char* name , TestDriver* driver ) :
  TestSuite(name)
{
  driver->add(this);
}

#define UT_ASSERT_EQUALS( X , Y ) \
  do {__result__.asserted(); \
    if (X!=Y) \
    { \
      __result__.failed(); \
      printf("%s::%s: assertion %s  == %s failed on line %d of file %s.\n",suite_->name().c_str(),name_.c_str(),#X,#Y,__LINE__,__FILE__); \
      std::cout << "--> expected " << Y << " but received " << X << std::endl; \
    } \
    else \
      __result__.passed(); \
  }while(0)

#define UT_ASSERT(X) \
  do {__result__.asserted(); \
    if (!(X)) \
    { \
      __result__.failed(); \
      printf("%s::%s: assertion %s failed on line %d of file %s.\n",suite_->name().c_str(),name_.c_str(),#X,__LINE__,__FILE__); \
    } \
    else \
      __result__.passed(); \
 }while(0)

#define UT_ASSERT_NEAR( X , Y , Z ) \
  do {__result__.asserted(); \
    if ( sqrt((X-Y)*(X-Y))>Z ) \
    { \
      __result__.failed(); \
      printf("%s::%s: assertion %s (%g) == %s +/- %s failed on line %d of file %s.\n",suite_->name().c_str(),__FUNCTION__,#X,X,#Y,#Z,__LINE__,__FILE__); \
    } \
    else \
      __result__.passed(); \
  }while(0)

#define UT_ASSERT_CLOSE( X , Y , Z , Z0 ) UT_ASSERT_NEAR(X,Y,Z)
#define UT_ASSERT_SMALL( X , Z ) UT_ASSERT_NEAR(X,0,Z)

#define UT_CATCH_EXCEPTION(X) \
  do {__result__.asserted(); \
    try { (X); __result__.failed(); \
      printf("%s::%s: expected exception on line %d of file %s.\n",suite_->name().c_str(),name_.c_str(),__LINE__,__FILE__); \
    } catch (...) { __result__.passed(); break; } \
 }while(0)

#ifdef STANDALONE
#define UT_TEST_SUITE(X) namespace suite_##X { TestSuite __suite__(#X);
#else
#define UT_TEST_SUITE(X) namespace suite_##X { static TestSuite __suite__(#X,__driver__);
#endif

#define UT_TEST_CASE(X) class X : public TestCase { public: X() : TestCase(#X,&__suite__) {suite_->add(this); } \
void run(TestResult& __result__)

#define UT_TEST_CASE_SKIP(X) class X : public TestCase { public: X() : TestCase(#X,&__suite__) {suite_->add(this); } \
void run(TestResult& __result__) {printf("skipping test!\n");} void run_skip(TestResult& __result__)

#define UT_TEST_CASE_END(X) }; static X instance_##X;

#ifdef STANDALONE
#define UT_TEST_SUITE_END(X) } void ut_pre(int,char**); void ut_post(); int main(int argc, char* argv[]) \
    { ut_pre(argc,argv); TestResult __result__; suite_##X::__suite__.run(__result__); __result__.summary(); ut_post(); if (__result__.nb_failed()!=0 || __result__.nb_exceptions()!=0) return 1; return 0; }
#else
#define UT_TEST_SUITE_END(X) }
#endif

#endif
