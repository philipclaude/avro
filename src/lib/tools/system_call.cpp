// Solution Adaptive Numerical Simulator (SANS)
// Copyright 2013-2019, Massachusetts Institute of Technology
// Licensed under The GNU Lesser General Public License, version 2.1
// See http://www.opensource.org/licenses/lgpl-2.1.php

#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <cstdlib> // system
#include <sys/wait.h>
#ifdef SANS_MPI
#include <unistd.h> // gethostname
#include <limits.h> // HOST_NAME_MAX
#ifndef HOST_NAME_MAX
#  if defined(_POSIX_HOST_NAME_MAX)
#    define HOST_NAME_MAX _POSIX_HOST_NAME_MAX
#  elif defined(_SC_HOST_NAME_MAX)
#    define HOST_NAME_MAX _SC_HOST_NAME_MAX
#  else
#    define HOST_NAME_MAX 255
#  endif
#endif // HOST_NAME_MAX
#endif // SANS_MPI

#include "system_call.h"

#include "SANSException.h"
#include "timer.h"

namespace numpack 
{

void system_call(const std::string& name, const std::string& exec, const std::string& args)
{
#ifdef SANS_MPI
  char hostname[HOST_NAME_MAX];
  gethostname(hostname, HOST_NAME_MAX);
  std::cout << std::endl;
  std::cout << "Executing " << name << " on host: " << hostname << std::endl;
  std::cout << std::endl;
#endif

  // Check if command processor exists
  int maxTry = 30;
  int iter = 0;
  while (std::system(NULL) == 0 && iter < maxTry)
  {
    std::cout << "Waiting for access to command processor..." << std::endl;
    std::this_thread::sleep_for( std::chrono::milliseconds(5000) );
    iter++;
  }
  if (iter == maxTry)
    SANS_RUNTIME_EXCEPTION("Could not access command processor!");

  std::string command = exec + " " + args;

  timer clock;

  int ret = 0;
  iter = 0;
  do
  {
    // wait first for slow hard-drives, and subsequently for system call problems
    std::this_thread::sleep_for( std::chrono::milliseconds(1000) );

    std::cout << "Calling " << name << " with the command: \n" << command << std::endl;
    ret = system(command.c_str());
    std::cout << name << " returned with status " << WEXITSTATUS(ret) << std::endl;
    iter++;
  }
  while (WEXITSTATUS(ret) != 0 && iter < maxTry);

  std::cout << name << " execution time : " << clock.elapsed() << " s" << std::endl;

  //
  // Details about exit codes can be found here:
  //
  // https://linux.die.net/man/3/system
  // https://linux.die.net/man/2/wait
  // http://www.tldp.org/LDP/abs/html/exitcodes.html
  //
  if (WEXITSTATUS(ret) == 127)
    SANS_RUNTIME_EXCEPTION("Command '%s' not found. Is it in PATH?.", exec.c_str());

  else if (WEXITSTATUS(ret) == 126)
    SANS_RUNTIME_EXCEPTION("Command '%s' does not have execute permissions.", exec.c_str());

  else if (WIFSIGNALED(ret) && (WTERMSIG(ret) == SIGINT || WTERMSIG(ret) == SIGQUIT))
    SANS_RUNTIME_EXCEPTION("%s was interrupted or quit unexpectedly.", name.c_str());

  else if (WEXITSTATUS(ret) != 0)
    SANS_RUNTIME_EXCEPTION("%s returned with error %d.", exec.c_str(), WEXITSTATUS(ret));
}

// useful function to make sure that files are available to read
void wait_for_file(const std::string& filename)
{
  int maxTry = 30;
  int iter = 0;

  // wait first for slow harddrives
  std::this_thread::sleep_for( std::chrono::milliseconds(1000) );

  std::ifstream file(filename.c_str());
  while (!file.is_open() && iter < maxTry)
  {
    file.close();
    std::cout << "Waiting for access to " << filename << "..." << std::endl;
    std::this_thread::sleep_for( std::chrono::milliseconds(1000) );
    file.open(filename.c_str());
    iter++;
  }

}

}
