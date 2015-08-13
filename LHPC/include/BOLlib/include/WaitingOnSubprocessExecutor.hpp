/*
 * WaitingOnSubprocessExecutor.hpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef WAITINGONSUBPROCESSEXECUTOR_HPP_
#define WAITINGONSUBPROCESSEXECUTOR_HPP_

#include <csignal>
#include <ctime>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include "StringParser.hpp"
#include "VectorlikeArray.hpp"

namespace BOL
{
  class WaitingOnSubprocessExecutor
  {
  public:
    WaitingOnSubprocessExecutor(
                         std::string const& executableStringIncludingArguments,
                                 bool const isVerbose,
                                 int const patienceMilliseconds = 10000 );
    ~WaitingOnSubprocessExecutor();

    void
    setExecutableName( std::string const& executableName );
    void
    setPatienceTicks( int const patienceTicks );
    /* the setArguments functions assign memory to make argumentListAsArray
     * point to an array of pointers to chars which are stored in the strings
     * of executableNameAndArguments & thus can be passed to execv to properly
     * execute the subprocess with its arguments.
     */
    void
    setArguments( std::string const& argumentsAsSingleString );
    bool
    forkAndExecvAndWait();


  protected:
    static int const sleepingMicroseconds;

    int processId;
    // this is for checking whether the process is the parent or the forked
    // child subprocess.
    bool const isVerbose;
    // this is to allow the suppression of the generated reports.
    VectorlikeArray< std::string > executableNameAndArguments;
    char const** argumentListAsArray;
    /* this is to hold the arguments for the executable. it is a constant
     * pointer to a pointer to a char because execv is ancient & requires
     * arguments in an irritatingly silly format.
     */
    clock_t patienceTicks;
    /* this is the number of clock ticks for which the parent process will
     * wait before losing patience & killing the subprocess, if the subprocess
     * hasn't finished.
     */
    clock_t killingTime;
    // this is the clock() value at which the parent process has run out of
    // patience.

    void
    setUpArgumentCharArray();
  };



  inline void
  WaitingOnSubprocessExecutor::setExecutableName(
                                            std::string const& executableName )
  {
    executableNameAndArguments.getFront().assign( executableName );
    setUpArgumentCharArray();
  }

  inline void
  WaitingOnSubprocessExecutor::setPatienceTicks( int const patienceTicks )
  {
    this->patienceTicks = patienceTicks;
  }


  inline void
  WaitingOnSubprocessExecutor::setArguments(
                                   std::string const& argumentsAsSingleString )
  /* the setArguments functions assign memory to make argumentListAsArray
   * point to an array of pointers to chars which are stored in the strings
   * of executableNameAndArguments & thus can be passed to execv to properly
   * execute the subprocess with its arguments.
   */
  {
    executableNameAndArguments.setSize( 1 );
    StringParser::parseByChar( argumentsAsSingleString,
                               executableNameAndArguments );
    setUpArgumentCharArray();
  }

  inline void
  WaitingOnSubprocessExecutor::setUpArgumentCharArray()
  {
    delete[] argumentListAsArray;
    argumentListAsArray
    = new char const*[ ( executableNameAndArguments.getSize() + 1 ) ];
    // create a char* const array which has room for the executable name,
    // then all the arguments, then NULL.
    for( int argumentCharArrayCounter( 0 );
         executableNameAndArguments.getSize() > argumentCharArrayCounter;
         ++argumentCharArrayCounter )
      // go through the list of arguments.
    {
      argumentListAsArray[ argumentCharArrayCounter ]
      = executableNameAndArguments[ argumentCharArrayCounter ].c_str();
      // write this argument as the (argumentCharArrayCounter)th element.
    }
    argumentListAsArray[ executableNameAndArguments.getSize() ] = NULL;
    // we fill the (argumentCharArrayCounter + 1)th element with NULL, since
    // that is what execv() wants.
  }

}

#endif /* WAITINGONSUBPROCESSEXECUTOR_HPP_ */
