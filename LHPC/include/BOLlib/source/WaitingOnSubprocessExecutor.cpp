/*
 * WaitingOnSubprocessExecutor.cpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "WaitingOnSubprocessExecutor.hpp"

namespace BOL
{
  int const WaitingOnSubprocessExecutor::sleepingMicroseconds( 1000 );

  WaitingOnSubprocessExecutor::WaitingOnSubprocessExecutor(
                         std::string const& executableStringIncludingArguments,
                                                          bool const isVerbose,
                                             int const patienceMilliseconds ) :
      processId( 0 ),
      isVerbose( isVerbose ),
      executableNameAndArguments( 0 ),
      argumentListAsArray( NULL ),
      patienceTicks( ( ( patienceMilliseconds * CLOCKS_PER_SEC ) / 1000 ) ),
      killingTime()
  {
    StringParser::parseByChar( executableStringIncludingArguments,
                               executableNameAndArguments );
    setUpArgumentCharArray();
  }

  WaitingOnSubprocessExecutor::~WaitingOnSubprocessExecutor()
  {
    delete[] argumentListAsArray;
  }


  bool
  WaitingOnSubprocessExecutor::forkAndExecvAndWait()
  /* this forks off a new process, then checks if this process is the parent
   * or the child subprocess. if it's the parent, it waits until the child
   * returns a value or until the "patience" time runs out, at which point it
   * kills the child subprocess if it has not returned a value. if it is the
   * child, it runs execv with its arguments & returns the return value from
   * the execv process.
   */
  {
    killingTime = ( clock() + patienceTicks );
    int processStatus( 0 );
    // this keeps track of the status of the process.
    int returnValue( 0 );
    // this keeps track of what processes return.
    bool subprocessReturned( true );
    // assume that the subprocess won't have to be killed.
    int forkReturnValue( fork() );
    if( 0 > forkReturnValue )
      // if fork failed to do what it should do...
    {
      std::cout
      << std::endl
      << "BOL::error! WaitingOnSubprocessExecutor::forkAndExecvAndWait() was"
      << " unable to create the subprocess; there was an error and fork()"
      << " returned " << forkReturnValue << std::endl;
      // print error message.
    }
    else if( 0 == forkReturnValue )
      // otherwise if this is the subprocess
      // [since it got the 0 from fork()]...
    {
      processId = getpid();
      if( isVerbose )
      {
        std::cout
        << std::endl
        << "now executing \"" << argumentListAsArray[ 0 ];
        for( int argumentWritingCounter( 1 );
             executableNameAndArguments.getSize() > argumentWritingCounter;
             ++argumentWritingCounter )
        {
          std::cout
          << " " << argumentListAsArray[ argumentWritingCounter ];
        }
        std::cout
        << "\", with process id " << getpid() << ", real user id " << getuid()
        << " & effective user id " << geteuid() << std::endl;
      }
      returnValue = execv( argumentListAsArray[ 0 ],
                           (char* const*)argumentListAsArray );
      // run the executable with its arguments.
      exit( returnValue );
    }
    else if( 0 < forkReturnValue )
      // otherwise this is the parent process, since it got given the process
      // id of the subprocess.
    {
      bool subprocessFinished( false );
      // the tracker of whether the subprocess has changed state or not.
      while( clock() < killingTime )
      {
        if( waitpid( forkReturnValue,
                     &processStatus,
                     WNOHANG ) == forkReturnValue )
          /* if the waitpid() function waited for the subprocess (with
           * process id fork_return) and the subprocess changed state, so
           * that waitpid() returned the subprocess's process id...
           */
        {
          killingTime = clock() - 1;
          // note that we don't have to wait any more.
          subprocessFinished = true;
          // note that the subprocess changed state.
        }
        else
        {
          usleep( sleepingMicroseconds );
          // wait 1000 microseconds = 1 millisecond.
        }
      }
      if( false == subprocessFinished )
        // if we have waited over (patienceTicks/1000) seconds & the
        // subprocess has not changed state...
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "killing forked process " << forkReturnValue
          << " due to it taking over " << patienceTicks << " milliseconds"
          << std::endl;
          // print warning about hanging subprocess.
        }
        kill( forkReturnValue,
              9 );
        // kill the hanging process.
        waitpid( -1,
                 &processStatus,
                 0 );
        // wait for all this killing to finish?  I'm not entirely sure.
        subprocessReturned = false;
        // note that the subprocess had to be killed because it did not
        // return a value.
      }
      returnValue = WEXITSTATUS( processStatus );
    }  // end of fork() statement.
    return subprocessReturned;
  }

}
