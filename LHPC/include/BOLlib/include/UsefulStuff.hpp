/*
 * UsefulStuff.hpp
 *
 *  Created on: Jan 6, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef USEFULSTUFF_HPP_
#define USEFULSTUFF_HPP_

#include <cstdlib>
#include <ctime>
#include <stdexcept>
#include <string>
#include <cmath>
#include <fstream>

namespace BOL
{
  class UsefulStuff
  {
  public:
    static double const notANumber;
    static std::string const nanString;
    static double const twicePi;

    static bool
    doublesCloseEnough( double const firstDouble,
                        double const secondDouble,
                        double maximumDifference );
    // this returns true if the difference between firstDouble & secondDouble
    // is less than maximumDifference.
    static double
    flatRandomDouble( double const lowerLimit,
                      double const upperLimit );
    /* this returns a double from a flat probability distribution from the
     * *inclusive* lower limit of lowerLimit to the *exclusive* upper limit
     * of upperLimit.
     */
    static int
    zeroOrOne();
    // this returns 0 or 1 with a 50:50 chance for each.
    static int
    plusOrMinusOne();
    // this returns -1 or +1 with a 50:50 chance for each.
    static bool
    fileExists( std::string const& fileName );
    // returns true if a file with the name fileName exists, false otherwise.
    static void
    runSystemCommand( std::string const& systemCommand );
    // this runs system( systemCommand.c_str() ) & throws an exception if it
    // returned -1.

  private:
    static bool randomSeedNotYetSet;

    static void
    ensureRandomSeedIsSet();
  };



  inline bool
  UsefulStuff::doublesCloseEnough( double const firstDouble,
                                   double const secondDouble,
                                   double maximumDifference )
  // this returns true if the difference between firstDouble & secondDouble
  // is less than maximumDifference.
  {
    return ( abs( maximumDifference ) >= abs( firstDouble - secondDouble ) );
  }

  inline double
  UsefulStuff::flatRandomDouble( double const lowerLimit,
                                 double const upperLimit )
  /* this returns a double from a flat probability distribution from the
   * *inclusive* lower limit of lower_limit to the *exclusive* upper limit
   * of upper_limit.
   */
  {
    ensureRandomSeedIsSet();
    return ( lowerLimit + ( (double)(rand()) / (double)RAND_MAX )
                          * ( upperLimit - lowerLimit ) );
  }

  inline int
  UsefulStuff::zeroOrOne()
  // this returns 0 or 1 with a 50:50 chance for each.
  {
    ensureRandomSeedIsSet();
    return ( (rand()) % 2 );
  }

  inline int
  UsefulStuff::plusOrMinusOne()
  // this returns -1 or +1 with a 50:50 chance for each.
  {
    if( 0 == zeroOrOne() )
    {
      return -1;
    }
    else
    {
      return 1;
    }
  }

  inline void
  UsefulStuff::ensureRandomSeedIsSet()
  // this sets the random seed if it had not already been set.
  {
    if( randomSeedNotYetSet )
    {
      // debugging:
      /**time_t currentTime = time( NULL );
      srand( currentTime );
      std::cout << std::endl << "currentTime = " << currentTime;**/
      srand( time( NULL ) );
      randomSeedNotYetSet = false;
    }
  }

  inline bool
  UsefulStuff::fileExists( std::string const& fileName )
  // returns true if a file with the name fileName exists, false otherwise.
  {
    std::ifstream fileStream( fileName.c_str() );
    if( fileStream.good() )
    {
      fileStream.close();
      return true;
    }
    else
    {
      return false;
    }
  }

  inline void
  UsefulStuff::runSystemCommand( std::string const& systemCommand )
  {
    int systemReturn( system( systemCommand.c_str() ) );
    if( -1 == systemReturn )
    {
      throw std::runtime_error( "system( \"" +  systemCommand
                                + "\" ) returned -1" );
    }
  }

}

#endif /* USEFULSTUFF_HPP_ */
