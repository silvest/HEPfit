/*
 * BasicTimer.hpp
 *
 *  Created on: Oct 23, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef BASICTIMER_HPP_
#define BASICTIMER_HPP_

#include <cstdlib>
#include <sys/time.h>

namespace BOL
{
  // this is a class to do some basic timing stuff.
  class BasicTimer
  {
  public:
    static double
    secondsSince( timeval const& referenceTimeval );
    // this returns the number of seconds since referenceTimeval as a double.

    BasicTimer( double givenDuration = 0.0,
                timeval const* const startTimeval = NULL );
    ~BasicTimer();

    void
    setStartTime( timeval const* const startTimeval = NULL );
    // if startTimeval is provided, that becomes the new reference "time = 0".
    // if a NULL pointer is provided, then the current time is taken and saved.
    double
    secondsSinceStart() const;
    void
    setDuration( double const givenDuration );
    bool
    withinDuration() const;
    // this returns true if the current time is not later than givenDuration
    // seconds after startTimeval.


  protected:
    double givenDuration;
    // this is a given duration in seconds for whether it is that much later
    // than the reference start time.
    timeval startTimeval;
    // this is the reference start time which time durations are calculated
    // from.
  };



  inline double
  BasicTimer::secondsSince( timeval const& referenceTimeval )
  // this returns the number of seconds since referenceTimeval as a double.
  {
    timeval currentTimeval;
    gettimeofday( &currentTimeval,
                  NULL );
    return ( (double)( currentTimeval.tv_sec - referenceTimeval.tv_sec )
             + 0.000001 * (double)( currentTimeval.tv_usec
                                    - referenceTimeval.tv_usec ) );
  }

  inline void
  BasicTimer::setStartTime( timeval const* const startTimeval )
  // if startTimeval is provided, that becomes the new reference "time = 0".
  // if a NULL pointer is provided, then the current time is taken and saved.
  {
    if( NULL == startTimeval )
    {
      gettimeofday( &(this->startTimeval),
                    NULL );
    }
    else
    {
      this->startTimeval = *startTimeval;
    }
  }

  inline double
  BasicTimer::secondsSinceStart() const
  {
    return secondsSince( startTimeval );
  }

  inline void
  BasicTimer::setDuration( double const givenDuration )
  {
    this->givenDuration = givenDuration;
  }

  inline bool
  BasicTimer::withinDuration() const
  // this returns true if the current time is not later than givenDuration
  // seconds after startTimeval.
  {
    return ( givenDuration >= secondsSinceStart() );
  }

} /* namespace BOL */
#endif /* BASICTIMER_HPP_ */
