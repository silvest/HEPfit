/*
 * BasicTimer.cpp
 *
 *  Created on: Oct 23, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "../include/BasicTimer.hpp"

namespace BOL
{

  BasicTimer::BasicTimer( double givenDuration,
                          timeval const* const startTimeval ):
    givenDuration( givenDuration ),
    startTimeval()
  {
    setStartTime( startTimeval );
  }

  BasicTimer::~BasicTimer()
  {
    // does nothing.
  }

} /* namespace BOL */
