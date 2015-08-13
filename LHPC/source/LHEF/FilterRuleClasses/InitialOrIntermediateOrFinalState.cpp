/*
 * InitialOrIntermediateOrFinalState.cpp
 *
 *  Created on: Jan 30, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "LHEF.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace FilterRuleClass
    {
      InitialOrIntermediateOrFinalState::InitialOrIntermediateOrFinalState(
                                               stateType const acceptableValue,
                                          bool const acceptRatherThanReject ) :
          FilterRule( acceptRatherThanReject ),
          acceptableValue( (int)acceptableValue )
      {
        // just an initialization list.
      }

      InitialOrIntermediateOrFinalState::InitialOrIntermediateOrFinalState(
                                                     int const acceptableValue,
                                          bool const acceptRatherThanReject ) :
          FilterRule( acceptRatherThanReject ),
          acceptableValue( acceptableValue )
      {
        // just an initialization list.
      }

      InitialOrIntermediateOrFinalState::~InitialOrIntermediateOrFinalState()
      {
        // does nothing.
      }

    }

  }

}
