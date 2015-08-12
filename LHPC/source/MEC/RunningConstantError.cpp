/*
 * RunningConstantError.cpp
 *
 *  Created on: Apr 1, 2012 (really!)
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "MEC.hpp"

namespace LHPC
{
  RunningConstantError::RunningConstantError() :
      minusUncertainty( BOL::UsefulStuff::notANumber ),
      plusUncertainty( BOL::UsefulStuff::notANumber ),
      schemeType( -1 ),
      evaluationScale( BOL::UsefulStuff::notANumber )
  {
    // just an initialization list.
  }

  RunningConstantError::RunningConstantError(
                                     RunningConstantError const& copySource ) :
      minusUncertainty( copySource.minusUncertainty ),
      plusUncertainty( copySource.plusUncertainty ),
      schemeType( copySource.schemeType ),
      evaluationScale( copySource.evaluationScale )
  {
    // just an initialization list.
  }

  RunningConstantError::~RunningConstantError()
  {
    // does nothing.
  }

}
