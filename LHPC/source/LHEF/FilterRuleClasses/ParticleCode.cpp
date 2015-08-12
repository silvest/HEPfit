/*
 * ParticleCode.cpp
 *
 *  Created on: Jan 26, 2012
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
      ParticleCode::ParticleCode( std::vector< int > const& acceptableValues,
                                  bool const acceptRatherThanReject ) :
          FilterRule( acceptRatherThanReject ),
          acceptableValues( acceptableValues )
      {
        // just an initialization list.
      }

      ParticleCode::ParticleCode( int const acceptableValue,
                                  bool const acceptRatherThanReject ) :
          FilterRule( acceptRatherThanReject ),
          acceptableValues( 1,
                            acceptableValue )
      {
        // just an initialization list.
      }

      ParticleCode::~ParticleCode()
      {
        // does nothing.
      }

    }

  }

}
