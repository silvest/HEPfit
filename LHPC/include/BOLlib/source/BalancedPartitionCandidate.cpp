/*
 * BalancedPartitionCandidate.cpp
 *
 *  Created on: Jul 13, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "BalancedPartitionCandidate.hpp"

namespace BOL
{
  BalancedPartitionCandidate::BalancedPartitionCandidate() :
      indexVector(),
      candidateWeight( 0.0 ),
      inverseWeight( UsefulStuff::notANumber )
  {
    // just lets everything initialize in the default way.
  }

  BalancedPartitionCandidate::~BalancedPartitionCandidate()
  {
    // does nothing.
  }

}
