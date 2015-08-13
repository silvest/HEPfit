/*
 * UsefulStuff.cpp
 *
 *  Created on: Jan 6, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "UsefulStuff.hpp"

namespace BOL
{
  double const UsefulStuff::notANumber( NAN );
  std::string const UsefulStuff::nanString( "NaN (\"Not a Number\")" );
  double const UsefulStuff::twicePi( 2.0 * M_PI );

  bool UsefulStuff::randomSeedNotYetSet( true );

}
