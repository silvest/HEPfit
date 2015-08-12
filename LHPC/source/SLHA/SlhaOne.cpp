/*
 * SlhaOne.cpp
 *
 *  Created on: Feb 6, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "SLHA.hpp"

namespace LHPC
{
  SlhaOne::SlhaOne( SlhaParser& fileParser,
                    bool const isVerbose ) :
      MODSEL( "MODSEL",
              BOL::UsefulStuff::notANumber,
              isVerbose ),
      SMINPUTS( "SMINPUTS",
                BOL::UsefulStuff::notANumber,
                isVerbose ),
      MINPAR( "MINPAR",
              BOL::UsefulStuff::notANumber,
              isVerbose ),
      EXTPAR( "EXTPAR",
              BOL::UsefulStuff::notANumber,
              isVerbose ),
      MASS( "MASS",
            BOL::UsefulStuff::notANumber,
            isVerbose,
            9 ),
      NMIX( "NMIX",
            0.0,
            isVerbose ),
      UMIX( "UMIX",
            0.0,
            isVerbose ),
      VMIX( "VMIX",
            0.0,
            isVerbose ),
      STOPMIX( "STOPMIX",
               0.0,
               isVerbose ),
      SBOTMIX( "SBOTMIX",
               0.0,
               isVerbose ),
      STAUMIX( "STAUMIX",
               0.0,
               isVerbose ),
      ALPHA( "ALPHA",
             BOL::UsefulStuff::notANumber,
             isVerbose ),
      HMIX( "HMIX",
            BOL::UsefulStuff::notANumber,
            isVerbose ),
      GAUGE( "GAUGE",
             BOL::UsefulStuff::notANumber,
             isVerbose ),
      MSOFT( "MSOFT",
             BOL::UsefulStuff::notANumber,
             isVerbose ),
      AU( "AU",
          0.0,
          isVerbose ),
      AD( "AD",
          0.0,
          isVerbose ),
      AE( "AE",
          0.0,
          isVerbose ),
      YU( "YU",
          0.0,
          isVerbose ),
      YD( "YD",
          0.0,
          isVerbose ),
      YE( "YE",
          0.0,
          isVerbose ),
      SPINFO( "SPINFO",
              isVerbose ),
      fileParser( fileParser ),
      isVerbose( isVerbose )
  {
    fileParser.registerBlock( MODSEL );
    fileParser.registerBlock( SMINPUTS );
    fileParser.registerBlock( MINPAR );
    fileParser.registerBlock( EXTPAR );
    fileParser.registerBlock( MASS );
    fileParser.registerBlock( NMIX );
    fileParser.registerBlock( UMIX );
    fileParser.registerBlock( VMIX );
    fileParser.registerBlock( STOPMIX );
    fileParser.registerBlock( SBOTMIX );
    fileParser.registerBlock( STAUMIX );
    fileParser.registerBlock( ALPHA );
    fileParser.registerBlock( HMIX );
    fileParser.registerBlock( GAUGE );
    fileParser.registerBlock( MSOFT );
    fileParser.registerBlock( AU );
    fileParser.registerBlock( AD );
    fileParser.registerBlock( AE );
    fileParser.registerBlock( YU );
    fileParser.registerBlock( YD );
    fileParser.registerBlock( YE );
    fileParser.registerBlock( SPINFO );
  }

  SlhaOne::~SlhaOne()
  {
    // does nothing.
  }

}
