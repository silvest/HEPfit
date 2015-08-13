/*
 * ParticleLine.cpp
 *
 *  Created on: Jan 25, 2012
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
    ParticleLine::ParticleLine() :
        particleCode( (int)BOL::UsefulStuff::notANumber ),
        initialOrIntermediateOrFinalState( (int)BOL::UsefulStuff::notANumber ),
        primaryMotherLineNumber( (int)BOL::UsefulStuff::notANumber ),
        secondaryMotherLineNumber( (int)BOL::UsefulStuff::notANumber ),
        colorCode( (int)BOL::UsefulStuff::notANumber ),
        anticolorCode( (int)BOL::UsefulStuff::notANumber ),
        xMomentum( BOL::UsefulStuff::notANumber ),
        yMomentum( BOL::UsefulStuff::notANumber ),
        zMomentum( BOL::UsefulStuff::notANumber ),
        tMomentum( BOL::UsefulStuff::notANumber ),
        restMass( BOL::UsefulStuff::notANumber ),
        properLifetime( BOL::UsefulStuff::notANumber ),
        twiceSpin( (int)BOL::UsefulStuff::notANumber ),
        ownLineNumber( (int)BOL::UsefulStuff::notANumber ),
        primaryMotherLinePointer( NULL ),
        secondaryMotherLinePointer( NULL ),
        daughterLines(),
        lineParser()
    {
      // just an initialization list.
    }

    ParticleLine::ParticleLine( ParticleLine const& copySource ) :
        particleCode( copySource.particleCode ),
        initialOrIntermediateOrFinalState(
                                copySource.initialOrIntermediateOrFinalState ),
        primaryMotherLineNumber( copySource.primaryMotherLineNumber ),
        secondaryMotherLineNumber( copySource.secondaryMotherLineNumber ),
        colorCode( copySource.colorCode ),
        anticolorCode( copySource.anticolorCode ),
        xMomentum( copySource.xMomentum ),
        yMomentum( copySource.yMomentum ),
        zMomentum( copySource.zMomentum ),
        tMomentum( copySource.tMomentum),
        restMass( copySource.restMass ),
        properLifetime( copySource.properLifetime ),
        twiceSpin( copySource.twiceSpin ),
        ownLineNumber( (int)BOL::UsefulStuff::notANumber ),
        primaryMotherLinePointer( NULL ),
        secondaryMotherLinePointer( NULL ),
        daughterLines(),
        lineParser()
    {
      /* just an initialization list. unfortunately there's no easy way to copy
       * the mother & daughter pointers. however, that should be taken care of
       * by the copy constructor of the LhefEvent that is copying its
       * ParticleLines.
       */
    }

    ParticleLine::~ParticleLine()
    {
      // does nothing.
    }


    double
    ParticleLine::IPUP( int const whichComponent ) const
    {
      if( 1 == whichComponent )
      {
        return getXMomentum();
      }
      else if( 2 == whichComponent )
      {
        return getYMomentum();
      }
      else if( 3 == whichComponent )
      {
        return getZMomentum();
      }
      else if( 4 == whichComponent )
      {
        return getEnergy();
      }
      else if( 5 == whichComponent )
      {
        return getMass();
      }
      else
      {
        return BOL::UsefulStuff::notANumber;
      }
    }

    bool
    ParticleLine::recordLine( int const ownLineNumber,
                              std::string const& lineAsString )
    /* this interprets a string as the data it is meant to correspond to. it
     * returns false if the wrong number of data were given, & in that case
     * it sets all its entries to BOL::UsefulStuff::reallyWrongValue. it
     * returns true if all the entries were filled as expected. this function
     * also clears daughterLines.
     */
    {
      this->ownLineNumber = ownLineNumber;
      daughterLines.clear();
      primaryMotherLinePointer = NULL;
      secondaryMotherLinePointer = NULL;
      lineParser.clearEntries();
      BOL::StringParser::parseByChar( lineAsString,
                                      lineParser,
                                      BOL::StringParser::whitespaceChars );
      if( 13 != lineParser.getSize() )
        // if the wrong number of data was given...
      {
        particleCode = (int)BOL::UsefulStuff::notANumber;
        initialOrIntermediateOrFinalState = (int)BOL::UsefulStuff::notANumber;
        primaryMotherLineNumber = (int)BOL::UsefulStuff::notANumber;
        secondaryMotherLineNumber = (int)BOL::UsefulStuff::notANumber;
        colorCode = (int)BOL::UsefulStuff::notANumber;
        anticolorCode = (int)BOL::UsefulStuff::notANumber;
        xMomentum = BOL::UsefulStuff::notANumber;
        yMomentum = BOL::UsefulStuff::notANumber;
        zMomentum = BOL::UsefulStuff::notANumber;
        tMomentum = BOL::UsefulStuff::notANumber;
        restMass = BOL::UsefulStuff::notANumber;
        properLifetime = BOL::UsefulStuff::notANumber;
        twiceSpin = (int)BOL::UsefulStuff::notANumber;
        return false;
      }
      else
      {
        particleCode = BOL::StringParser::stringToInt( lineParser[ 0 ] );
        initialOrIntermediateOrFinalState
        = BOL::StringParser::stringToInt( lineParser[ 1 ] );
        primaryMotherLineNumber
        = BOL::StringParser::stringToInt( lineParser[ 2 ] );
        secondaryMotherLineNumber
        = BOL::StringParser::stringToInt( lineParser[ 3 ] );
        colorCode = BOL::StringParser::stringToInt( lineParser[ 4 ] );
        anticolorCode = BOL::StringParser::stringToInt( lineParser[ 5 ] );
        xMomentum = BOL::StringParser::stringToDouble( lineParser[ 6 ] );
        yMomentum = BOL::StringParser::stringToDouble( lineParser[ 7 ] );
        zMomentum = BOL::StringParser::stringToDouble( lineParser[ 8 ] );
        tMomentum = BOL::StringParser::stringToDouble( lineParser[ 9 ] );
        restMass = BOL::StringParser::stringToDouble( lineParser[ 10 ] );
        properLifetime = BOL::StringParser::stringToDouble( lineParser[ 11 ] );
        twiceSpin = BOL::StringParser::stringToInt( lineParser[ 12 ] );
        return true;
      }
    }

  }

}
