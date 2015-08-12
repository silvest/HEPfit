/*
 * ObjectLine.cpp
 *
 *  Created on: Jun 25, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "LHCO.hpp"

namespace LHPC
{
  namespace LHCO
  {
   int const ObjectLine::minimumNumberOfEntries( 9 );
   // this is the minimum number of numbers expected for a valid line.

    ObjectLine::ObjectLine() :
        lineNumber( -1 ),
        objectType( -1 ),
        valueVector( 11,
                     BOL::UsefulStuff::notANumber )
    {
      // just an initialization list.
    }

    ObjectLine::ObjectLine( ObjectLine const& copySource ) :
        lineNumber( copySource.lineNumber ),
        objectType( copySource.objectType ),
        valueVector( copySource.valueVector )
    {
      // just an initialization list.
    }

    ObjectLine::~ObjectLine()
    {
      // does nothing.
    }


    ObjectLine const*
    ObjectLine::recordLine( int const lineNumber,
                     BOL::VectorlikeArray< std::string > const& lineAsStrings )
    // this interprets a string as the data it is meant to correspond to. it
    // assumes that it has been given enough numbers.
    {
      if( minimumNumberOfEntries <= lineAsStrings.getSize() )
      {
        valueVector.assign( lineAsStrings.getSize(),
                            BOL::UsefulStuff::notANumber );
      }
      else
      {
        valueVector.assign( minimumNumberOfEntries,
                            BOL::UsefulStuff::notANumber );
      }
      this->lineNumber = lineNumber;
      valueVector[ 0 ] = (double)lineNumber;
      objectType = BOL::StringParser::stringToInt( lineAsStrings[ 1 ] );
      valueVector[ 1 ] = (double)objectType;
      for( int whichElement( 2 );
           lineAsStrings.getLastIndex() >= whichElement;
           ++whichElement )
      {
        valueVector[ whichElement ]
        = BOL::StringParser::stringToDouble( lineAsStrings[ whichElement ] );
      }
      return this;
    }

  }

}
