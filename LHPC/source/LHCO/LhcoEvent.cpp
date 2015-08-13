/*
 * LhcoEvent.cpp
 *
 *  Created on: Jun 29, 2012
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
    std::string const LhcoEvent::trimmingChars( " \t\r\n" );
    bool const LhcoEvent::trueForVerbosity( true );
    std::string const LhcoEvent::eventCommentLine(
     "#   typ     eta    phi       pt  jmass  ntrk  btag   had/em  dummy dummy"
                                                                             );
    int const LhcoEvent::charactersForEventNumber( 14 );
    int const LhcoEvent::charactersForTriggerWord( 6 );

    LhcoEvent::LhcoEvent( bool const isVerbose ) :
        eventNumber( 0 ),
        triggerWord( 0 ),
        nextEventNumber( 0 ),
        nextTriggerWord( 0 ),
        objectLines(),
        objectLineNumber( -1 ),
        missingEnergyLinePointer( NULL ),
        objectLists( (unsigned int)missingEnergyObject + 1 ),
        isVerbose( isVerbose ),
        eventAsString( "" ),
        interpretingObjectLine( NULL ),
        lineParser()
    {
      // just an initialization list.
    }

    LhcoEvent::LhcoEvent( LhcoEvent const& copySource ) :
        eventNumber( copySource.eventNumber ),
        triggerWord( copySource.triggerWord ),
        nextEventNumber( copySource.nextEventNumber ),
        nextTriggerWord( copySource.nextTriggerWord ),
        objectLines( copySource.objectLines,
                     &ObjectLine::copyObjectLine ),
        objectLineNumber( -1 ),
        missingEnergyLinePointer( NULL ),
        objectLists( copySource.objectLists.size() ),
        isVerbose( trueForVerbosity ),
        eventAsString( "" ),
        interpretingObjectLine( NULL ),
        lineParser()
    {
      for( int whichLine( 0 );
           objectLines.getLastIndex() >= whichLine;
           ++whichLine )
      {
        interpretingObjectLine = objectLines.getPointer( whichLine );
        getObjectList( interpretingObjectLine->getObjectType()
                                         ).push_back( interpretingObjectLine );
      }
      if( !(getObjectList( (unsigned int)missingEnergyObject ).empty()) )
      {
        missingEnergyLinePointer
        = getObjectList( (unsigned int)missingEnergyObject ).front();
      }
    }

    LhcoEvent::~LhcoEvent()
    {
      // does nothing.
    }


    int
    LhcoEvent::recordLine( std::string const& lineAsString )
    // this parses the 1st word of lineAsString and returns it, first adding
    // a new ObjectLine if appropriate, noting its pointer in the appropriate
    // std::list. if lineAsString is the start of a new event (as the line
    // begins with 0, which is also what this function will return), it
    // prepares nextEventNumber & nextTriggerWord based on this line. -1 is
    // returned in the case of invalid input.
    {
      eventAsString.append( "\n" );
      eventAsString.append( lineAsString );

      lineParser.clearEntries();
      BOL::StringParser::parseByChar( lineAsString,
                                      lineParser,
                                      BOL::StringParser::whitespaceChars );

      if( 2 >= lineParser.getSize() )
        // if the wrong number of data was given (at least 2 numbers are
        // required)...
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "LHPC::warning! \"" << lineAsString
          << "\" is not a valid LHCO line!";
          std::cout << std::endl;
        }
        objectLineNumber = -1;
      }
      else
      {
        objectLineNumber
        = BOL::StringParser::stringToInt( lineParser[ 0 ] );

        if( 0 == objectLineNumber )
        {
          nextEventNumber
          = BOL::StringParser::stringToInt( lineParser[ 1 ] );
          nextTriggerWord
          = BOL::StringParser::stringToInt( lineParser[ 2 ] );
        }
        else
        {
          interpretingObjectLine
          = objectLines.newEnd().recordLine( objectLineNumber,
                                             lineParser );
          if( (int)missingEnergyObject
              == interpretingObjectLine->getObjectType() )
          {
            missingEnergyLinePointer = interpretingObjectLine;
          }
          getObjectList( interpretingObjectLine->getObjectType()
                                         ).push_back( interpretingObjectLine );
        }
      }
      // end of if-else checking that there were enough numbers to form a
      // valid line.
      return objectLineNumber;
    }

  }

}
