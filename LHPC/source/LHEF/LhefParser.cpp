/*
 * LhefParser.cpp
 *
 *  Created on: Jan 11, 2012
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
  std::string const LhefParser::eventTag( "event" );

  LhefParser::LhefParser( std::string const eventFileName,
                          bool const isVerbose ) :
      isVerbose( isVerbose ),
      fileParser( isVerbose ),
      currentEvent( isVerbose ),
      automaticFilters(),
      fileIsOpen( false ),
      eventAsString( "" ),
      eventIsValid( false )
  {
    if( !(eventFileName.empty()) )
    {
      openFile( eventFileName );
    }
  }

  LhefParser::~LhefParser()
  {
    // does nothing.
  }


  bool
  LhefParser::readNextEvent()
  // this reads in the next event in the event file, & returns true if
  // successful.
  {
    if( fileIsOpen )
    {
      fileIsOpen = fileParser.readNextElement();
      while( !(fileParser.currentElementNameMatches( eventTag )) )
      {
        fileIsOpen = fileParser.readNextElement();
        if( !fileIsOpen )
        {
          fileParser.closeFile();
          eventAsString.assign( "" );
          return false;
        }
      }
      eventAsString.assign( fileParser.getCurrentElementContent() );
      eventIsValid = currentEvent.recordEvent( eventAsString );
      for( int filterIndex( automaticFilters.size() - 1 );
           0 <= filterIndex;
           --filterIndex )
      {
        automaticFilters[ filterIndex ]->updateForNewEvent( currentEvent );
      }
    }
    return fileIsOpen;
  }

}
