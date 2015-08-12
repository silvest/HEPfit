/*
 * LhcoParser.cpp
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
  LhcoParser::LhcoParser( std::string const eventFileName,
                          bool const isVerbose ) :
      isVerbose( isVerbose ),
      fileParser( "#",
                  this->isVerbose ),
      currentEvent( this->isVerbose ),
      fileIsOpenAndNotYetAtEndOfFile( false ),
      currentLine( "" ),
      lineNumber( -1 )
  {
    if( !(eventFileName.empty()) )
    {
      openFile( eventFileName );
    }
  }

  LhcoParser::~LhcoParser()
  {
    // does nothing.
  }


  bool
  LhcoParser::readNextEvent()
  // this reads in the next event in the event file, & returns true if
  // successful.
  {
    if( fileIsOpenAndNotYetAtEndOfFile )
    {
      currentEvent.prepareForNextEvent();
      lineNumber = -1;
      // this is so that the following while loop doesn't terminate until the
      // *next* "line 0", which indicates the start of the next event.
      while( fileIsOpenAndNotYetAtEndOfFile
             &&
             ( 0 != lineNumber ) )
      {
        fileIsOpenAndNotYetAtEndOfFile
        = fileParser.readJustNextValidLine( currentLine );
        if( fileIsOpenAndNotYetAtEndOfFile )
        {
          lineNumber = currentEvent.recordLine( currentLine );
          if( ( 0 > lineNumber )
              &&
              isVerbose )
          {
            std::cout
            << std::endl
            << "LHPC::warning! LHCO parser ignoring invalid line \""
            << currentLine << "\"";
            std::cout << std::endl;
          }
        }
      }
      // at this point, either a new "line 0" has been found, priming the next
      // event, or the last line of the file has been read. if the file has
      // been read to the end, it is closed:
      if( !fileIsOpenAndNotYetAtEndOfFile )
      {
        fileParser.closeFile();
      }
      return true;
    }
    else
    {
      return false;
    }
  }

}
