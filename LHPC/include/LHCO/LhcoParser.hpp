/*
 * LhcoParser.hpp
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

#ifndef LHCOPARSER_HPP_
#define LHCOPARSER_HPP_

#include "BOLlib/include/BOLlib.hpp"
#include "LhcoEvent.hpp"

namespace LHPC
{
  // this is a class for reading in a LHCO format file & parsing the events
  // from it.
  class LhcoParser
  {
  public:
    LhcoParser( std::string const eventFileName = "",
                bool const isVerbose = true );
    ~LhcoParser();

    bool
    openFile( std::string const& eventFileName );
    // this opens the file with the given name as the source of its events &
    // returns true, unless the file could not be opened.
    bool
    readNextEvent();
    // this reads in the next event in the event file, & returns true if
    // successful.
    LHCO::LhcoEvent const&
    getEvent() const;
    // this returns the last parsed event.


  protected:
    bool isVerbose;
    BOL::CommentedTextParser fileParser;
    LHCO::LhcoEvent currentEvent;
    bool fileIsOpenAndNotYetAtEndOfFile;
    std::string currentLine;
    int lineNumber;

    bool
    lookForFirstEvent();
    // this reads in to the 1st non-commented line that begins with '0' with
    // optional whitespace before it & at least 1 whitespace-type character
    // following it. it returns true if it found such a line, which it also
    // assigns to currentLine.
  };





  inline bool
  LhcoParser::openFile( std::string const& eventFileName )
  // this opens the file with the given name as the source of its events &
  // returns true, unless the file could not be opened. it also reads in up to
  // the 1st event (preparing currentEvent with the 1st "line 0" of an event).
  {
    fileIsOpenAndNotYetAtEndOfFile = fileParser.openFile( eventFileName );
    if( fileIsOpenAndNotYetAtEndOfFile )
    {
      fileIsOpenAndNotYetAtEndOfFile = lookForFirstEvent();
      if( fileIsOpenAndNotYetAtEndOfFile )
      {
        lineNumber = currentEvent.recordLine( currentLine );
      }
    }
    return fileIsOpenAndNotYetAtEndOfFile;
  }

  inline LHCO::LhcoEvent const&
  LhcoParser::getEvent() const
  // this returns the last parsed event.
  {
    return currentEvent;
  }

  inline bool
  LhcoParser::lookForFirstEvent()
  // this reads in to the 1st non-commented line that begins with '0' with
  // optional whitespace before it & at least 1 whitespace-type character
  // following it. it returns true if it found such a line, which it also
  // assigns to currentLine.
  {
    if( fileIsOpenAndNotYetAtEndOfFile )
    {
      bool notYetFoundFirstLineOfAnEvent( true );
      while( fileIsOpenAndNotYetAtEndOfFile
             &&
             notYetFoundFirstLineOfAnEvent )
      {
        fileIsOpenAndNotYetAtEndOfFile
        = fileParser.readJustNextValidLine( currentLine );
        if( fileIsOpenAndNotYetAtEndOfFile )
        {
          notYetFoundFirstLineOfAnEvent
          = ( 0.0 != BOL::StringParser::stringToDouble( currentLine ) );
          // stringToDouble(...) takes the 1st number found in the string.
        }
      }
      // this loop ends as soon as a line starting with 0 was found, or if the
      // end of the file was reached.
    }
    return fileIsOpenAndNotYetAtEndOfFile;
  }

}

#endif /* LHCOPARSER_HPP_ */
