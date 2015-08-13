/*
 * LhefEvent.cpp
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
    std::string const LhefEvent::trimmingChars( " \t\r\n" );
    bool const LhefEvent::trueForVerbosity( true );

    LhefEvent::LhefEvent( bool const isVerbose ) :
        numberOfParticles( (int)BOL::UsefulStuff::notANumber ),
        eventId( (int)BOL::UsefulStuff::notANumber ),
        eventWeight( BOL::UsefulStuff::notANumber ),
        eventScale( BOL::UsefulStuff::notANumber ),
        alphaQed( BOL::UsefulStuff::notANumber ),
        alphaQcd( BOL::UsefulStuff::notANumber ),
        particleLines(),
        eventNumberInFile( 0 ),
        eventAsLines( 5 ),
        lineAsNumbersAsStrings( 6 ),
        recordingSucceeded( false ),
        headerLine( "" ),
        motherLineNumber( 0 ),
        primaryMotherLinePointer( NULL ),
        secondaryMotherLinePointer( NULL ),
        isVerbose( isVerbose ),
        eventAsString( "" ),
        optionalInformation( "" )
    {
      // just an initialization list.
    }

    LhefEvent::LhefEvent( LhefEvent const& copySource ) :
        numberOfParticles( copySource.numberOfParticles ),
        eventId( copySource.eventId ),
        eventWeight( copySource.eventWeight ),
        eventScale( copySource.eventScale ),
        alphaQed( copySource.alphaQed ),
        alphaQcd( copySource.alphaQcd ),
        particleLines( copySource.particleLines,
                       &ParticleLine::copyParticleLine ),
        eventNumberInFile( copySource.eventNumberInFile ),
        eventAsLines( 5 ),
        lineAsNumbersAsStrings( 6 ),
        recordingSucceeded( copySource.recordingSucceeded ),
        headerLine( "" ),
        motherLineNumber( 0 ),
        primaryMotherLinePointer( NULL ),
        secondaryMotherLinePointer( NULL ),
        isVerbose( trueForVerbosity ),
        eventAsString( copySource.eventAsString ),
        optionalInformation( copySource.optionalInformation )
    {
      /* the copy constructors of particleLines cannot get their mother &
       * daughter pointers by themselves, so this constructor has to sort them
       * out:
       */
      setUpPointersBetweenParticleLines();
    }

    LhefEvent::~LhefEvent()
    {
      // does nothing.
    }


    bool
    LhefEvent::recordEvent( std::string const& eventAsString )
    /* this interprets a string as the data it is meant to correspond to, &
     * ensures that the ParticleLine instances correctly point to each other.
     * if any of the lines had the wrong number of data entries, or if the
     * header number of particles didn't correspond to the number of
     * particle lines, setAsInvalid() is called & then false is returned.
     * true is returned if all the entries were filled as expected.
     */
    {
      this->eventAsString.assign( BOL::StringParser::trimFromFrontAndBack(
                                                                 eventAsString,
                                                             trimmingChars ) );
      ++eventNumberInFile;
      eventAsLines.clearEntries();
      size_t startOfOptionalInformation( this->eventAsString.find( '#' ) );
      if( std::string::npos != startOfOptionalInformation )
      {
        optionalInformation.assign(
                      this->eventAsString.substr( startOfOptionalInformation ) );
        BOL::StringParser::parseByChar(
                                       BOL::StringParser::trimFromFrontAndBack(
                                                 this->eventAsString.substr( 0,
                                                  startOfOptionalInformation ),
                                                               trimmingChars ),
                                        eventAsLines,
                                        BOL::StringParser::newlineChars );
      }
      else
      {
        optionalInformation.assign( "" );
        BOL::StringParser::parseByChar( this->eventAsString,
                                        eventAsLines,
                                        BOL::StringParser::newlineChars );
      }
      // now lineParser should be 1 string for the header line &
      // ( lineParser.getSize() - 1 ) strings for the particle lines.
      if( 1 >= eventAsLines.getSize() )
        // if the wrong number of data was given (header line & at least 1
        // particle line are required)...
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "LHPC::warning! recording event " << eventNumberInFile
          << " as invalid because only " << eventAsLines.getSize()
          << " lines were found in the event, which is not enough for a header"
          << " line + at least 1 particle line.";
          std::cout << std::endl;
        }
        setAsInvalid();
        return false;
      }
      else
      {
        // 1st the header line is set up:
        headerLine.assign( eventAsLines.getFront() );
        lineAsNumbersAsStrings.clearEntries();
        BOL::StringParser::parseByChar( headerLine,
                                        lineAsNumbersAsStrings,
                                        BOL::StringParser::whitespaceChars );

        if( 6 != lineAsNumbersAsStrings.getSize() )
          // if the wrong number of data was given (a header line & at least 1
          // particle line are required)...
        {
          if( isVerbose )
          {
            std::cout
            << std::endl
            << "LHPC::warning! recording event " << eventNumberInFile
            << " as invalid because " << lineAsNumbersAsStrings.getSize()
            << " numbers were found for the header line, instead of exactly 6"
            << " numbers, as required.";
            std::cout << std::endl;
          }
          setAsInvalid();
          return false;
        }
        else
        {
          numberOfParticles
          = BOL::StringParser::stringToInt( lineAsNumbersAsStrings[ 0 ] );
          if( numberOfParticles != ( eventAsLines.getSize() - 1 ) )
            // if the header is wrong about how many particle lines there
            // are...
          {
            if( isVerbose )
            {
              std::cout
              << std::endl
              << "LHPC::warning! recording event " << eventNumberInFile
              << " as invalid because the header line declared "
              << numberOfParticles << " particles in the event, but "
              << ( eventAsLines.getSize() - 1 ) << " lines were found for the"
              << " particles.";
              std::cout << std::endl;
            }
            setAsInvalid();
            return false;
          }
          else
          {
            eventId
            = BOL::StringParser::stringToInt( lineAsNumbersAsStrings[ 1 ] );
            eventWeight
            = BOL::StringParser::stringToDouble( lineAsNumbersAsStrings[ 2 ] );
            eventScale
            = BOL::StringParser::stringToDouble( lineAsNumbersAsStrings[ 3 ] );
            alphaQed
            = BOL::StringParser::stringToDouble( lineAsNumbersAsStrings[ 4 ] );
            alphaQcd
            = BOL::StringParser::stringToDouble( lineAsNumbersAsStrings[ 5 ] );

            // then the ParticleLines are set up:
            particleLines.setSize( numberOfParticles );
            recordingSucceeded = true;
            for( int lineCount( 1 );
                 ( recordingSucceeded
                 &&
                 ( lineCount <= numberOfParticles ) );
                 ++lineCount )
            {
              recordingSucceeded
              = particleLines[ lineCount - 1 ].recordLine( lineCount,
                                               eventAsLines[ lineCount ] );
              if( !recordingSucceeded
                  &&
                  isVerbose )
              {
                BOL::StringParser::parseByChar( eventAsLines[ lineCount ],
                                                lineAsNumbersAsStrings,
                                          BOL::StringParser::whitespaceChars );
                std::cout
                << std::endl
                << "LHPC::warning! recording event " << eventNumberInFile
                << " as invalid because particle line " << lineCount
                << " consisted of " << lineAsNumbersAsStrings.getSize()
                << " numbers rather than the required 13 numbers.";
                std::cout << std::endl;
              }
            }
            if( !recordingSucceeded )
            {
              // the warning about an invalid event is given in the loop rather
              // than here so that it can state which line was problematic.
              setAsInvalid();
              return false;
            }
            else
            {
              setUpPointersBetweenParticleLines();
              return true;
            }
            // end of if-else to finally set up the ParticleLines after
            // checking that they all recorded properly.

          }
          // end of if-else checking that the number of particles in the header
          // matches the number of lines for particles in the event.
        }  // end of if-else checking that the header is correctly formed.
      }
      // end of if-else checking that there were enough lines to form an event.
    }

    void
    LhefEvent::setUpPointersBetweenParticleLines()
    // this sets up the mother & daughter pointers for each ParticleLine in
    // particleLines.
    {
      for( int lineCount( 0 );
           lineCount < numberOfParticles;
           ++lineCount )
      {
        motherLineNumber
        = particleLines[ lineCount ].getPrimaryMotherLineNumber();
        if( motherLineNumber > 0 )
        {
          primaryMotherLinePointer
          = particleLines.getPointer( motherLineNumber - 1 );
        }
        else
        {
          primaryMotherLinePointer = NULL;
        }
        motherLineNumber
        = particleLines[ lineCount ].getSecondaryMotherLineNumber();
        if( motherLineNumber > 0 )
        {
          secondaryMotherLinePointer
          = particleLines.getPointer( motherLineNumber - 1 );
        }
        else
        {
          secondaryMotherLinePointer = NULL;
        }
        particleLines[ lineCount ].setMotherLines( primaryMotherLinePointer,
                                                  secondaryMotherLinePointer );
      }
    }

  }

}
