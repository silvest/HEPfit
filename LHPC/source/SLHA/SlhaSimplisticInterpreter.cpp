/*
 * SlhaSimplisticInterpreter.cpp
 *
 *  Created on: Sep 13, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "SLHA.hpp"

namespace LHPC
{
  SlhaSimplisticInterpreter::SlhaSimplisticInterpreter(
                                            std::string const& slhaFilename ) :
      slhaParser( true,
                  false ),
      stringParser()
  {
    slhaParser.readFile( slhaFilename );
  }

  SlhaSimplisticInterpreter::~SlhaSimplisticInterpreter()
  {
    // does nothing.
  }


  std::string
  SlhaSimplisticInterpreter::operator()( std::string blockNameAndIndices )
  {
    std::string returnString( "" );
    BOL::StringParser::substituteCharacterWith( blockNameAndIndices,
                                                '[',
                                                ' ' );
    BOL::StringParser::substituteCharacterWith( blockNameAndIndices,
                                                ']',
                                                ' ' );
    BOL::StringParser::substituteCharacterWith( blockNameAndIndices,
                                                '(',
                                                ' ' );
    BOL::StringParser::substituteCharacterWith( blockNameAndIndices,
                                                ')',
                                                ' ' );
    std::string indicesString( "" );
    std::string blockName( BOL::StringParser::substringToFirst(
                         BOL::StringParser::trimFromFront( blockNameAndIndices,
                                BOL::StringParser::whitespaceAndNewlineChars ),
                                                                " ",
                                                            &indicesString ) );
    std::vector< int >
    indicesVector( BOL::StringParser::stringToIntVector( indicesString ) );
    SLHA::SameNameBlockSet*
    blockPointer( slhaParser.getBlockAsStrings( blockName ) );
    if( NULL == blockPointer )
    {
      return returnString;
    }
    SLHA::BlockClass::BaseStringBlock const&
    blockAsStrings( (*blockPointer)[ 0 ] );
    BOL::VectorlikeArray< std::string > blockLine;

    for( int whichLine( blockAsStrings.getNumberOfBodyLines() );
         0 < whichLine;
         --whichLine )
    {
      blockLine.clearEntries();
      BOL::StringParser::parseByChar( blockAsStrings[ whichLine ].first,
                                      blockLine,
                                      BOL::StringParser::whitespaceChars );
      bool indicesMatch( false );
      if( blockLine.getSize() > (int)(indicesVector.size()) )
      {
        indicesMatch = true;
        for( unsigned int whichIndex( 0 );
             indicesVector.size() > whichIndex;
             ++whichIndex )
        {
          if( indicesVector[ whichIndex ]
              != BOL::StringParser::stringToInt( blockLine[ whichIndex ] ) )
          {
            indicesMatch = false;
            break;
          }
        }
      }
      if( indicesMatch )
      {
        for( int whichReturnWord( indicesVector.size() );
             blockLine.getSize() > whichReturnWord;
             ++whichReturnWord )
        {
          returnString.append( blockLine[ whichReturnWord ] );
          if( blockLine.getLastIndex() > whichReturnWord )
          {
            returnString.append( "   " );
          }
        }
        break;
      }
    }
    return returnString;
  }

} /* namespace LHPC */
