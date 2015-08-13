/*
 * CommentedTextParser.cpp
 *
 *  Created on: Jan 22, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "CommentedTextParser.hpp"

namespace BOL
{
  std::string const CommentedTextParser::trimmingChars( " \t\r" );

  CommentedTextParser::CommentedTextParser( std::string const& commentMarker,
                                            bool const isVerbose ) :
      isVerbose( isVerbose ),
      commentMarkerSet( 1 ),
      parsedText(),
      parsedLine( "",
                  "" ),
      textAsLines(),
      lineBeingRead( "" ),
      linesOfFileRemain( false ),
      inputFile()
  {
    commentMarkerSet.getFront().assign( commentMarker );
  }

  CommentedTextParser::CommentedTextParser(
                        VectorlikeArray< std::string > const& commentMarkerSet,
                                            bool const isVerbose ) :
      isVerbose( isVerbose ),
      commentMarkerSet( commentMarkerSet ),
      parsedText(),
      parsedLine( "",
                  "" ),
      textAsLines(),
      lineBeingRead( "" ),
      linesOfFileRemain( false ),
      inputFile()
  {
    // just an initialization list.
  }

  CommentedTextParser::~CommentedTextParser()
  {
    if( inputFile.is_open() )
    {
      inputFile.close();
    }
  }

}
