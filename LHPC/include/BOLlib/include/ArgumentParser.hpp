/*
 * ArgumentParser.hpp
 *
 *  Created on: Sep 13, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef ARGUMENTPARSER_HPP_
#define ARGUMENTPARSER_HPP_

#include <string>
#include <vector>
#include <fstream>
#include "AsciiXmlParser.hpp"
#include "StringParser.hpp"

namespace BOL
{
  // this is a class to search the command line arguments for given strings.
  class ArgumentParser
  {
  public:
    ArgumentParser( int argumentCount,
                    char** argumentCharArrays,
                    std::string const inputTag = "",
                    std::string const fallbackInputFilename = "" );
    ~ArgumentParser();

    std::string
    fromLiteral( std::string const& argumentName,
                 std::string const defaultReturnString = "" ) const;
    // this looks for argumentName as a substring (from the starting char) of
    // any of the elements of argumentStrings, & returns the rest of the 1st
    // matching string (or defaultReturnString if no match was found).
    std::string
    fromTag( std::string const& tagString,
             std::string const defaultReturnString = "" );
    // this calls fromLiteral( "--" + argumentName + "=" ), and if it gets an
    // empty string, it looks for an element with argumentName as an XML tag
    // from the input file (if the input file was readable). if there is still
    // no string to return, it returns defaultReturnString.


  protected:
    std::vector< std::string > argumentStrings;
    AsciiXmlParser inputXmlParser;
  };





  inline std::string
  ArgumentParser::fromLiteral( std::string const& argumentName,
                               std::string const defaultReturnString ) const
  // this looks for argumentName as a substring (from the starting char) of
  // any of the elements of argumentStrings, & returns the rest of the 1st
  // matching string (or defaultReturnValue if no match was found).
  {
    for( std::vector< std::string >::const_iterator
         whichArgument( argumentStrings.begin() );
         argumentStrings.end() > whichArgument;
         ++whichArgument )
    {
      if( 0 == whichArgument->compare( 0,
                                       argumentName.size(),
                                       argumentName ) )
      {
        return whichArgument->substr( argumentName.size() );
      }
    }
    return defaultReturnString;
  }

  inline std::string
  ArgumentParser::fromTag( std::string const& tagString,
                           std::string const defaultReturnString )
  // this calls fromLiteral( "--" + argumentName + "=" ), and if it gets an
  // empty string, it looks for an element with argumentName as an XML tag
  // from the input file (if the input file was readable). if there is still
  // no string to return, it returns defaultReturnString.
  {
    std::string argumentString( "--" );
    argumentString.append( tagString );
    argumentString.append( "=" );
    argumentString.assign( fromLiteral( argumentString ) );
    if( !(argumentString.empty()) )
    {
      return argumentString;
    }
    inputXmlParser.returnToBeginningOfText();
    while( inputXmlParser.readNextElement() )
    {
      if( inputXmlParser.currentElementNameMatches( tagString ) )
      {
        return StringParser::trimFromFrontAndBack(
                                     inputXmlParser.getCurrentElementContent(),
                                     StringParser::whitespaceAndNewlineChars );
      }
    }
    return defaultReturnString;
  }

} /* namespace BOL */
#endif /* ARGUMENTPARSER_HPP_ */
