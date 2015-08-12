/*
 * ArgumentParser.cpp
 *
 *  Created on: Sep 13, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "ArgumentParser.hpp"

namespace BOL
{
  ArgumentParser::ArgumentParser( int argumentCount,
                                  char** argumentCharArrays,
                                  std::string const inputTag,
                                  std::string const fallbackInputFilename ) :
    argumentStrings( ( argumentCount - 1 ),
                     "" ),
    inputXmlParser()
  {
    inputXmlParser.loadString( "" );
    for( int whichArgument( 1 );
         argumentCount > whichArgument;
         ++whichArgument )
    {
      argumentStrings[ whichArgument - 1 ].assign(
                                         argumentCharArrays[ whichArgument ] );
    }
    if( !(inputTag.empty()) )
    {
      std::string inputFilename( fromTag( inputTag ) );
      if( inputFilename.empty() )
      {
        inputFilename.assign( fallbackInputFilename );
      }
      // at this point, if no input filename was provided either by the command
      // line arguments or by the constructor, I think that it's obvious that
      // the user didn't want to look for arguments in an XML file, so the
      // error message should not be shown for failing to open it.
      if( !(inputFilename.empty()) )
      {
        bool successfullyRead( inputXmlParser.readAllOfRootElementOfFile(
                                                                inputFilename )
                               &&
                               inputXmlParser.loadString(
                                 inputXmlParser.getCurrentElementContent() ) );
        if( !successfullyRead )
        {
          std::cout
          << std::endl
          << "BOL::ArgumentParser constructor failed to open root element of"
          << " \"" << inputFilename << "\"!";
          std::cout << std::endl;
        }
      }
    }
  }

  ArgumentParser::~ArgumentParser()
  {
    // does nothing.
  }

} /* namespace BOL */
