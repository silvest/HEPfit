/*
 * InputFileToOutputFileClaimer.cpp
 *
 *  Created on: Aug 14, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "InputFileToOutputFileClaimer.hpp"

namespace BOL
{
  std::string const InputFileToOutputFileClaimer::claimedSuffix( ".claimed" );

  InputFileToOutputFileClaimer::InputFileToOutputFileClaimer() :
    matchedFilePairs(),
    fileNumber( -1 ),
    claimedFilename( "" ),
    copyCommand( "cp" ),
    removeCommand( "rm" )
  {
    // just an initialization list.
  }

  InputFileToOutputFileClaimer::InputFileToOutputFileClaimer(
                    std::vector< std::pair< std::string, std::string > > const&
                                                           matchedFilePairs ) :
    matchedFilePairs( matchedFilePairs ),
    fileNumber( -1 ),
    claimedFilename( "" ),
    copyCommand( "cp" ),
    removeCommand( "rm" )
  {
    // just an initialization list.
  }

  InputFileToOutputFileClaimer::~InputFileToOutputFileClaimer()
  {
    // does nothing.
  }


  bool
  InputFileToOutputFileClaimer::claimNextUnclaimedFile(
                                                std::string const claimString )
  // this increments fileNumber, & returns false if there is no input file
  // corresponding to it. otherwise, it updates currentOutputFilename &
  // sets claimedFilename as the same, then appends claimedFilename with
  // claimedSuffix. it then creates a file with claimedFilename as its name,
  // fills it with claimString, & returns true.
  {
    // the conditional advances through the list of filename pairs until either
    // the loop is broken by the early return of true because the next output
    // file was claimed, or the early return of false because the claimed file
    // was not openable.
    while( matchedFilePairs.size() > (unsigned int)(++fileNumber) )
    {
      if( !(UsefulStuff::fileExists( matchedFilePairs[ fileNumber ].second )) )
      {
        claimedFilename.assign( matchedFilePairs[ fileNumber ].second );
        claimedFilename.append( claimedSuffix );
        if( !(UsefulStuff::fileExists( claimedFilename )) )
        {
          std::ofstream outputFile( claimedFilename.c_str() );
          if( outputFile.good() )
          {
            outputFile << claimString << std::endl;
            outputFile.close();
            return true;
          }
          else
          {
            std::cout
            << std::endl
            << "BOL::error!"
            << " InputFileToOutputFileClaimer::claimNextUnclaimedFile( \""
            << claimString << "\" ) unable to create file \""
            << claimedFilename << "\"!";
            std::cout << std::endl;
            return false;
          }
        }
      }
    }
    // if we get to the end of the loop, all the output filenames are taken or
    // claimed. (I've decided that SESE is a bit redundant these days.)
    return false;
  }

} /* namespace BOL */
