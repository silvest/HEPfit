/*
 * InputFileToOutputFileClaimer.hpp
 *
 *  Created on: Aug 14, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef INPUTFILETOOUTPUTFILECLAIMER_HPP_
#define INPUTFILETOOUTPUTFILECLAIMER_HPP_

#include <string>
#include <fstream>
#include "UsefulStuff.hpp"
#include "StringParser.hpp"

namespace BOL
{
  class InputFileToOutputFileClaimer
  {
  public:
    static std::string const claimedSuffix;

    InputFileToOutputFileClaimer();
    InputFileToOutputFileClaimer(
                    std::vector< std::pair< std::string, std::string > > const&
                                                            matchedFilePairs );
    ~InputFileToOutputFileClaimer();

    void
    setCopyCommand( std::string const& copyCommand );
    void
    setRemoveCommand( std::string const& removeCommand );
    void
    addFilenamePair( std::string const& inputFilename,
                     std::string const outputFilename );
    std::vector< std::pair< std::string, std::string > >&
    getFilenamePairs();
    bool
    claimNextUnclaimedFile(
                   std::string const claimString = "claimed, but not ready!" );
    // this looks for the next output filename in matchedFilePairs that does
    // not match an existing file's name & also has no corresponding claimed
    // file. if there is no such filename, false is returned. otherwise it
    // creates the corresponding claimed file & returns true.
    int
    getLastClaimedFileIndex() const;
    // this returns the index in matchedFilePairs that corresponds to the last
    // file claimed by this InputFileToOutputFileClaimer.
    std::pair< std::string, std::string > const&
    getFilenamesForLastClaim() const;
    // this returns the filenames in matchedFilePairs that correspond to the
    // last file claimed by this InputFileToOutputFileClaimer.
    std::string const&
    getInputFilenameForLastClaim() const;
    // this returns the input filename in matchedFilePairs that corresponds to
    // the last file claimed by this InputFileToOutputFileClaimer.
    std::string const&
    getOutputFilenameForLastClaim() const;
    // this returns the output filename in matchedFilePairs that corresponds to
    // the last file claimed by this InputFileToOutputFileClaimer.
    void
    copyToClaimedOutputFile( std::string const& filenameToCopy );
    // this copies the file called filenameToCopy to currentOutputFilename.


  protected:
    std::vector< std::pair< std::string, std::string > > matchedFilePairs;
    int fileNumber;
    std::string claimedFilename;
    std::string copyCommand;
    std::string removeCommand;
  };





  inline void
  InputFileToOutputFileClaimer::setCopyCommand(
                                               std::string const& copyCommand )
  {
    this->copyCommand.assign( copyCommand );
  }

  inline void
  InputFileToOutputFileClaimer::setRemoveCommand(
                                             std::string const& removeCommand )
  {
    this->removeCommand.assign( removeCommand );
  }

  inline void
  InputFileToOutputFileClaimer::addFilenamePair(
                                              std::string const& inputFilename,
                                             std::string const outputFilename )
  {
    matchedFilePairs.push_back( std::pair< std::string, std::string >(
                                                                 inputFilename,
                                                            outputFilename ) );
  }

  inline std::vector< std::pair< std::string, std::string > >&
  InputFileToOutputFileClaimer::getFilenamePairs()
  {
    return matchedFilePairs;
  }

  inline int
  InputFileToOutputFileClaimer::getLastClaimedFileIndex() const
  {
    return fileNumber;
  }

  inline std::pair< std::string, std::string > const&
  InputFileToOutputFileClaimer::getFilenamesForLastClaim() const
  // this returns the filenames in matchedFilePairs that correspond to the last
  // file claimed by this InputFileToOutputFileClaimer.
  {
    return matchedFilePairs[ fileNumber ];
  }

  inline std::string const&
  InputFileToOutputFileClaimer::getInputFilenameForLastClaim() const
  // this returns the input filename in matchedFilePairs that corresponds to
  // the last file claimed by this InputFileToOutputFileClaimer.
  {
    return matchedFilePairs[ fileNumber ].first;
  }

  inline std::string const&
  InputFileToOutputFileClaimer::getOutputFilenameForLastClaim() const
  // this returns the output filename in matchedFilePairs that corresponds to
  // the last file claimed by this InputFileToOutputFileClaimer.
  {
    return matchedFilePairs[ fileNumber ].second;
  }

  inline void
  InputFileToOutputFileClaimer::copyToClaimedOutputFile(
                                            std::string const& filenameToCopy )
  // this copies the file called filenameToCopy to currentOutputFilename.
  {
    std::string systemCommand( copyCommand );
    systemCommand.append( " " );
    systemCommand.append( filenameToCopy );
    systemCommand.append( " " );
    systemCommand.append( matchedFilePairs[ fileNumber ].second );
    system( systemCommand.c_str() );
    systemCommand.assign( removeCommand );
    systemCommand.append( " " );
    systemCommand.append( claimedFilename );
    system( systemCommand.c_str() );
  }

} /* namespace BOL */
#endif /* INPUTFILETOOUTPUTFILECLAIMER_HPP_ */
