/*
 * FilePlaceholderManager.cpp
 *
 *  Created on: Oct 10, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../include/FilePlaceholderManager.hpp"
#include <cstdio>
#include <iostream>

namespace BOL
{

  FilePlaceholderManager::FilePlaceholderManager(
                                                 std::string const inputSuffix,
                                           std::string const placeholderSuffix,
                                              std::string const outputSuffix ):
    inputSuffix( inputSuffix ),
    placeholderSuffix( placeholderSuffix ),
    outputSuffix( outputSuffix ),
    filenameTriples(),
    whichTriple(),
    currentTriple(),
    currentSuffix( "" ),
    lastPlaceholder( "" )
  {
    // just an initialization list.
  }

  FilePlaceholderManager::~FilePlaceholderManager()
  {
    // does nothing.
  }


  void
  FilePlaceholderManager::prepareFilenames( std::string const& inputDirectory,
                                       std::string const& placeholderDirectory,
                                           std::string const& outputDirectory )
  /* this takes all the names of all the files in the directory given by
   * inputDirectory which end in inputSuffixToStrip, and places each in a
   * FilenameTriple with the placeholder filename given by replacing
   * inputSuffixToStrip with placeholderSuffixToAdd, replacing the directory
   * path with placeholderDirectory if placeholderDirectory is not empty, and
   * the same for the output file using outputSuffixToAdd and outputDirectory
   * if not empty.
   */
  {
    std::vector< std::string > validBaseFilenames;
    DIR* directoryPointer( opendir( inputDirectory.c_str() ) );
    if( NULL == directoryPointer )
    {
      throw std::runtime_error( "Could not read directory!" );
    }
    std::string currentFile;
    struct dirent* structPointer( readdir( directoryPointer ) );
    while( NULL != structPointer )
    {
      currentFile.assign( structPointer->d_name );
      structPointer = readdir( directoryPointer );
      if( ( 0 == currentFile.compare( "." ) )
          ||
          ( 0 == currentFile.compare( ".." ) )
          ||
          !( inputSuffix.size() < currentFile.size() ) )
      {
        continue;
      }
      currentSuffix.assign( currentFile.substr( currentFile.size()
                                                   - inputSuffix.size() ) );
      if( 0 == currentSuffix.compare( inputSuffix ) )
      {
        validBaseFilenames.push_back( currentFile.substr( 0,
                               ( currentFile.size() - inputSuffix.size() ) ) );
      }
    }
    prepareFilenames( validBaseFilenames,
                      inputDirectory,
                      placeholderDirectory,
                      outputDirectory );
  }

  void
  FilePlaceholderManager::prepareFilenames(
                               std::vector< std::string > const& baseFilenames,
                                            std::string const inputDirectory,
                                            std::string placeholderDirectory,
                                            std::string outputDirectory )
  /* this takes all the names of all the files in baseFilenames and places
   * each in a FilenameTriple with inputSuffixToAdd appended as the input
   * filename, with placeholderSuffixToAdd appended as the placeholder
   * filename, replacing the directory path with placeholderDirectory if
   * placeholderDirectory is not empty, and with outputSuffixToAdd appended
   * as the output file, replacing the directory path with outputDirectory if
   * outputDirectory is not empty.
   */
  {
    filenameTriples.clear();
    if( placeholderDirectory.empty() )
    {
      placeholderDirectory.assign( inputDirectory );
    }
    else
    {
      UsefulStuff::runSystemCommand( "mkdir -p " + placeholderDirectory );
    }
    if( outputDirectory.empty() )
    {
      outputDirectory.assign( inputDirectory );
    }
    else
    {
      UsefulStuff::runSystemCommand( "mkdir -p " + outputDirectory );
    }
    for( std::vector< std::string >::const_iterator
         whichBaseFilename( baseFilenames.begin() );
         baseFilenames.end() > whichBaseFilename;
         ++whichBaseFilename )
    {
      currentTriple.inputFile.assign( inputDirectory + "/"
                                      + (*whichBaseFilename) + inputSuffix );
      currentTriple.placeholderFile.assign( placeholderDirectory + "/"
                                  + (*whichBaseFilename) + placeholderSuffix );
      currentTriple.outputFile.assign( outputDirectory + "/"
                                       + (*whichBaseFilename) + outputSuffix );
      filenameTriples.push_back( currentTriple );
    }
    whichTriple = filenameTriples.begin();
  }

  bool
  FilePlaceholderManager::holdNextPlace( bool const deleteLastPlaceholder )
  /* this looks to find the first FilenameTriple in filenameTriples which
   * has both a placeholder filename and an output filename which do not yet
   * exist in the file system, and returns true if there was such a triple.
   * if deleteLastPlaceholder is true, the previous placeholder is also
   * deleted, obviously.
   */
  {
    while( ( filenameTriples.end() > whichTriple )
           &&
           ( UsefulStuff::fileExists( whichTriple->outputFile )
             ||
             UsefulStuff::fileExists( whichTriple->placeholderFile )
             ||
             !(UsefulStuff::fileExists( whichTriple->inputFile )) ) )
    {
      ++whichTriple;
    }
    if( deleteLastPlaceholder )
    {
      UsefulStuff::runSystemCommand( "rm " + lastPlaceholder );
    }
    if( !( filenameTriples.end() > whichTriple ) )
    {
      return false;
    }
    lastPlaceholder.assign( whichTriple->placeholderFile );
    std::ofstream placeholderStream( lastPlaceholder.c_str() );
    placeholderStream << "placeholder" << std::endl;
    placeholderStream.close();
    return true;
  }

} /* namespace BOL */
