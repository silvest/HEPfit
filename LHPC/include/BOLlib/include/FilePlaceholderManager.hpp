/*
 * FilePlaceholderManager.hpp
 *
 *  Created on: Oct 10, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef FILEPLACEHOLDERMANAGER_HPP_
#define FILEPLACEHOLDERMANAGER_HPP_

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <dirent.h>
#include "UsefulStuff.hpp"

namespace BOL
{
  // this struct is just a convenient grouping of data for the following class.
  struct FilenameTriple
  {
    std::string inputFile;
    std::string placeholderFile;
    std::string outputFile;
  };

  /* this class manages triples of filenames with the intention that each
   * triple consists of a name for an input file, a name for a placeholder
   * file to indicate that the output file is currently being worked on, and a
   * name for an output file.
   */
  class FilePlaceholderManager
  {
  public:
    FilePlaceholderManager( std::string const inputSuffix = "",
                            std::string const placeholderSuffix = "",
                            std::string const outputSuffix = "" );
    ~FilePlaceholderManager();

    void
    prepareFilenames( std::string const& inputDirectory,
                      std::string const& placeholderDirectory = "",
                      std::string const& outputDirectory = "" );
    /* this takes all the names of all the files in the directory given by
     * inputDirectory which end in inputSuffixToStrip, and places each in a
     * FilenameTriple with the placeholder filename given by replacing
     * inputSuffixToStrip with placeholderSuffixToAdd, replacing the directory
     * path with placeholderDirectory if placeholderDirectory is not empty, and
     * the same for the output file using outputSuffixToAdd and outputDirectory
     * if not empty.
     */
    void
    prepareFilenames( std::vector< std::string > const& baseFilenames,
                      std::string const inputDirectory = "",
                      std::string placeholderDirectory = "",
                      std::string outputDirectory = "" );
    /* this takes all the names of all the files in baseFilenames and places
     * each in a FilenameTriple with inputSuffixToAdd appended as the input
     * filename, with placeholderSuffixToAdd appended as the placeholder
     * filename, replacing the directory path with placeholderDirectory if
     * placeholderDirectory is not empty, and with outputSuffixToAdd appended
     * as the output file, replacing the directory path with outputDirectory if
     * outputDirectory is not empty.
     */
    bool
    holdNextPlace( bool const deleteLastPlaceholder = true );
    /* this looks to find the first FilenameTriple in filenameTriples which
     * has both a placeholder filename and an output filename which do not yet
     * exist in the file system, and returns true if there was such a triple.
     * if deleteLastPlaceholder is true, the previous placeholder is also
     * deleted, obviously.
     */
    std::string const&
    getInput() const;
    std::string const&
    getPlaceholder() const;
    std::string const&
    getOutput() const;


  protected:
    std::string const inputSuffix;
    std::string const placeholderSuffix;
    std::string const outputSuffix;
    std::vector< FilenameTriple > filenameTriples;
    std::vector< FilenameTriple >::const_iterator whichTriple;
    FilenameTriple currentTriple;
    std::string currentSuffix;
    std::string lastPlaceholder;
  };



  inline std::string const&
  FilePlaceholderManager::getInput() const
  {
    return whichTriple->inputFile;
  }

  inline std::string const&
  FilePlaceholderManager::getPlaceholder() const
  {
    return whichTriple->placeholderFile;
  }

  inline std::string const&
  FilePlaceholderManager::getOutput() const
  {
    return whichTriple->outputFile;
  }

} /* namespace BOL */
#endif /* FILEPLACEHOLDERMANAGER_HPP_ */
