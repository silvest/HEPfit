/*
 * CommentedTextParser.hpp
 *
 *  Created on: Jan 22, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef COMMENTEDTEXTPARSER_HPP_
#define COMMENTEDTEXTPARSER_HPP_

#include <fstream>
#include "StringParser.hpp"
#include "VectorlikeArray.hpp"

namespace BOL
{
  /* this class breaks up lines of text into pairs of "line before comment
   * marking character" with "rest of line (including comment marker)". it can
   * take in a block of text as a string (including newline characters) or it
   * can open a file to read in text.
   */
  class CommentedTextParser
  {
  public:
    CommentedTextParser( std::string const& commentMarker,
                         bool const isVerbose );
    CommentedTextParser(
                      VectorlikeArray< std::string > const& commentMarkerSet,
                         bool const isVerbose );
    virtual
    ~CommentedTextParser();

    VectorlikeArray< std::pair< std::string, std::string > > const&
    parseString( std::string const& textToParse );
    /* this breaks up textToParse into lines based on '\n' & '\r', then breaks
     * up each line based on commentMarker, stores everything in parsedText &
     * returns a reference to parsedText. N.B. entirely empty lines are skipped
     * over (i.e. any uninterrupted sequence of '\n' chars is treated as a
     * single newline)!
     */
    bool
    openFile( std::string const& fileName );
    /* this opens a file based on the given name. if this CommentedTextParser
     * instance was already reading a file, the old file is closed. this
     * returns true if the file was opened successfully.
     */
    bool
    atEndOfFile();
    // this returns true if all the lines of the file have been read, either by
    // parseNextLineOfFile() or readNextNonEmptyLineOfFileWithoutComment(...).
    CommentedTextParser&
    closeFile();
    bool
    parseNextLineOfFile( std::string& stringForData,
                         std::string& stringForComment );
    /* this reads the next line of the file (i.e. up to the next newline char),
     * splits it into a bit before the comment marker & the rest (including the
     * comment marker characters), & puts the 2 parts into the given strings
     * (the comment in the 2nd part), & returns true if the line was
     * successfully parsed.
     */
    bool
    readNextNonEmptyLineOfFileWithoutComment( std::string& stringToFill );
    /* this reads in lines into stringToFill until, after the following
     * manipulations, stringToFill is not empty, or the end of inputFile is
     * reached. the manipulations are that any comment it has is removed, &
     * then any leading & trailing whitespace characters are removed. if the
     * end of the file is reached before stringToFill is not empty after the
     * manipulations, stringToFill is left empty & false is returned.
     */
    bool
    readJustNextValidLine( std::string& stringToFill )
    { return readNextNonEmptyLineOfFileWithoutComment( stringToFill ); }
    // this is just to shorten the name without breaking backwards
    // compatibility.


  protected:
    static std::string const trimmingChars;

    bool const isVerbose;
    VectorlikeArray< std::string > commentMarkerSet;
    VectorlikeArray< std::pair< std::string, std::string > > parsedText;
    std::pair< std::string, std::string > parsedLine;
    VectorlikeArray< std::string > textAsLines;
    std::string lineBeingRead;
    bool linesOfFileRemain;
    // this is used as a safer bool than inputFile.eof() because inputFile is
    // not necessarily open or otherwise.
    std::ifstream inputFile;
    // this is the file being read in.

    bool
    readInNextLine();
    /* this reads in the next line of inputFile into lineBeingRead, unless
     * inputFile has been read to its end, & if reading in this line brings
     * inputFile to its EOF, this function sets parsedLine & linesOfFileRemain
     * appropriately, as well as closing inputFile. it returns true if a line
     * was successfully read in.
     */
    void
    reportStateOfFile();
    // this reports on whether this instance will no longer read in any more
    // lines because the end of the file was reached or because of some error.
  };



  inline VectorlikeArray< std::pair< std::string, std::string > > const&
  CommentedTextParser::parseString( std::string const& textToParse )
  /* this breaks up textToParse into lines based on '\n' & '\r', then breaks
   * up each line based on commentMarker, stores everything in parsedText &
   * returns a reference to parsedText. N.B. entirely empty lines are skipped
   * over (i.e. any uninterrupted sequence of '\n' chars is treated as a single
   * newline)!
   */
  {
    BOL::StringParser::parseByChar( textToParse,
                                    textAsLines,
                                    BOL::StringParser::newlineChars );
    // now textAsLines is textToParse broken up into strings corresponding to
    // each line of textToParse. parsedText has to fit all of textAsLines:
    parsedText.setSize( textAsLines.getSize() );
    for( int lineIndex( textAsLines.getLastIndex() );
         0 <= lineIndex;
         --lineIndex )
    {
      parsedText[ lineIndex ].first.assign(
                 BOL::StringParser::substringToFirst( textAsLines[ lineIndex ],
                                                      commentMarkerSet,
                                         &(parsedText[ lineIndex ].second) ) );
      /* this puts all of textAsLines[ lineIndex ] before the comment marker
       * into parsedText[ lineIndex ].first & the remainder into
       * parsedText[ lineIndex ].second (including the comment characters).
       */
    }
    return parsedText;
  }

  inline bool
  CommentedTextParser::openFile( std::string const& fileName )
  /* this opens a file based on the given name. if this CommentedTextParser
   * instance was already reading a file, the old file is closed. this returns
   * true if the file was opened successfully.
   */
  {
    closeFile();
    // this sets notYetAtEndOfFile to false & also closes inputFile if it was
    // open.
    inputFile.open( fileName.c_str() );
    if( inputFile.is_open()
        &&
        !(inputFile.eof()) )
    {
      linesOfFileRemain = true;
    }
    else
    {
      std::cout
      << std::endl
      << "BOL::error! CommentedTextParser tried to open "
      << fileName << " but could not (or " << fileName << " is empty)!";
      std::cout << std::endl;
      linesOfFileRemain = false;
    }
    return linesOfFileRemain;
  }

  inline bool
  CommentedTextParser::atEndOfFile()
  // this returns true if all the lines of the file have been read, either by
  // parseNextLineOfFile() or readNextNonEmptyLineOfFileWithoutComment(...).
  {
    return !linesOfFileRemain;
  }

  inline CommentedTextParser&
  CommentedTextParser::closeFile()
  {
    linesOfFileRemain = false;
    if( inputFile.is_open() )
    {
      inputFile.close();
    }
    return *this;
  }

  inline bool
  CommentedTextParser::parseNextLineOfFile( std::string& stringForData,
                                            std::string& stringForComment )
  /* this reads the next line of the file (i.e. up to the next newline char),
   * splits it into a bit before the comment marker & the rest (including the
   * comment marker characters), & puts the 2 parts into the given strings
   * (the comment in the 2nd part), & returns true if the line was
   * successfully parsed.
   */
  {
    if( readInNextLine() )
    {
      stringForData.assign( BOL::StringParser::substringToFirst(
                                BOL::StringParser::trimFromBack( lineBeingRead,
                                             BOL::StringParser::newlineChars ),
                                                              commentMarkerSet,
                                                         &stringForComment ) );
      /* this puts all of the line just read before the comment marker into
       * stringForData & the remainder into stringForComment (including the
       * comment characters).
       */
      BOL::StringParser::trimFromBack( stringForData,
                                       BOL::StringParser::newlineChars );
      return true;
    }
    else
    {
      return false;
    }
  }

  inline bool
  CommentedTextParser::readNextNonEmptyLineOfFileWithoutComment(
                                                    std::string& stringToFill )
  /* this reads in lines into stringToFill until, after the following
   * manipulations, stringToFill is not empty, or the end of inputFile is
   * reached. the manipulations are that any comment it has is removed, &
   * then any leading & trailing whitespace characters are removed. if the
   * end of the file is reached before stringToFill is not empty after the
   * manipulations, stringToFill is left empty & false is returned.
   */
  {
    stringToFill.clear();
    while( stringToFill.empty()
           &&
           readInNextLine() )
      // the order of these conditionals is important: it has to check if it
      // needs to read in the next line before it tries to read it in!
    {
      stringToFill.assign( BOL::StringParser::substringToFirst( lineBeingRead,
                                                          commentMarkerSet ) );
      /* this puts all of the line just read before the comment marker into
       * stringToFill & throws away the remainder. now all trimmingChars chars
       * are trimmed from the front & back of stringToFill:
       */
      stringToFill.assign( BOL::StringParser::trimFromFrontAndBack(
                                                                  stringToFill,
                                                             trimmingChars ) );
    }
    // now either stringToFill is not empty, or inputFile ran out of lines.
    return !(stringToFill.empty());
  }

  inline bool
  CommentedTextParser::readInNextLine()
  /* this reads in the next line of inputFile into lineBeingRead, unless
   * inputFile has been read to its end, & if reading in this line brings
   * inputFile to its EOF, this function sets parsedLine & linesOfFileRemain
   * appropriately, as well as closing inputFile. it returns true if a line
   * was successfully read in.
   */
  {
    if( linesOfFileRemain )
    {
      std::getline( inputFile,
                    lineBeingRead );
      linesOfFileRemain = inputFile.good();
      if( !linesOfFileRemain
          &&
          isVerbose )
      {
        reportStateOfFile();
      }
      return true;
    }
    else
    {
      return false;
    }
  }

  inline void
  CommentedTextParser::reportStateOfFile()
  // this reports on whether this instance will no longer read in any more
  // lines because the end of the file was reached or because of some error.
  {
    if( inputFile.eof() )
      // if reading in this line brought inputFile to its end, note that.
    {
      std::cout
      << std::endl
      << "BOL::CommentedTextParser reached the end of the file.";
      std::cout << std::endl;
    }
    else
    {
      std::cout
      << std::endl
      << "BOL::CommentedTextParser cannot read in any more of the file after"
      << " a bad reading operation.";
      std::cout << std::endl;
    }
  }

}

#endif /* COMMENTEDTEXTPARSER_HPP_ */
