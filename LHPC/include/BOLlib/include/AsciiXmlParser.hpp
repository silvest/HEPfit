/*
 * AsciiXmlParser.hpp
 *
 *  Created on: Jan 23, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef ASCIIXMLPARSER_HPP_
#define ASCIIXMLPARSER_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include "StringParser.hpp"

namespace BOL
{
  /* this parses out blocks of ASCII text from a string between XML opening &
   * closing tags, & returns the text between the tags (without interpreting
   * it further as XML).
   */
  class AsciiXmlParser
  {
  public:
    AsciiXmlParser( bool const isVerbose = false );
    ~AsciiXmlParser();

    bool
    openRootElementOfFile( std::string const& fileName );
    /* this loads the file with the given name into the internal ifstream for
     * parsing, then opens the root element. if there was a problem loading the
     * file, or no root element could be opened, false is returned. this closes
     * the previously-loaded file, if it was open.
     */
    bool
    readAllOfRootElementOfFile( std::string const& fileName );
     /* this loads the file with the given name into the internal ifstream for
      * parsing, then reads in all the root element. if there was a problem
      * loading the file, or no root element could be opened, false is returned.
      * this closes the previously-loaded file, if it was open.
      */
    bool
    loadString( std::string const stringToParse );
    /* this loads stringToParse into the internal stringstream for parsing. if
      * there was a problem loading the file, false is returned. if there was a
      * file open, it is now closed.
      */
    void
    closeFile();
    // this closes fileToParse, if it was open.
    void
    returnToBeginningOfText();
    // this sets whichever stream textToParse is pointing at to start reading
    // again from its start.
    bool
    readNextElement();
    /* this reads in the entire next XML element. false is returned if the end
     * of the text was reached before a new XML element could be found (or if
     * no end for the next element could be found).
     */
    std::string const&
    getCurrentElementName() const;
    bool
    currentElementNameMatches( std::string const& comparisonString ) const;
    std::string const&
    getCurrentElementContent() const;
    std::string
    getTrimmedCurrentElementContent() const;
    // this returns a string that is a copy of the current XML, element, with
    // whitespace & newline characters trimmed from the front & back.
    std::string const&
    getCurrentOpeningTagAsFound() const;
    std::map< std::string, std::string > const&
    getCurrentElementAttributes();
    std::pair< int, int > const&
    getLineRangeOfCurrentElement() const;
    /* this returns a pair of ints where the 1st is 1 plus the number of
     * newline characters read in before the '<' of the current element's
     * opening tag was found, & the 2nd is the 1st plus the number of newline
     * characters read in before the '>' of the current element's closing tag
     * was found.
     */


  protected:
    static char const markupOpener;
    static char const markupCloser;
    static char const tagCloser;
    static std::string const allowedXmlWhitespaceChars;
    static std::string const allowedXmlQuoteChars;
    static std::pair< std::string, std::string > const commentDelimiter;
    static std::pair< std::string, std::string > const piDelimiter;
    static std::pair< std::string, std::string > const doctypeDelimiter;
    static std::pair< std::string, std::string > const cdataDelimiter;

    bool const isVerbose;
    std::ifstream fileParsingStream;
    std::istringstream stringParsingStream;
    std::istream* textStream;
    std::string rootTag;
    std::map< std::string, std::string > rootAttributeMap;
    std::pair< int, int > rootLineRange;
    std::string elementName;
    std::string fullElementContentAsFound;
    // this includes any whitespace & child elements (as text); this does not
    // include the opening or closing tags (& so does not include attributes).
    std::string fullOpeningTagAsFound;
    // this includes any attributes; this does not include '<' or '>'.
    std::map< std::string, std::string > elementAttributeMap;
    std::pair< int, int > elementLineRange;

    std::string markupString;
    bool streamIsGood;
    char currentChar;
    int readNewlines;
    int newlinesBeforeMarkup;
    size_t parseStart;
    size_t parseEnd;
    size_t previousLength;
    char currentQuoteChar;
    std::pair< std::string, std::string > currentAttribute;
    std::string closingTag;
    int closingTagsToFind;
    std::string currentTagName;

    void
    resetContent();
    // this sets the various recording data to the values they should have
    // before reading in some text.
    bool
    readCharacter();
    /* this reads the next character from textStream into currentChar &
     * increments readNewlines if currentChar was '\n', returning
     * textStream->good().
     */
    bool
    recordTo( std::string& recordingString,
              char const endChar );
    /* this records characters from textStream by appending them to
     * recordingStream, up to the 1st instance of endChar. it returns false if
     * no further instances of endChar were found in textStream.
     */
    bool
    closeMarkup( size_t const startPosition = 0 );
    /* this records characters from textStream by appending them to
     * markupString, up to the 1st instance of markupCloser that is not
     * enclosed in quotes, but only looking for quote characters from
     * startPosition onwards. an exception is made if the markup was a comment:
     * if markupString begins with commentDelimiter.first, all characters up to
     * the next found ( commentDelimiter.first + markupCloser ) are discarded
     * & markupString is emptied, then true is returned.
     */
    bool
    discardToNextMarkup();
    /* this skips characters up to the next XML markup opening character, then
     * the characters up to the next markup closing character are stored in
     * markupString. this also records the number of newline characters up to
     * this markup.
     */
    bool
    discardToNextTag();
    /* this calls discardToNextMarkup() & ignoreDelimited() on the non-tag
     * markups until a markup corresponding to a tag is found. this also
     * records the number of newline characters up to this tag.
     */
    bool
    recordToNextMarkup();
    /* this appends characters from textStream to fullElementContentAsFound up
     * to the next XML markup opening character, then the characters up to the
     * next markup closing character are stored in markupString.
     */
    bool
    recordToNextTag();
    /* this calls recordToNextMarkup() & recordDelimited() on the non-tag
     * markups (ignoring comments however) until a markup corresponding to a
     * tag is found.
     */
    bool
    compareMarkupStart( std::string const& comparisonString ) const;
    // this returns true if the start of "<" + markupString matches all of
    // comparisonString.
    bool
    compareMarkupEnd( std::string const& comparisonString ) const;
    // this returns true if the end of markupString + ">" matches all of
    // comparisonString.
    bool
    ignoreDelimited(
              std::pair< std::string, std::string > const& delimitingStrings );
    /* this checks to see if markupString begins with delimitingStrings.first,
     * & if so, ensures that markupString ends with delimitingStrings.second,
     * appending to markupString if necessary, then empties it. true is then
     * returned, unless the end of the text was reached before this could
     * happen.
     */
    bool
    recordDelimited(
              std::pair< std::string, std::string > const& delimitingStrings );
    /* this checks to see if markupString begins with delimitingStrings.first,
     * & if so, ensures that markupString ends with delimitingStrings.second,
     * appending to markupString if necessary, then appends
     * markupOpener + markupString + markupCloser to fullElementContentAsFound,
     * then empties markupString. true is then returned, unless the end of the
     * text was reached before this could happen.
     */
    bool
    eraseQuotedStringsInMarkup( size_t const startPosition );
    /* this erases any quoted text in markupString starting from startPosition,
     * extending markupString from textStream to the next unquoted
     * markupCloser.
     */
    void
    recordTagTo( std::string& recordingString );
    // this appends markupOpener + markupString + markupCloser to
    // recordingString.
    bool
    parseTagName( std::string& nameDestination );
    // this parses markupString as the tag of an XML element, leaving parseEnd
    // as the end of the name string.
    bool
    parseAttributes();
    /* this parses any attributes in markupString, assuming that parseEnd is at
     * the end of the tag's name. false is returned if a malformed attribute is
     * found.
     */
    bool
    recordToEndOfElement();
    /* this parses & stores any attributes & then stores the characters between
     * the opening tag & the corresponding closing tag in
     * fullElementContentAsFound, returning false if the end of the text was
     * reached before finding the closing tag.
     */
    bool
    skipPrologAndOpenRootElement();
  };





  inline bool
  AsciiXmlParser::readAllOfRootElementOfFile( std::string const& fileName )
  /* this loads the file with the given name into the internal ifstream for
   * parsing, then reads in all the root element. if there was a problem
   * loading the file, or no root element could be opened, false is returned.
   * this closes the previously-loaded file, if it was open.
   */
  {
    streamIsGood = ( openRootElementOfFile( fileName )
                     &&
                     recordToEndOfElement() );
    if( !streamIsGood
        &&
        isVerbose )
    {
      std::cout
      << std::endl
      << "BOL::error! AsciiXmlParser::readAllOfRootElementOfFile( " << fileName
      << " ) could not find a well-formed root element!";
      std::cout << std::endl;
    }
    rootLineRange.second = readNewlines;
    return streamIsGood;
  }

  inline bool
  AsciiXmlParser::loadString( std::string const stringToParse )
  /* this loads stringToParse into the internal stringstream for parsing. if
    * there was a problem loading the file, false is returned. if there was a
    * file open, it is now closed.
    */
  {
    resetContent();
    stringParsingStream.clear();
    stringParsingStream.str( stringToParse );
    textStream = &stringParsingStream;
    return stringParsingStream.good();
  }

  inline void
  AsciiXmlParser::closeFile()
  // this closes fileToParse, if it was open.
  {
    if( fileParsingStream.is_open() )
    {
      fileParsingStream.close();
    }
    fileParsingStream.clear();
    textStream = &stringParsingStream;
  }

  inline void
  AsciiXmlParser::returnToBeginningOfText()
  // this sets whichever stream textToParse is pointing at to start reading
  // again from its start.
  {
    resetContent();
    textStream->clear();
    textStream->seekg( std::ios::beg );
  }

  inline bool
  AsciiXmlParser::readNextElement()
  /* this reads in the entire next XML element. false is returned if the end
   * of the text was reached before a new XML element could be found (or if
   * no end for the next element could be found).
   */
  {
    return ( discardToNextTag()
             &&
             parseTagName( elementName )
             &&
             parseAttributes()
             &&
             recordToEndOfElement() );
  }

  inline std::string const&
  AsciiXmlParser::getCurrentElementName() const
  {
    return elementName;
  }

  inline bool
  AsciiXmlParser::currentElementNameMatches(
                                    std::string const& comparisonString ) const
  {
    return ( 0 == elementName.compare( comparisonString ) );
  }

  inline std::string const&
  AsciiXmlParser::getCurrentElementContent() const
  {
    return fullElementContentAsFound;
  }

  inline std::string
  AsciiXmlParser::getTrimmedCurrentElementContent() const
  {
    return StringParser::trimFromFrontAndBack( fullElementContentAsFound,
                                               allowedXmlWhitespaceChars );
  }

  inline std::string const&
  AsciiXmlParser::getCurrentOpeningTagAsFound() const
  {
    return fullOpeningTagAsFound;
  }

  inline std::map< std::string, std::string > const&
  AsciiXmlParser::getCurrentElementAttributes()
  {
    return elementAttributeMap;
  }

  inline std::pair< int, int > const&
  AsciiXmlParser::getLineRangeOfCurrentElement() const
  /* this returns a pair of ints where the 1st is 1 plus the number of
   * newline characters read in before the '<' of the current element's
   * opening tag was found, & the 2nd is the 1st plus the number of newline
   * characters read in before the '>' of the current element's closing tag
   * was found.
   */
  {
    return elementLineRange;
  }

  inline void
  AsciiXmlParser::resetContent()
  // this sets the various recording data to the values they should have
  // before reading in some text.
  {
    rootTag.assign( "" );
    rootAttributeMap.clear();
    rootLineRange.first = -1;
    rootLineRange.second = -1;
    elementName.assign( "" );
    fullElementContentAsFound.assign( "" );
    fullOpeningTagAsFound.assign( "" );
    elementAttributeMap.clear();
    elementLineRange.first = -1;
    elementLineRange.second = -1;
    markupString.assign( "" );
    readNewlines = 0;
  }

  inline bool
  AsciiXmlParser::readCharacter()
  /* this reads the next character from textStream into currentChar &
   * increments readNewlines if currentChar was '\n', returning
   * textStream->good().
   */
  {
    streamIsGood = textStream->get( currentChar ).good();
    if( streamIsGood
        &&
        ( '\n' == currentChar ) )
    {
      ++readNewlines;
    }
    return streamIsGood;
  }

  inline bool
  AsciiXmlParser::recordTo( std::string& recordingString,
                            char const endChar )
  /* this records characters from textStream by appending them to
   * recordingStream, up to the 1st instance of endChar. it returns false if
   * no further instances of endChar were found in textStream.
   */
  {
    while( readCharacter()
           &&
           endChar != currentChar )
    {
      recordingString.append( 1,
                              currentChar );
    }
    return streamIsGood;
  }

  inline bool
  AsciiXmlParser::discardToNextMarkup()
  /* this skips characters up to the next XML markup opening character, then
   * the characters up to the next markup closing character are stored in
   * markupString. this also records the number of newline characters up to
   * this markup.
   */
  {
    while( readCharacter()
           &&
           ( markupOpener != currentChar ) )
    {
      // the conditional does the work of reading up to the next tag.
    }
    newlinesBeforeMarkup = readNewlines;
    return ( streamIsGood
             &&
             closeMarkup() );
  }

  inline bool
  AsciiXmlParser::discardToNextTag()
  // this calls discardToNextMarkup() & ignoreDelimited() on the non-tag
  // markups until a markup corresponding to a tag is found.
  {
    markupString.assign( "" );
    while( markupString.empty() )
    {
      streamIsGood = ( discardToNextMarkup()
                       &&
                       ignoreDelimited( commentDelimiter )
                       &&
                       ignoreDelimited( piDelimiter )
                       &&
                       ignoreDelimited( cdataDelimiter ) );
      if( !streamIsGood )
      {
        return false;
      }
    }
    return true;
  }

  inline bool
  AsciiXmlParser::recordToNextMarkup()
  /* this appends characters from textStream to fullElementContentAsFound up
   * to the next XML markup opening character, then the characters up to the
   * next markup closing character are stored in markupString.
   */
  {
    return ( recordTo( fullElementContentAsFound,
                       markupOpener )
             &&
             closeMarkup() );
  }

  inline bool
  AsciiXmlParser::recordToNextTag()
  /* this calls recordToNextMarkup() & recordDelimited() on the non-tag
   * markups (ignoring comments however) until a markup corresponding to a tag
   * is found.
   */
  {
    markupString.assign( "" );
    while( markupString.empty() )
    {
      streamIsGood = ( recordToNextMarkup()
                       &&
                       recordDelimited( piDelimiter )
                       &&
                       recordDelimited( cdataDelimiter ) );
      if( !streamIsGood )
      {
        return false;
      }
    }
    return true;
  }

  inline bool
  AsciiXmlParser::compareMarkupStart(
                                    std::string const& comparisonString ) const
  // this returns true if the start of "<" + markupString matches all of
  // comparisonString.
  {
    return ( 0 == comparisonString.compare( 1,
                                            ( comparisonString.size() - 1 ),
                                            markupString,
                                            0,
                                           ( comparisonString.size() - 1 ) ) );
  }

  inline bool
  AsciiXmlParser::compareMarkupEnd( std::string const& comparisonString ) const
  // this returns true if the end of markupString + ">" matches all of
  // comparisonString.
  {
    return ( 0 == comparisonString.compare( 0,
                                            ( comparisonString.size() - 1 ),
                                            markupString,
                        ( markupString.size() - comparisonString.size() + 1 ),
                                           ( comparisonString.size() - 1 ) ) );
  }

  inline void
  AsciiXmlParser::recordTagTo( std::string& recordingString )
  // this appends markupOpener + markupString + markupCloser to
  // recordingString.
  {
    recordingString.append( 1,
                            markupOpener );
    recordingString.append( markupString );
    recordingString.append( 1,
                            markupCloser );
  }

  inline bool
  AsciiXmlParser::parseTagName( std::string& nameDestination )
  // this parses markupString as the tag of an XML element, leaving parseEnd
  // as the end of the name string.
  {
    parseStart = markupString.find_first_not_of( allowedXmlWhitespaceChars );
    if( std::string::npos == parseStart )
    {
      // empty markup is a sign of malformed XML:
      if( isVerbose )
      {
        std::cout
        << std::endl
        << "BOL::error! AsciiXmlParser::parseOpeningTag(...) found empty"
        << " markup!";
        std::cout << std::endl;
      }
      return false;
    }
    parseEnd = markupString.find_first_of( allowedXmlWhitespaceChars,
                                           parseStart );
    nameDestination.assign( markupString.substr( parseStart,
                                                 ( parseEnd - parseStart ) ) );
    return true;
  }

}

#endif /* ASCIIXMLPARSER_HPP_ */
