/*
 * AsciiXmlStringParser.cpp
 *
 *  Created on: Jan 23, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "AsciiXmlParser.hpp"

namespace BOL
{
  char const AsciiXmlParser::markupOpener( '<' );
  char const AsciiXmlParser::markupCloser( '>' );
  char const AsciiXmlParser::tagCloser( '/' );
  std::string const AsciiXmlParser::allowedXmlWhitespaceChars( " \t\r\n" );
  std::string const AsciiXmlParser::allowedXmlQuoteChars( "\'\"" );
  std::pair< std::string, std::string > const
  AsciiXmlParser::commentDelimiter( "<!--",
                                    "-->" );
  std::pair< std::string, std::string > const
  AsciiXmlParser::piDelimiter( "<?",
                               "?>" );
  std::pair< std::string, std::string > const
  AsciiXmlParser::doctypeDelimiter( "<!DOCTYPE",
                                    ">" );
  std::pair< std::string, std::string > const
  AsciiXmlParser::cdataDelimiter( "<![CDATA[",
                                  "]]>" );


  AsciiXmlParser::AsciiXmlParser( bool const isVerbose ) :
      isVerbose( isVerbose ),
      fileParsingStream(),
      stringParsingStream(),
      textStream( &stringParsingStream ),
      rootTag( "" ),
      rootAttributeMap(),
      rootLineRange( -1,
                     -1 ),
      elementName( "" ),
      fullElementContentAsFound( "" ),
      fullOpeningTagAsFound( "" ),
      elementAttributeMap(),
      elementLineRange( -1,
                        -1 ),
      markupString( "" ),
      streamIsGood( false ),
      currentChar( ' ' ),
      readNewlines( 0 ),
      newlinesBeforeMarkup( 0 ),
      parseStart( 0 ),
      parseEnd( 0 ),
      previousLength( 0 ),
      currentQuoteChar( '\'' ),
      currentAttribute( "",
                        "" ),
      closingTag( "" ),
      closingTagsToFind( 0 ),
      currentTagName( "" )
  {
    // just an initialization list.
  }

  AsciiXmlParser::~AsciiXmlParser()
  {
    if( fileParsingStream.is_open() )
    {
      fileParsingStream.close();
    }
  }


  bool
  AsciiXmlParser::openRootElementOfFile( std::string const& fileName )
  /* this loads the file with the given name into the internal ifstream for
   * parsing, then opens the root element. if there was a problem loading the
   * file, or no root element could be opened, false is returned. this closes
   * the previously-loaded file, if it was open.
   */
  {
    closeFile();
    fileParsingStream.open( fileName.c_str() );
    if( !(fileParsingStream.good()) )
    {
      if( isVerbose )
      {
        std::cout
        << std::endl
        << "BOL::error! AsciiXmlParser::openRootElementOfFile( " << fileName
        << " ) could not open the file!";
        std::cout << std::endl;
      }
      return false;
    }
    textStream = &fileParsingStream;
    streamIsGood = skipPrologAndOpenRootElement();
    if( !streamIsGood
        &&
        isVerbose )
    {
      std::cout
      << std::endl
      << "BOL::error! AsciiXmlParser::openRootElementOfFile( " << fileName
      << " ) could not find a root element!";
      std::cout << std::endl;
    }
    return streamIsGood;
  }

  bool
  AsciiXmlParser::closeMarkup( size_t const startPosition )
  /* this records characters from textStream by appending them to markupString,
   * up to the 1st instance of markupCloser that is not enclosed in quotes,
   * but only looking for quote characters from startPosition onwards. an
   * exception is made if the markup was a comment: if markupString begins with
   * commentDelimiter.first, all characters up to the next found
   * ( commentDelimiter.first + markupCloser ) are discarded & markupString is
   * emptied, then true is returned.
   */
  {
    streamIsGood = ( recordTo( markupString,
                               markupCloser )
                     &&
                     ignoreDelimited( commentDelimiter ) );
    if( !streamIsGood )
    {
      if( isVerbose )
      {
        std::cout
        << std::endl
        << "BOL::error! AsciiXmlParser::closeMarkup( " << startPosition
        << " ) could not find the end of a markup!";
        std::cout << std::endl;
      }
      return false;
    }
    parseStart = markupString.find_first_of( allowedXmlQuoteChars,
                                             startPosition );
    while( std::string::npos != parseStart )
    {
      currentQuoteChar = markupString[ parseStart ];
      parseEnd = markupString.find( currentQuoteChar,
                                    ( parseStart + 1 ) );
      while( std::string::npos == parseEnd )
      {
        markupString.append( 1,
                             markupCloser );
        previousLength = markupString.size();
        streamIsGood = recordTo( markupString,
                                 markupCloser );
        if( !streamIsGood )
        {
          if( isVerbose )
          {
            std::cout
            << std::endl
            << "BOL::error! AsciiXmlParser::closeMarkup( " << startPosition
            << ") could not find the end of a markup!";
            std::cout << std::endl;
          }
          return false;
        }
        parseEnd = markupString.find( currentQuoteChar,
                                      previousLength );
      }
      parseStart = markupString.find_first_of( allowedXmlQuoteChars,
                                               ( parseEnd + 1 ) );
    }
    return true;
  }

  bool
  AsciiXmlParser::ignoreDelimited(
               std::pair< std::string, std::string > const& delimitingStrings )
  /* this checks to see if markupString begins with delimitingStrings.first,
   * & if so, ensures that markupString ends with delimitingStrings.second,
   * appending to markupString if necessary, then empties it. true is then
   * returned, unless the end of the text was reached before this could
   * happen.
   */
  {
    if( !(compareMarkupStart( delimitingStrings.first )) )
    {
      return true;
    }
    while( !(compareMarkupEnd( delimitingStrings.second )) )
    {
      streamIsGood = recordTo( markupString,
                               markupCloser );
      if( !streamIsGood )
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "BOL::error! AsciiXmlParser::ignoreDelimited( \""
          << delimitingStrings.first << ", " << delimitingStrings.second
          << "\" ) could not find the ending delimiter!";
          std::cout << std::endl;
        }
        return false;
      }
    }
    markupString.assign( "" );
    return true;
  }

  bool
  AsciiXmlParser::recordDelimited(
               std::pair< std::string, std::string > const& delimitingStrings )
  /* this checks to see if markupString begins with delimitingStrings.first,
   * & if so, ensures that markupString ends with delimitingStrings.second,
   * appending to markupString if necessary, then appends
   * markupOpener + markupString + markupCloser to fullElementContentAsFound,
   * then empties markupString. if markupString does not begin with
   * delimitingStrings.first, no change is made to either markupString or
   * fullElementContentAsFound. true is then returned, unless the end of the
   * text was reached before this could happen.
   */
  {
    if( !(compareMarkupStart( delimitingStrings.first )) )
    {
      return true;
    }
    while( !(compareMarkupEnd( delimitingStrings.second )) )
    {
      markupString.append( 1,
                           markupCloser );
      streamIsGood = recordTo( markupString,
                               markupCloser );
      if( !streamIsGood )
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "BOL::error! AsciiXmlParser::ignoreDelimited( \""
          << delimitingStrings.first << ", " << delimitingStrings.second
          << "\" ) could not find the ending delimiter!";
          std::cout << std::endl;
        }
        return false;
      }
    }
    fullElementContentAsFound.append( 1,
                                      markupOpener );
    fullElementContentAsFound.append( markupString );
    fullElementContentAsFound.append( 1,
                                      markupCloser );
    markupString.assign( "" );
    return true;
  }

  bool
  AsciiXmlParser::eraseQuotedStringsInMarkup( size_t const startPosition )
  /* this erases any quoted text in markupString starting from startPosition,
   * extending markupString from textStream to the next unquoted
   * markupCloser.
   */
  {
    // at this point, every opening quote in markupString from startPosition is
    // matched by a closing quote (properly nested).
    parseStart = markupString.find_first_of( allowedXmlQuoteChars,
                                             startPosition );
    while( std::string::npos != parseStart )
    {
      currentQuoteChar = markupString[ parseStart ];
      parseEnd = markupString.find( currentQuoteChar,
                                    ( parseStart + 1 ) );
      markupString.erase( parseStart,
                          ( parseEnd - parseStart + 1 ) );
      // we have to erase the closing quote character as well.
      parseStart = markupString.find_first_of( allowedXmlQuoteChars,
                                               parseStart );
    }
    // at this point, all quoted strings have been removed from markupString,
    // & now it ends at the 1st unquoted '>' in the text.
    return true;
  }

  bool
  AsciiXmlParser::parseAttributes()
  /* this parses any attributes in markupString, assuming that parseEnd is at
   * the end of the tag's name. false is returned if a malformed attribute is
   * found.
   */
  {
    elementAttributeMap.clear();
    parseStart = markupString.find_first_not_of( allowedXmlWhitespaceChars,
                                                 parseEnd );
    while( std::string::npos != parseStart )
    {
      if( ( ( markupString.size() - 1 ) == parseStart )
          &&
          ( tagCloser == markupString[ parseStart ] ) )
      {
        // if there is only whitespace left in the markup or the indicator of
        // an empty element, the parsing is done:
        return true;
      }
      // at this point, we have found a new attribute.
      parseEnd = markupString.find( '=',
                                    parseStart );
      if( std::string::npos == parseEnd )
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "BOL::error! AsciiXmlParser::parseOpeningTag() found an attribute"
          << " (\"" << markupString.substr( parseStart )
          << "\") without a value!";
          std::cout << std::endl;
        }
        return false;
      }
      currentAttribute.first.assign( markupString.substr( parseStart,
                                                 ( parseEnd - parseStart ) ) );
      currentQuoteChar = markupString[ ( ++parseEnd ) ];
      // this really should be ' or " for valid XML, but we won't bother
      // checking...
      parseStart = ( ++parseEnd );
      // the attribute's value begins after the '=' & the quote character (by
      // this point, parseEnd has been incremented twice before parseStart gets
      // set).
      parseEnd = markupString.find( currentQuoteChar,
                                    parseStart );
      if( std::string::npos == parseEnd )
        // this shouldn't ever happen, because such cases should have already
        // been caught by closeMarkup().
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "BOL::error! AsciiXmlParser::parseOpeningTag() found an attribute"
          << " (\"" << markupString.substr( parseStart )
          << "\") without a well-formed value (no closing quote mark)!";
          std::cout << std::endl;
        }
        return false;
      }
      currentAttribute.second.assign( markupString.substr( parseStart,
                                                 ( parseEnd - parseStart ) ) );
      elementAttributeMap.insert( currentAttribute );
      parseStart = markupString.find_first_not_of( allowedXmlWhitespaceChars,
                                                   ( ++parseEnd ) );
      // parseEnd has to be incremented so that parseStart doesn't just sit on
      // the closing quote.
    }
    return true;
  }

  bool
  AsciiXmlParser::recordToEndOfElement()
  /* this stores the characters between the opening tag & the corresponding
   * closing tag in fullElementContentAsFound, returning false if the end of
   * the text was reached before finding the closing tag.
   */
  {
    fullOpeningTagAsFound.assign( "" );
    recordTagTo( fullOpeningTagAsFound );
    fullElementContentAsFound.assign( "" );
    // we note which line the opening tag is on:
    elementLineRange.first = ( readNewlines + 1 );
    if( tagCloser == markupString[ markupString.size() - 1 ] )
    {
      // empty tags need no recording:
      elementLineRange.second = ( readNewlines + 1 );
      return true;
    }
    closingTag.assign( 1,
                       tagCloser );
    closingTag.append( elementName );
    closingTagsToFind = 1;
    // there could be nested elements of the same name (not necessarily nested
    // directly).
    while( 0 < closingTagsToFind )
    {
      streamIsGood = ( recordToNextTag()
                       &&
                       parseTagName( currentTagName ) );
      if( !streamIsGood )
      {
        return false;
      }
      if( 0 == currentTagName.compare( elementName ) )
      {
        // if a nested child element of the same name is found, the child
        // element is recorded too.
        ++closingTagsToFind;
      }
      else if( 0 == currentTagName.compare( closingTag ) )
      {
        --closingTagsToFind;
      }
      if( 0 < closingTagsToFind )
      {
        // if unless this is the closing tag of element, the markup must be
        // recorded, & here, because textStream has already gone past it.
        recordTagTo( fullElementContentAsFound );
      }
      else
      {
        // if this is the closing tag, we note which line it is on:
        elementLineRange.second = ( readNewlines + 1 );
      }
    }
    return true;
  }

  bool
  AsciiXmlParser::skipPrologAndOpenRootElement()
  /* this skips any XML prolog & then opens the root element. the prolog is an
   * optional XML declaration & an optional document type declaration, with any
   * number of comments & processing instructions (though not before the XML
   * declaration if there is one).
   */
  {
    resetContent();
    // the XML declaration markup is skipped as it is a special case of
    // processing instruction markup.
    while( markupString.empty() )
    {
      streamIsGood = ( discardToNextMarkup()
                       &&
                       ignoreDelimited( commentDelimiter )
                       &&
                       ignoreDelimited( piDelimiter ) );
      if( !streamIsGood )
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "BOL::error! AsciiXmlParser::skipPrologAndOpenRootElement() could"
          << " not find any valid root element tag markup!";
          std::cout << std::endl;
        }
        return false;
      }
    }
    // at this point, either the current markup is the root element's opening
    // tag, or is the document declaration if there is one, either only up to
    // the 1st '>' regardless of if it actually is the end of the markup or
    // not (e.g. if it is within quotes).
    if( ( doctypeDelimiter.first.size() < ( markupString.size() + 2 ) )
        &&
        compareMarkupStart( doctypeDelimiter.first ) )
      // if the markup is a document type declaration...
    {
      // since we will discard the document type declaration, it doesn't matter
      // if we mangle markupString in doing so, since it will be over-written
      // the opening tag of the root element anyway.
      streamIsGood = eraseQuotedStringsInMarkup( 0 );
      if( !streamIsGood )
      {
        if( isVerbose )
        {
          std::cout
          << std::endl
          << "BOL::error! AsciiXmlParser::skipPrologAndOpenRootElement() could"
          << " not find any valid root element tag markup!";
          std::cout << std::endl;
        }
        return false;
      }
      size_t unclosedSubmarkupOpener( markupString.find( markupOpener ) );
      while( std::string::npos != unclosedSubmarkupOpener )
      {
        previousLength = markupString.size();
        streamIsGood = ( recordTo( markupString,
                                   markupCloser )
                         &&
                         eraseQuotedStringsInMarkup( previousLength ) );
        if( !streamIsGood )
        {
          if( isVerbose )
          {
            std::cout
            << std::endl
            << "BOL::error! AsciiXmlParser::skipPrologAndOpenRootElement()"
            << " could not find any valid root element tag markup!";
            std::cout << std::endl;
          }
          return false;
        }
        unclosedSubmarkupOpener = markupString.find( markupOpener,
                                                     unclosedSubmarkupOpener );
      }
      markupString.assign( "" );
      // at this point, we have read in an unquoted '>' for every '<'. now we
      // can forget the mangled document type declaration & we only have to
      // discard any more comments or processing instructions before the
      // opening tag of the root element.
      while( markupString.empty() )
      {
        streamIsGood = ( discardToNextMarkup()
                         &&
                         ignoreDelimited( commentDelimiter )
                         &&
                         ignoreDelimited( piDelimiter ) );
        if( !streamIsGood )
        {
          if( isVerbose )
          {
            std::cout
            << std::endl
            << "BOL::error! AsciiXmlParser::skipPrologAndOpenRootElement()"
            << " could not find any valid root element tag markup!";
            std::cout << std::endl;
          }
          return false;
        }
      }
    }
    // at this point, the current markup is the root element's opening tag,
    // though only up to the 1st '>' regardless of if it actually is the end of
    // the markup or not (e.g. if it is within quotes).
    rootLineRange.first = ( readNewlines + 1 );
    streamIsGood = ( parseTagName( elementName )
                     &&
                     parseAttributes() );
    if( !streamIsGood )
    {
      if( isVerbose )
      {
        std::cout
        << std::endl
        << "BOL::error! AsciiXmlParser::skipPrologAndOpenRootElement() could"
        << " not find any valid root element tag markup!";
        std::cout << std::endl;
      }
      return false;
    }
    rootAttributeMap = elementAttributeMap;
    rootTag.assign( elementName );
    return true;
  }

}
