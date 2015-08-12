/*
 * SlhaValueColoredSegments.cpp
 *
 *  Created on: Mar 24, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2015 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../../include/SSP/SlhaValueColoredSegments.hpp"

namespace LHPC
{
  namespace SLHA
  {

    SlhaValueColoredSegments::SlhaValueColoredSegments(
                                         std::string const& colorDefinitionXml,
                                                 double const coloredLineWidth,
                                LHPC::SlhaSimplisticInterpreter& slhaParser ) :
        SlhaValueLineColoring(),
        colorsWithWidths(),
        cumulativeWeightTotal( 0.0 )
    {
      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.loadString( colorDefinitionXml ) )
      {
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "LineSegment" ) )
          {
            addSegment( xmlParser.getTrimmedCurrentElementContent(),
                        slhaParser );
          }
        }

        // Now that all the weights have been read, their widths can be set.
        if( cumulativeWeightTotal == 0.0 )
        {
          double const equalWeight( coloredLineWidth
                          / static_cast< double >( colorsWithWidths.size() ) );
          for( std::vector< std::pair< std::string, double > >::iterator
               colorWithWidth( colorsWithWidths.begin() );
               colorWithWidth < colorsWithWidths.end();
               ++colorWithWidth )
          {
            colorWithWidth->second = equalWeight;
          }
        }
        else
        {
          double const
          normalizationFactor( coloredLineWidth / cumulativeWeightTotal );

          for( std::vector< std::pair< std::string, double > >::iterator
               colorWithWidth( colorsWithWidths.begin() );
               colorWithWidth < colorsWithWidths.end();
               ++colorWithWidth )
          {
            colorWithWidth->second *= normalizationFactor;
          }
        }
      }
      else
      {
        throw std::runtime_error( "Could not parse <LineColor>." );
      }
    }

    SlhaValueColoredSegments::~SlhaValueColoredSegments()
    {
      // This does nothing.
    }


    // This adds a color with an unnormalized line segment length weight to
    // colorsWithWidths based on the XML in combination with the given SLHA
    // values, as well as updating the weight total. The widths of the segments
    // will be correctly normalized later in the constructor.
    void SlhaValueColoredSegments::addSegment( std::string segmentXml,
                                 LHPC::SlhaSimplisticInterpreter& slhaParser )
    {
      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.loadString( segmentXml ) )
      {
        std::string segmentColor( "" );
        std::string segmentValue( "0.0" );
        double powerOfValueForWeight( 1.0 );
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "SegmentColor" ) )
          {
            segmentColor.assign( xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "ColorWeight" ) )
          {
            std::map< std::string, std::string > const&
            weightAttributes( xmlParser.getCurrentElementAttributes() );
            std::map< std::string, std::string >::const_iterator
            attributeFinder( weightAttributes.find("RaiseToPower") );
            if( attributeFinder != weightAttributes.end() )
            {
              powerOfValueForWeight
              = BOL::StringParser::stringToDouble( attributeFinder->second );
            }
            segmentValue.assign( xmlParser.getTrimmedCurrentElementContent() );
          }
        }
        double interpretedValue( 0.0 );

        // The evaluation of the conditional will set interpretedValue
        // correctly if it is just a number (i.e. if stringIsDouble returns
        // true).
        if( !(BOL::StringParser::stringIsDouble( segmentValue,
                                                 interpretedValue )) )
        {
          interpretedValue = slhaParser.getDouble( segmentValue );
        }

        // Now that the value has been correctly parsed, we can turn it into a
        // weight.
        double segmentWeight = pow( interpretedValue, powerOfValueForWeight );
        cumulativeWeightTotal += segmentWeight;
        colorsWithWidths.push_back( std::make_pair( segmentColor,
                                                    segmentWeight ) );
      }
      else
      {
        throw std::runtime_error( "Could not parse <LineSegment>." );
      }
    }

  } /* namespace SLHA */
} /* namespace LHPC */
