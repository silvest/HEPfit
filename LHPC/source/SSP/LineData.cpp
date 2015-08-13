/*
 * LineData.cpp
 *
 *  Created on: Feb 26, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "SSP/LineData.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace SpectrumPlotting
    {
      BOL::StringParser const LineData::overlargeMassPrinter( 1,
                                                              '0',
                                                              3,
                                                              1,
                                                              "",
                                                              "-",
                                                              "",
                                                              "-",
                                                              "e" );

      LineData::LineData() :
          columnIndex( 0 ),
          whichJustification( centerJustified ),
          massValue( 0.0 ),
          labelPosition( -1.0 ),
          labelString( "error" ),
          colorString( "black" ),
          remainderString( "" )
      {
        // just an initialization list (with default unset values).
      }

      LineData::LineData( LineData const& copySource ) :
          columnIndex( copySource.columnIndex ),
          whichJustification( copySource.whichJustification ),
          massValue( copySource.massValue ),
          labelPosition( copySource.labelPosition ),
          labelString( copySource.labelString ),
          colorString( copySource.colorString ),
          remainderString( "" )
      {
        // just an initialization list.
      }

      LineData::~LineData()
      {
        // does nothing.
      }


      void
      LineData::setValues( std::string const& dataString,
                           double const massValue )
      {
        if( 0.0 > massValue )
        {
          this->massValue = -massValue;
        }
        else
        {
          this->massValue = massValue;
        }
        labelPosition = this->massValue;
        columnIndex
        = BOL::StringParser::stringToInt( BOL::StringParser::firstWordOf(
                                                                    dataString,
                                                              &remainderString,
                                        BOL::StringParser::whitespaceChars ) );
        if( 0 == ( columnIndex % 2 ) )
        {
          whichJustification = leftJustified;
          /* columns 2, 4, ... are the right-hand columns of their pairs & get
           * their labels flushed up against their right side, thus the labels'
           * left sides.
           */
        }
        else
        {
          whichJustification = rightJustified;
          /* columns 2, 4, ... are the left-hand columns of their pairs & get
           * their labels flushed up against their left side, thus the labels'
           * right sides.
           */
        }
        colorString.assign( BOL::StringParser::firstWordOf( remainderString,
                                                            &labelString,
                                        BOL::StringParser::whitespaceChars ) );
        labelString.assign( BOL::StringParser::trimFromFrontAndBack(
                                                                   labelString,
                                        BOL::StringParser::whitespaceChars ) );
      }

    }

  }

}
