/*
 * SlhaValuePlotLine.cpp
 *
 *  Created on: Mar 22, 2015
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2015 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "../../include/SSP/SlhaValuePlotLine.hpp"

namespace LHPC
{
  namespace SLHA
  {
    std::string SlhaValuePlotLine::joinerColor( "black" );
    BOL::StringParser const SlhaValuePlotLine::overlargeValuePrinter( 1,
                                                                      '0',
                                                                      3,
                                                                      1,
                                                                      "",
                                                                      "-",
                                                                      "",
                                                                      "-",
                                                                      "e" );

    SlhaValuePlotLine::SlhaValuePlotLine( double const slhaValue,
                                          double const columnCenterPosition,
                                          double const columnLineWidth,
                                          double const columnPairOffset,
                                          double const labelJoinerWidth,
                                          std::string const& labelString,
                                          bool const labelLeftOfColumn,
                                  SlhaValueLineColoring const* lineColoring ) :
        slhaValue( slhaValue ),
        leftEndHorizontalPosition( columnCenterPosition
                                   - (0.5 * ( columnLineWidth
                                              + columnPairOffset ) ) ),
        rightEndHorizontalPosition( leftEndHorizontalPosition
                                    + columnLineWidth),
        labelVerticalPosition( slhaValue ),
        labelHorizontalPosition( leftEndHorizontalPosition
                                 - labelJoinerWidth ),
        labelString( labelString ),
        labelLeftOfColumn( labelLeftOfColumn ),
        valueAbovePlotRange( false ),
        lineColoring( lineColoring )
    {
      if( !labelLeftOfColumn )
      {
        leftEndHorizontalPosition += columnPairOffset;
        rightEndHorizontalPosition += columnPairOffset;
        labelHorizontalPosition
        = ( rightEndHorizontalPosition + labelJoinerWidth );
      }
    }

    SlhaValuePlotLine::SlhaValuePlotLine(
                                        SlhaValuePlotLine const& copySource ) :
        slhaValue( copySource.slhaValue ),
        leftEndHorizontalPosition( copySource.leftEndHorizontalPosition ),
        rightEndHorizontalPosition( copySource.rightEndHorizontalPosition ),
        labelVerticalPosition( copySource.labelVerticalPosition ),
        labelHorizontalPosition( copySource.labelHorizontalPosition ),
        labelString( copySource.labelString ),
        labelLeftOfColumn( copySource.labelLeftOfColumn ),
        valueAbovePlotRange( copySource.valueAbovePlotRange ),
        lineColoring( copySource.lineColoring )
    {
      // This constructor is just an initialization list.
    }

    SlhaValuePlotLine::~SlhaValuePlotLine()
    {
      // The lineColoring pointer is not deleted by this object, as then copies
      // of this object (and there is a lot of copying and deleting of old
      // SlhaValuePlotLine objects in the course of std::list objects being
      // copied during std::vector re-allocations) could not use the pointer.
    }


    // This writes input to have gnuplot print the line label correctly.
    void SlhaValuePlotLine::writeLabel( std::ofstream& gnuplotCommandFile )
    {
      std::string gnuplotLabelString( "{" );
      // This case only should be center-justified.
      if( valueAbovePlotRange )
      {
        gnuplotLabelString.append( labelString );
      }
      else if( labelLeftOfColumn )
      {
        gnuplotLabelString.append( "{" );
        gnuplotLabelString.append( labelString );
        gnuplotLabelString.append( "}{\\phantom{" );
        gnuplotLabelString.append( labelString );
        gnuplotLabelString.append( "}}" );
      }
      else
      {
        gnuplotLabelString.append( "{\\phantom{" );
        gnuplotLabelString.append( labelString );
        gnuplotLabelString.append( "}{" );
        gnuplotLabelString.append( labelString );
        gnuplotLabelString.append( "}}" );
      }
      gnuplotLabelString.append( "}" );

      gnuplotCommandFile
      << "set label '" << gnuplotLabelString << "' at "
      << labelHorizontalPosition << ", " << labelVerticalPosition << std::endl;
    }

    // This writes the line as a segment or series of segments defined as
    // lines themselves in gnuplotDataFile, using gnuplotLineIndex to keep
    // handles on the lines and their styles.
    void SlhaValuePlotLine::writeLineData( int& gnuplotLineIndex,
                                           std::ofstream& gnuplotDataFile,
                                           std::ofstream& gnuplotCommandFile )
    {
      // The "joiner" from the label to the end of the line is a separate line.
      if( labelLeftOfColumn )
      {
        SlhaValueLineColoring::AddLine( gnuplotLineIndex,
                                        gnuplotDataFile,
                                        gnuplotCommandFile,
                                        labelHorizontalPosition,
                                        labelVerticalPosition,
                                        leftEndHorizontalPosition,
                                        slhaValue,
                                        joinerColor );
      }
      else
      {
        SlhaValueLineColoring::AddLine( gnuplotLineIndex,
                                        gnuplotDataFile,
                                        gnuplotCommandFile,
                                        rightEndHorizontalPosition,
                                        slhaValue,
                                        labelHorizontalPosition,
                                        labelVerticalPosition,
                                        joinerColor );
      }

      if( !valueAbovePlotRange )
      {
        lineColoring->writeLineData( gnuplotLineIndex,
                                     gnuplotDataFile,
                                     gnuplotCommandFile,
                                     leftEndHorizontalPosition,
                                     slhaValue,
                                     rightEndHorizontalPosition );
      }

    }

  } /* namespace SLHA */
} /* namespace LHPC */
