/*
 * SlhaValuePlotLine.hpp
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

#ifndef SLHAVALUEPLOTLINE_HPP_
#define SLHAVALUEPLOTLINE_HPP_

#include <string>
#include <map>
#include <fstream>
#include "../SLHA/SlhaSimplisticInterpreter.hpp"
#include "../BOLlib/include/AsciiXmlParser.hpp"
#include "../BOLlib/include/StringParser.hpp"
#include "SlhaValueLineColoring.hpp"
#include "SlhaValueSingleColor.hpp"
#include "SlhaValueColoredSegments.hpp"

namespace LHPC
{
  namespace SLHA
  {

    class SlhaValuePlotLine
    {
    public:
      static bool
      lowToHigh( SlhaValuePlotLine const& firstLine,
                 SlhaValuePlotLine const& secondLine )
      { return ( firstLine.getValue() <= secondLine.getValue() ); }

      SlhaValuePlotLine( double const slhaValue,
                         double const columnCenterPosition,
                         double const columnLineWidth,
                         double const columnPairOffset,
                         double const labelJoinerWidth,
                         std::string const& labelString,
                         bool const labelLeftOfColumn,
                         SlhaValueLineColoring const* lineColoring );
      SlhaValuePlotLine( SlhaValuePlotLine const& copySource );
      ~SlhaValuePlotLine();

      double getValue() const { return slhaValue; }
      void setValue( double const slhaValue ) { this->slhaValue = slhaValue; }

      double getLabelPosition() const { return labelVerticalPosition; }
      void setLabelPosition( double const labelPosition )
      { this->labelVerticalPosition = labelPosition; }

      bool hasLabelOnLeftOfColumn() const { return labelLeftOfColumn; }

      // This places the label just below the upper edge of the plot, and
      // appends the value in brackets to the label.
      void relabelForAboveRange( double const verticalAxisUpperRange );

      // This writes the line as a segment or series of segments defined as
      // lines themselves in gnuplotDataFile, using gnuplotLineIndex to keep
      // handles on the lines and their styles.
      void writeLineData( int& gnuplotLineIndex,
                          std::ofstream& gnuplotDataFile,
                          std::ofstream& gnuplotCommandFile );

      // This writes input to have gnuplot print the line label correctly.
      void writeLabel( std::ofstream& gnuplotCommandFile );


    protected:
      static std::string joinerColor;
      static BOL::StringParser const overlargeValuePrinter;

      double slhaValue;
      double leftEndHorizontalPosition;
      double rightEndHorizontalPosition;
      double labelVerticalPosition;
      double labelHorizontalPosition;
      std::string labelString;
      bool labelLeftOfColumn;
      bool valueAbovePlotRange;
      SlhaValueLineColoring const* lineColoring;
    };





    // This places the label just below the upper edge of the plot, and
    // appends the value in brackets to the label.
    inline void
    SlhaValuePlotLine::relabelForAboveRange( double const scaleMaximum )
    {
      labelLeftOfColumn = true;
      valueAbovePlotRange = true;
      labelString.append( "{\\footnotesize (" );
      labelString.append( overlargeValuePrinter.doubleToString( slhaValue ) );
      labelString.append( ")}" );
      labelVerticalPosition = ( 0.99 * scaleMaximum );
      leftEndHorizontalPosition
      = ( 0.5 * ( leftEndHorizontalPosition + rightEndHorizontalPosition ) );
      rightEndHorizontalPosition = leftEndHorizontalPosition;
    }

  } /* namespace SLHA */
} /* namespace LHPC */

#endif /* SLHAVALUEPLOTLINE_HPP_ */
