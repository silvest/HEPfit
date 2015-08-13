/*
 * SlhaValueColoredSegments.hpp
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

#ifndef SLHAVALUECOLOREDSEGMENTS_HPP_
#define SLHAVALUECOLOREDSEGMENTS_HPP_

#include <string>
#include <map>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include "../SLHA/SlhaSimplisticInterpreter.hpp"
#include "../BOLlib/include/AsciiXmlParser.hpp"
#include "../BOLlib/include/StringParser.hpp"
#include "SlhaValueLineColoring.hpp"

namespace LHPC
{
  namespace SLHA
  {

    class SlhaValueColoredSegments : public SlhaValueLineColoring
    {
    public:
      SlhaValueColoredSegments( std::string const& colorDefinitionXml,
                                 double const coloredLineWidth,
                                 LHPC::SlhaSimplisticInterpreter& slhaParser );
      virtual ~SlhaValueColoredSegments();

      // This writes the sets of colors as segments with lengths given by their
      // fractional weighting of width of the line for the value given to the
      // constructor.
      virtual void writeLineData( int& gnuplotLineIndex,
                                  std::ofstream& gnuplotDataFile,
                                  std::ofstream& gnuplotCommandFile,
                                  double const leftEndHorizontalPosition,
                                  double const lineVerticalPositon,
                               double const rightEndHorizontalPosition ) const;


    protected:
      std::vector< std::pair< std::string, double > > colorsWithWidths;
      double cumulativeWeightTotal;

      // This adds a color with an unnormalized line segment length weight to
      // colorsWithWidths based on the XML in combination with the given SLHA
      // values, as well as updating the weight total. The widths of the
      // segments will be correctly normalized later in the constructor.
      void addSegment( std::string segmentXml,
                       LHPC::SlhaSimplisticInterpreter& slhaParser );
    };





    // This writes the sets of colors as segments with lengths given by their
    // fractional weighting of width of the line for the value given to the
    // constructor.
    inline void SlhaValueColoredSegments::writeLineData( int& gnuplotLineIndex,
                                                std::ofstream& gnuplotDataFile,
                                             std::ofstream& gnuplotCommandFile,
                                        double const leftEndHorizontalPosition,
                                              double const lineVerticalPositon,
                                double const rightEndHorizontalPosition ) const
    {
      double segmentLeftEnd( leftEndHorizontalPosition );
      for( std::vector< std::pair< std::string, double > >::const_iterator
           colorWithWidth( colorsWithWidths.begin() );
           colorWithWidth < colorsWithWidths.end();
           ++colorWithWidth )
      {
        double segmentRightEnd( segmentLeftEnd + colorWithWidth->second );
        AddLine( gnuplotLineIndex,
                 gnuplotDataFile,
                 gnuplotCommandFile,
                 segmentLeftEnd,
                 lineVerticalPositon,
                 segmentRightEnd,
                 lineVerticalPositon,
                 colorWithWidth->first );
        segmentLeftEnd = segmentRightEnd;
      }
    }

  } /* namespace SLHA */
} /* namespace LHPC */

#endif /* SLHAVALUECOLOREDSEGMENTS_HPP_ */
