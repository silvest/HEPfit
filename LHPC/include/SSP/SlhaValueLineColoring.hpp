/*
 * SlhaValueLineColoring.hpp
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

#ifndef SLHAVALUELINECOLORING_HPP_
#define SLHAVALUELINECOLORING_HPP_

#include <iostream>
#include <string>
#include <fstream>

namespace LHPC
{
  namespace SLHA
  {

    class SlhaValueLineColoring
    {
    public:
      // This adds the line segment from the given co-ordinates to the gnuplot
      // input files, incrementing gnuplotLineIndex such that the data in
      // gnuplotDataFile with "index" = gnuplotLineIndex (before increment) is
      // the line, while it has line style ls gnuplotLineIndex (after
      // increment) - e.g. index 2 data has line style 3.
      // It does not plot the line if lineColor is "transparent".
      static void AddLine( int& gnuplotLineIndex,
                           std::ofstream& gnuplotDataFile,
                           std::ofstream& gnuplotCommandFile,
                           double const leftEndHorizontalPosition,
                           double const leftEndVerticalPosition,
                           double const rightEndHorizontalPosition,
                           double const rightEndVerticalPosition,
                           std::string const& lineColor );

      SlhaValueLineColoring();
      virtual ~SlhaValueLineColoring();

      // This should write all the line segments required with the correct
      // line styles. SlhaValueLineColoring::AddLine should be used to keep
      // the indexing with gnuplotLineIndex consistent.
      virtual void writeLineData( int& gnuplotLineIndex,
                                  std::ofstream& gnuplotDataFile,
                                  std::ofstream& gnuplotCommandFile,
                                  double const leftEndHorizontalPosition,
                                  double const lineVerticalPositon,
                           double const rightEndHorizontalPosition ) const = 0;
    };






    // This adds the line segment from the given co-ordinates to the gnuplot
    // input files, incrementing gnuplotLineIndex such that the data in
    // gnuplotDataFile with "index" = gnuplotLineIndex (before increment) is
    // the line, while it has line style ls gnuplotLineIndex (after
    // increment) - e.g. index 2 data has line style 3.
    inline void SlhaValueLineColoring::AddLine( int& gnuplotLineIndex,
                                                std::ofstream& gnuplotDataFile,
                                             std::ofstream& gnuplotCommandFile,
                                        double const leftEndHorizontalPosition,
                                          double const leftEndVerticalPosition,
                                       double const rightEndHorizontalPosition,
                                         double const rightEndVerticalPosition,
                                                std::string const& lineColor )
    {
      if( lineColor.compare( "transparent" ) != 0 )
      {
        gnuplotDataFile
        << leftEndHorizontalPosition << " " << leftEndVerticalPosition
        << std::endl
        << rightEndHorizontalPosition << " " << rightEndVerticalPosition
        << std::endl;
        // Double blank lines are treated by gnuplot as meaning that the data
        // above has a separate "index":
        gnuplotDataFile << std::endl << std::endl;
        gnuplotCommandFile
        << "set style line " << (++gnuplotLineIndex) << " lt rgb \""
        << lineColor << "\" lw 3"
        << std::endl;
      }
    }

  } /* namespace SLHA */
} /* namespace LHPC */

#endif /* SLHAVALUELINECOLORING_HPP_ */
