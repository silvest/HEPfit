/*
 * SlhaValueSingleColor.hpp
 *
 *  Created on: Mar 24, 2015
 *      Author: bol
 */

#ifndef SLHAVALUESINGLECOLOR_HPP_
#define SLHAVALUESINGLECOLOR_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include "SlhaValueLineColoring.hpp"

namespace LHPC
{
  namespace SLHA
  {

    class SlhaValueSingleColor : public SlhaValueLineColoring
    {
    public:
      SlhaValueSingleColor( std::string const& colorDefinitionXml );
      virtual ~SlhaValueSingleColor();

      // This writes a single line segment in a single color.
      virtual void writeLineData( int& gnuplotLineIndex,
                                  std::ofstream& gnuplotDataFile,
                                  std::ofstream& gnuplotCommandFile,
                                  double const leftEndHorizontalPosition,
                                  double const lineVerticalPositon,
                                double const rightEndHorizontalPosition ) const
      { AddLine( gnuplotLineIndex,
                 gnuplotDataFile,
                 gnuplotCommandFile,
                 leftEndHorizontalPosition,
                 lineVerticalPositon,
                 rightEndHorizontalPosition,
                 lineVerticalPositon,
                 lineColor ); }

    protected:
      std::string const lineColor;
    };

  } /* namespace SLHA */
} /* namespace LHPC */

#endif /* SLHAVALUESINGLECOLOR_HPP_ */
