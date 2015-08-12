/*
 * LineData.hpp
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

#ifndef LINEPLOTTINGDATA_HPP_
#define LINEPLOTTINGDATA_HPP_

#include <string>
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace SpectrumPlotting
    {
      // this class holds the co-ordinates of the line to plot to visualize the
      // mass of a particle from an SLHA file, along with the label.
      class LineData
      {
      public:
        enum JustificationStyle
        {
          leftJustified,
          centerJustified,
          rightJustified
        };

        static bool
        lowToHigh( LineData const& firstLineData,
                   LineData const& secondLineData );
        // this returns true if firstLineData has a valueDouble lower than or
        // equal to that of secondLineData.

        LineData();
        LineData( LineData const& copySource );
        ~LineData();

        void
        setValues( std::string const& dataString,
                   double const massValue );
        JustificationStyle
        getJustification() const;
        void
        setJustification( JustificationStyle const whichJustification );
        int
        getColumn() const;
        double
        getMass() const;
        double
        getLabelPosition() const;
        void
        setLabelPosition( double const labelPosition );
        std::string const&
        getLabelString() const;
        void
        relabelForOverlargeMass( double const scaleMaximum );
        std::string const&
        getColor() const;


      protected:
        static BOL::StringParser const overlargeMassPrinter;

        int columnIndex;
        JustificationStyle whichJustification;
        double massValue;
        double labelPosition;
        std::string labelString;
        std::string colorString;
        std::string remainderString;
      };





      inline bool
      LineData::lowToHigh( LineData const& firstLineData,
                           LineData const& secondLineData )
      // this returns true if firstLineData has a valueDouble lower than or
      // equal to that of secondLineData.
      {
        if( firstLineData.getMass() > secondLineData.getMass() )
        {
          return false;
        }
        else
        {
          return true;
        }
      }

      inline LineData::JustificationStyle
      LineData::getJustification() const
      {
        return whichJustification;
      }

      inline void
      LineData::setJustification(
                        LineData::JustificationStyle const whichJustification )
      {
        this->whichJustification = whichJustification;
      }

      inline int
      LineData::getColumn() const
      {
        return columnIndex;
      }

      inline double
      LineData::getMass() const
      {
        return massValue;
      }

      inline double
      LineData::getLabelPosition() const
      {
        return labelPosition;
      }

      inline void
      LineData::setLabelPosition( double const labelPosition )
      {
        this->labelPosition = labelPosition;
      }

      inline std::string const&
      LineData::getLabelString() const
      {
        return labelString;
      }

      inline void
      LineData::relabelForOverlargeMass( double const scaleMaximum )
      {
        whichJustification = centerJustified;
        labelString.append( "{\\footnotesize (" );
        labelString.append( overlargeMassPrinter.doubleToString( massValue ) );
        labelString.append( ")}" );
        labelPosition = ( 0.99 * scaleMaximum );
      }

      inline std::string const&
      LineData::getColor() const
      {
        return colorString;
      }

    }

  }

}

#endif /* LINEPLOTTINGDATA_HPP_ */
