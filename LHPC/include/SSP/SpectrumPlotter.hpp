/*
 * SpectrumPlotter.hpp
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

#ifndef SPECTRUMPLOTTER_HPP_
#define SPECTRUMPLOTTER_HPP_

#include <cstdlib>
#include <fstream>
#include <string>
#include <list>
#include <map>
#include "BOLlib/include/BOLlib.hpp"
#include "../SLHA/BlockTypes.hpp"
#include "../SSP/LineData.hpp"

namespace LHPC
{
  namespace SLHA
  {
    // this class does the job of reading in data from a SLHA MASS or FMASS
    // block & plotting it according to a LHPCPLOTTING block.
    class SpectrumPlotter
    {
    public:
      typedef SparseSinglyIndexedBlock< std::string > StringBlock;
      typedef std::map< int, std::string > LineMap;
      typedef std::list< SpectrumPlotting::LineData > LineList;
      typedef SparseSinglyIndexedBlock< RunningConstant > FmassBlock;
      typedef SparseSinglyIndexedBlock< double > MassBlock;

      SpectrumPlotter( StringBlock const& plotControlBlock,
                       StringBlock const& linePlottingBlock,
                       FmassBlock const* const fmassPointer,
                       MassBlock const* const massPointer = NULL );
      virtual
      ~SpectrumPlotter();

      bool
      plotSpectrum( std::string const& plotFileName,
                    bool const shouldCleanUp = true );


    protected:
      static int const unitIndex;
      static int const scaleIndex;
      static int const labelYSizeIndex;
      static int const labelXSizeIndex;
      static int const gnuplotIndex;
      static int const latexIndex;
      static int const dvipsIndex;
      static int const ps2epsIndex;
      static int const rmIndex;
      static int const mvIndex;
      static std::string const gnuplotDataFileName;
      static std::string const gnuplotCommandFileName;
      static std::string const gnuplotTexBaseName;
      static std::string const fullLatexBaseName;
      static double const automaticScaleFactor;
      static double const labelSeparationShuffleFactor;
      static double const marginWidth;
      static double const joinerWidth;
      static double const flatBitWidth;
      static double const columnPairOffset;
      static int const maximumLabelFloatingShuffles;

      StringBlock const& plotControlBlock;
      StringBlock const& linePlottingBlock;
      FmassBlock const* const fmassPointer;
      MassBlock const* const massPointer;
      std::string unitString;
      double unitFactor;
      double scaleMaximum;
      double largestMass;
      BOL::VectorlikeArray< LineList > columnSet;
      LineList* columnPointer;
      LineMap const* plotLineMap;
      LineMap::const_iterator lineIterator;
      SpectrumPlotting::LineData lineAdder;
      std::list< SpectrumPlotting::LineData >::iterator lowerMassIterator;
      std::list< SpectrumPlotting::LineData >::iterator upperMassIterator;
      bool notYetFinishedShuffling;
      int remainingShuffles;
      int whichMassEigenstate;
      double massValue;
      double labelRoomWidth;
      double labelLatexWidth;
      double fullColumnWidth;
      double labelSeparation;
      double labelAverage;
      bool lastOperationSuccessful;
      int systemCallReturn;
      std::string fullLatexFilename;
      std::string gnuplotCommand;
      std::string latexCommand;
      std::string dvipsCommand;
      std::string ps2epsCommand;
      bool epsiInstead;
      std::string mainCleanupCommand;
      std::string moveCommand;
      bool leftColumnRatherThanRight;
      double leftLineXValue;
      double middleLineXValue;
      double rightLineXValue;
      int gnuplotLineIndex;
      std::string gnuplotLabelString;

      void
      loadCommands( std::string const& plotFileName );
      void
      loadLines();
      void
      sortMasses();
      // this sorts all the masses in the columns, & then sets the scale range.
      void
      floatLabels();
      // this tries to move all the labels up or down until they are all
      // separated by ( labelSeparation * scaleMaximum ).
      bool
      writeGnuplotFiles();
    };

  }

}

#endif /* SPECTRUMPLOTTER_HPP_ */
