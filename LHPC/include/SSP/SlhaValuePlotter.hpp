/*
 * SlhaValuePlotter.hpp
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

#ifndef SLHAVALUEPLOTTER_HPP_
#define SLHAVALUEPLOTTER_HPP_

#include <cstdlib>
#include <exception>
#include <fstream>
#include <string>
#include <vector>
#include <list>
#include <map>
#include "BOLlib/include/BOLlib.hpp"
#include "../SLHA/BlockTypes.hpp"
#include "../SSP/LineData.hpp"
#include "SlhaValuePlotLine.hpp"
#include "SlhaValueLineColoring.hpp"
#include "SlhaValueSingleColor.hpp"
#include "SlhaValueColoredSegments.hpp"

namespace LHPC
{
  namespace SLHA
  {
    class SlhaValuePlotter
    {
    public:
      SlhaValuePlotter(std::string const& xmlInputFilename);
      ~SlhaValuePlotter();

      // This parses the child elements from the XML elements parsed from the
      // file given to the constructor, and then plots the values from the SLHA
      // file given by the XML file.
      void plotValues( std::string slhaFilename,
                       std::string plotFilename,
                       bool const shouldCleanUp );


    protected:
      static std::string const gnuplotDataFileName;
      static std::string const gnuplotCommandFileName;
      static std::string const gnuplotTexBaseName;
      static std::string const fullLatexBaseName;
      static double const flatBitWidth;
      static double const automaticScaleFactor;
      static double const labelSeparationShuffleFactor;
      static int const maximumLabelFloatingShuffles;

      std::string generalControlsXml;
      std::string columnDefinitionsXml;

      double marginWidth;
      double joinerWidth;
      double columnPairOffset;
      bool useAbsoluteValues;
      std::string verticalAxisLabel;
      double verticalAxisScaling;
      double verticalAxisUpperRange;
      double labelHeightFractionOfRange;
      std::string horizontalAxisLabel;
      double labelWidth;

      std::string gnuplotCommand;
      std::string latexCommand;
      std::string dvipsCommand;
      std::string ps2epsCommand;
      std::string deleteCommand;
      std::string moveCommand;

      std::vector< std::list< SlhaValuePlotLine > > columnSet;
      double labelHeight;
      double labelShuffleDistance;
      double fullColumnWidth;

      std::vector< SlhaValueLineColoring* > lineColoringPointers;

      // This parses the SLHA file and the plot output file names from
      // specificInputXml.
      void parseFilenames();

      // This parses the label settings and commands from generalControlsXml.
      void parseGeneralControls();

      // This parses the column definitions into columnSet.
      void createColumns( std::string slhaFilename );

      // This sorts every column in columnSet to arrange the lines in order,
      // and returns the largest vertical value found.
      double sortColumns();

      // This sets up the maximum vertical range, and takes care of labels for
      // lines which would be above the upper limit.
      void setUpVerticalRange( double const largestValueFound );

      // This tries to move all the labels up or down until they are all
      // separated by ( labelHeightFractionOfRange * verticalAxisUpperRange ).
      void floatLabels();

      // This prepares input files for gnuplot to plot the lines.
      void prepareTemporaryInputFiles();

      // This executes the commands to run gnuplot etc. on the prepared input
      // files.
      void executeCommands( std::string plotFilename,
                            bool const shouldCleanUp );

      // This adds a new column to columnSet based on the given definition.
      void addNewColumn( std::string columnDefinitionXml,
                         LHPC::SlhaSimplisticInterpreter& slhaParser );

      // This adds a new column line to the last column in columnSet based on
      // the given definition.
      void addColumnLine( std::string lineDefinitionXml,
                          LHPC::SlhaSimplisticInterpreter& slhaParser );

      // This creates an instance of a class derived from SlhaValueLineColoring
      // with "new" and returns a pointer to it, also storing it in
      // lineColoringPointers so that it can be correctly deleted later.
      inline SlhaValueLineColoring*
      getNewLineColoring( std::string const& drawType,
                          std::string const& xmlBody,
                          LHPC::SlhaSimplisticInterpreter& slhaParser );

      // This tries to move all the labels up or down until they are all
      // separated by ( labelHeightFractionOfRange * verticalAxisUpperRange ).
      void
      floatLabelsFromOneSide( std::vector< SlhaValuePlotLine* >& labelSet );
    };





    // This parses the child elements from the XML elements parsed from the
    // file given to the constructor, and then plots the values from the SLHA
    // file given by the XML file.
    inline void SlhaValuePlotter::plotValues( std::string slhaFilename,
                                              std::string plotFilename,
                                              bool const shouldCleanUp )
    {
      parseGeneralControls();
      createColumns( slhaFilename );
      double const largestValue( sortColumns() );
      setUpVerticalRange( largestValue );
      floatLabels();
      prepareTemporaryInputFiles();
      executeCommands( plotFilename,
                       shouldCleanUp );
    }

    // This parses the column definitions into columnSet.
    inline void SlhaValuePlotter::createColumns( std::string slhaFilename )
    {
      LHPC::SlhaSimplisticInterpreter slhaParser( slhaFilename );

      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.loadString( columnDefinitionsXml ) )
      {
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "ColumnDefinition" ) )
          {
            addNewColumn( xmlParser.getTrimmedCurrentElementContent(),
                          slhaParser );
          }
        }
      }
      else
      {
        throw std::runtime_error( "Could not parse <ColumnDefinitions>." );
      }
    }

    // This sorts every column in columnSet to arrange the lines in order,
    // and then notes the largest vertical value found.
    inline double SlhaValuePlotter::sortColumns()
    {
      double largestValue( 0.0 );
      for( std::vector< std::list< SlhaValuePlotLine > >::iterator
           whichColumn( columnSet.begin() );
           whichColumn < columnSet.end();
           ++whichColumn )
      {
        whichColumn->sort(&(SlhaValuePlotLine::lowToHigh));
        if( whichColumn->back().getValue() > largestValue )
        {
          largestValue =  whichColumn->back().getValue();
        }
      }
      return largestValue;
    }

    // This creates an instance of a class derived from SlhaValueLineColoring
    // with "new" and returns a pointer to it, also storing it in
    // lineColoringPointers so that it can be correctly deleted later.
    inline SlhaValueLineColoring*
    SlhaValuePlotter::getNewLineColoring( std::string const& drawType,
                                          std::string const& xmlBody,
                                  LHPC::SlhaSimplisticInterpreter& slhaParser )
    {
      if( drawType.compare( "SingleColor" ) == 0 )
      {
        lineColoringPointers.push_back( new SlhaValueSingleColor( xmlBody ) );
      }
      else if( drawType.compare( "ColoredSegments" ) == 0 )
      {
        lineColoringPointers.push_back( new SlhaValueColoredSegments( xmlBody,
                                                                  flatBitWidth,
                                                                slhaParser ) );
      }
      else
      {
        throw std::runtime_error( "Unknown \"DrawType\" in <LineColor> ("
                                  + drawType + ")." );
      }
      return lineColoringPointers.back();
    }

  } /* namespace SLHA */
} /* namespace LHPC */

#endif /* SLHAVALUEPLOTTER_HPP_ */
