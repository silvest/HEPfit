/*
 * TwoDimensionalDataPlotter.cpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "TwoDimensionalDataPlotter.hpp"

namespace BOL
{
  TwoDimensionalDataPlotter::TwoDimensionalDataPlotter(
                                    std::string const& pathToGnuplotExecutable,
                                            std::string const& plotFileName ) :
      plotFileName( plotFileName ),
      gnuplotDataFileName( "gnuplot_data.txt" ),
      gnuplotCommandFileName( "gnuplot_input.txt" ),
      patienceTicks( 10000 ),
      dataPoints( 0 ),
      outputStream(),
      gnuplotExecutor( pathToGnuplotExecutable,
                       patienceTicks )
  {
    gnuplotExecutor.setArguments( gnuplotCommandFileName );
  }

  TwoDimensionalDataPlotter::~TwoDimensionalDataPlotter()
  {
    if( outputStream.is_open() )
    {
      outputStream.close();
    }
  }


  void
  TwoDimensionalDataPlotter::plotData( std::string const legendString,
                                       std::string const xAxisLabel,
                                       std::string const yAxisLabel )
  {
    outputStream.open( gnuplotDataFileName.c_str() );
    if( outputStream.is_open() )
      // if the file was successfully opened...
    {
      for( int dataIndex( 0 );
           dataPoints.getSize() > dataIndex;
           ++dataIndex )
      {
        outputStream
        << dataPoints[ dataIndex ].first << " "
        << dataPoints[ dataIndex ].second << std::endl;
      }
      outputStream.close();

      outputStream.open( gnuplotCommandFileName.c_str() );
      if( outputStream.is_open() )
        // if the file was successfully opened...
      {
        outputStream
        << "set term postscript eps enhanced color solid" << std::endl
        << "set output \"" << plotFileName << "\"" << std::endl
        << "set style line 1 lt rgb \"red\" lw 3" << std::endl
        << "set ylabel \"" << yAxisLabel << "\"" << std::endl
        << "set xlabel \"" << xAxisLabel << "\"" << std::endl
        << "plot '" << gnuplotDataFileName << "' index 0 title \""
        << legendString << "\" with lines ls 1" << std::endl;
        outputStream.close();

        gnuplotExecutor.forkAndExecvAndWait();
      }
      else
      {
        std::cout
        << std::endl
        << "BOL::error! TwoDimensionalDataPlotter::plotData() could not open"
        << gnuplotCommandFileName;
        std::cout << std::endl;
      }
    }
    else
    {
      std::cout
      << std::endl
      << "BOL::error! TwoDimensionalDataPlotter::plotData() could not open"
      << gnuplotDataFileName;
      std::cout << std::endl;
    }
  }

  void
  TwoDimensionalDataPlotter::plotData(
                std::vector< DoublePairVectorWithString > const& dataAndColors,
                                       std::string const xAxisLabel,
                                       std::string const yAxisLabel )
  {
    std::string dataFileName;
    size_t numberOfFileNumberDigits( StringParser::numberOfDigitsInInt(
                                                      dataAndColors.size() ) );

    std::vector< std::pair< std::string, std::string > > colorsWithFilenames;
    for( size_t whichEntry( 0 );
         dataAndColors.size() > whichEntry;
         ++whichEntry )
    {
      dataFileName.assign( gnuplotDataFileName );
      dataFileName.append( StringParser::intToString( whichEntry,
                                                      numberOfFileNumberDigits,
                                                      "" ) );
      colorsWithFilenames.push_back( std::pair< std::string, std::string >(
                                            dataAndColors[ whichEntry ].second,
                                                              dataFileName ) );

      outputStream.open( dataFileName.c_str() );
      if( outputStream.is_open() )
        // if the file was successfully opened...
      {
        for( unsigned int dataIndex( 0 );
             dataAndColors[ whichEntry ].first.size() > dataIndex;
             ++dataIndex )
        {
          outputStream
          << dataAndColors[ whichEntry ].first[ dataIndex ].first << " "
          << dataAndColors[ whichEntry ].first[ dataIndex ].second
          << std::endl;
        }
        outputStream.close();
      }
      else
      {
        std::cout
        << std::endl
        << "BOL::error! TwoDimensionalDataPlotter::plotData() could not open"
        << dataFileName;
        std::cout << std::endl;
      }
    }

    outputStream.open( gnuplotCommandFileName.c_str() );
    if( outputStream.is_open() )
      // if the file was successfully opened...
    {
      outputStream
      << "set term postscript eps enhanced color solid" << std::endl
      << "set ylabel \"" << yAxisLabel << "\"" << std::endl
      << "set xlabel \"" << xAxisLabel << "\"" << std::endl
      << "set output \"" << plotFileName << "\"" << std::endl;
      for( size_t whichEntry( 0 );
           colorsWithFilenames.size() > whichEntry;
           ++whichEntry )
      {
        outputStream
        << "set style line " << ( whichEntry + 1 ) << " lt rgb \""
        << colorsWithFilenames[ whichEntry ].first << "\" lw 3"
        << std::endl;
      }
      outputStream
      << "plot ";
      for( size_t whichEntry( 0 );
           colorsWithFilenames.size() > whichEntry;
           ++whichEntry )
      {
        if( 0 < whichEntry )
        {
          outputStream << ", ";
        }
        outputStream
        << "'" << colorsWithFilenames[ whichEntry ].second
        << "' index 0 notitle with lines ls " << ( whichEntry + 1 );
      }
      outputStream.close();

      gnuplotExecutor.forkAndExecvAndWait();
    }
    else
    {
      std::cout
      << std::endl
      << "BOL::error! TwoDimensionalDataPlotter::plotData() could not open"
      << gnuplotCommandFileName;
      std::cout << std::endl;
    }
  }

  void
  TwoDimensionalDataPlotter::plotData(
                                  PlotDataVector const& dataAndColorsAndLabels,
                                       std::string const xAxisLabel,
                                       std::string const yAxisLabel )
  {
    std::string dataFileName;
    size_t numberOfFileNumberDigits( StringParser::numberOfDigitsInInt(
                                             dataAndColorsAndLabels.size() ) );
    size_t numberOfFiles( dataAndColorsAndLabels.size() );

    std::vector< std::string > dataFileNames;
    for( size_t whichEntry( 0 );
         whichEntry < numberOfFiles;
         ++whichEntry )
    {
      dataFileName.assign( gnuplotDataFileName );
      dataFileName.append( StringParser::intToString( whichEntry,
                                                      numberOfFileNumberDigits,
                                                      "" ) );
      dataFileNames.push_back( dataFileName );
      outputStream.open( dataFileName.c_str() );
      if( outputStream.is_open() )
        // if the file was successfully opened...
      {
        for( size_t dataIndex( 0 );
             dataAndColorsAndLabels[ whichEntry ].first.size() > dataIndex;
             ++dataIndex )
        {
          outputStream
          << dataAndColorsAndLabels[ whichEntry ].first[ dataIndex ].first
          << " "
          << dataAndColorsAndLabels[ whichEntry ].first[ dataIndex ].second
          << std::endl;
        }
        outputStream.close();
      }
      else
      {
        std::cout
        << std::endl
        << "BOL::error! TwoDimensionalDataPlotter::plotData() could not open"
        << dataFileName;
        std::cout << std::endl;
      }
    }

    outputStream.open( gnuplotCommandFileName.c_str() );
    if( outputStream.is_open() )
      // if the file was successfully opened...
    {
      outputStream
      << "set term postscript eps enhanced color solid" << std::endl
      << "set ylabel \"" << yAxisLabel << "\"" << std::endl
      << "set xlabel \"" << xAxisLabel << "\"" << std::endl
      << "set output \"" << plotFileName << "\"" << std::endl;
      for( size_t whichEntry( 0 );
           whichEntry < numberOfFiles;
           ++whichEntry )
      {
        outputStream
        << "set style line " << ( whichEntry + 1 ) << " lt rgb \""
        << dataAndColorsAndLabels[ whichEntry ].second.first << "\" lw 3"
        << std::endl;
      }
      outputStream
      << "plot ";
      for( size_t whichEntry( 0 );
           whichEntry < numberOfFiles;
           ++whichEntry )
      {
        if( 0 < whichEntry )
        {
          outputStream << ", ";
        }
        outputStream
        << "'" << dataFileNames[ whichEntry ]
        << "' index 0 title \""
        << dataAndColorsAndLabels[ whichEntry ].second.second
        << "\" with lines ls " << ( whichEntry + 1 );
      }
      outputStream.close();

      gnuplotExecutor.forkAndExecvAndWait();
    }
    else
    {
      std::cout
      << std::endl
      << "BOL::error! TwoDimensionalDataPlotter::plotData() could not open"
      << gnuplotCommandFileName;
      std::cout << std::endl;
    }
  }
}
