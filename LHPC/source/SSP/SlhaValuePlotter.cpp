/*
 * SlhaValuePlotter.cpp
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

#include "../../include/SSP/SlhaValuePlotter.hpp"

namespace LHPC
{
  namespace SLHA
  {
    std::string const
    SlhaValuePlotter::gnuplotDataFileName( "LHPC_SlhaPlotter_gnuplot.dat" );
    std::string const SlhaValuePlotter::gnuplotCommandFileName(
                                            "LHPC_SlhaPlotter_gnuplot.input" );
    std::string const
    SlhaValuePlotter::gnuplotTexBaseName( "LHPC_SlhaPlotter_gnuplot_TeX" );
    std::string const
    SlhaValuePlotter::fullLatexBaseName( "LHPC_SlhaPlotter_LaTeX" );
    double const SlhaValuePlotter::automaticScaleFactor( 1.1 );
    double const SlhaValuePlotter::labelSeparationShuffleFactor( 0.51 );
    double const SlhaValuePlotter::flatBitWidth( 100.0 );

    int const SlhaValuePlotter::maximumLabelFloatingShuffles( 100 );

    SlhaValuePlotter::SlhaValuePlotter(std::string const& xmlInputFilename) :
        generalControlsXml( "" ),
        columnDefinitionsXml( "" ),
        marginWidth( 10.0 ),
        joinerWidth( 15.0 ),
        columnPairOffset( 10.0 ),
        useAbsoluteValues( false ),
        verticalAxisLabel( "" ),
        verticalAxisScaling( -1.0 ),
        verticalAxisUpperRange( -1.0 ),
        labelHeightFractionOfRange( -1.0 ),
        horizontalAxisLabel( "" ),
        labelWidth( -1.0 ),
        gnuplotCommand( "" ),
        latexCommand( "" ),
        dvipsCommand( "" ),
        ps2epsCommand( "" ),
        deleteCommand( "" ),
        moveCommand( "" ),
        columnSet(),
        labelHeight( -1.0 ),
        labelShuffleDistance( -1.0 ),
        fullColumnWidth( -1.0 ),
        lineColoringPointers()
    {
      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.openRootElementOfFile( xmlInputFilename ) )
      {
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "GeneralControls" ) )
          {
            generalControlsXml.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "ColumnDefinitions" ) )
          {
            columnDefinitionsXml.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
        }
        xmlParser.closeFile();
      }
      else
      {
        throw std::runtime_error( "Could not open "
                                  + xmlInputFilename
                                  + " as valid XML file." );
      }
    }

    SlhaValuePlotter::~SlhaValuePlotter()
    {
      for( size_t deletionIndex( 0 );
           deletionIndex < lineColoringPointers.size();
           ++deletionIndex )
      {
        delete lineColoringPointers[ deletionIndex ];
      }
    }


    // This parses the label settings and commands from generalControlsXml.
    void SlhaValuePlotter::parseGeneralControls()
    {
      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.loadString( generalControlsXml ) )
      {
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "UseAbsoluteValues" ) )
          {
            std::string
            inputXml( xmlParser.getTrimmedCurrentElementContent() );
            BOL::StringParser::transformToLowercase( inputXml );
            useAbsoluteValues = ( inputXml.compare( "true" ) == 0 );
          }
          else if( xmlParser.currentElementNameMatches( "VerticalAxisLabel" ) )
          {
            verticalAxisLabel.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches(
                                                      "VerticalAxisScaling" ) )
          {
            verticalAxisScaling
            = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches(
                                                   "VerticalAxisUpperRange" ) )
          {
            verticalAxisUpperRange
            = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "LabelHeight" ) )
          {
            labelHeightFractionOfRange
            = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches(
                                                      "HorizontalAxisLabel" ) )
          {
            horizontalAxisLabel.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "LabelWidth" ) )
          {
            labelWidth
            = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "GnuplotExecutable" ) )
          {
            gnuplotCommand.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "LatexExecutable" ) )
          {
            latexCommand.assign( xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "DvipsExecutable" ) )
          {
            dvipsCommand.assign( xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "Ps2epsExecutable" ) )
          {
            ps2epsCommand.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "RemoveCommand" ) )
          {
            deleteCommand.assign(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "MoveCommand" ) )
          {
            moveCommand.assign( xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "MarginWidth" ) )
          {
            marginWidth = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "JoinerWidth" ) )
          {
            joinerWidth = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
          else if( xmlParser.currentElementNameMatches( "ColumnPairOffset" ) )
          {
            columnPairOffset = BOL::StringParser::stringToDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
          }
        }
        gnuplotCommand.append( " " );
        gnuplotCommand.append( gnuplotCommandFileName );
        // "gnuplot LHPC_SlhaPlotter_gnuplot.input"
        latexCommand.append( " " );
        latexCommand.append( fullLatexBaseName );
        latexCommand.append( ".tex" );
        // "latex LHPC_SlhaPlotter_LaTeX.tex"
        dvipsCommand.append( " " );
        dvipsCommand.append( fullLatexBaseName );
        dvipsCommand.append( ".dvi" );
        // "dvips LHPC_SlhaPlotter_LaTeX.dvi"
        if( 0 == ps2epsCommand.compare( ( ps2epsCommand.size() - 7 ),
                                        7,
                                        "ps2epsi" ) )
        {
          ps2epsCommand.append( " " );
          ps2epsCommand.append( fullLatexBaseName );
          ps2epsCommand.append( ".ps" );
          ps2epsCommand.append( " " );
          ps2epsCommand.append( fullLatexBaseName );
          ps2epsCommand.append( ".eps" );
        }
        else
        {
          ps2epsCommand.append( " -f " );
          ps2epsCommand.append( fullLatexBaseName );
          ps2epsCommand.append( ".ps" );
        }
        // "ps2eps - f LHPC_SlhaPlotter_LaTeX.ps" or
        // "ps2epsi LHPC_SlhaPlotter_LaTeX.ps LHPC_SlhaPlotter_LaTeX.eps"
        deleteCommand.append( " " );
        deleteCommand.append( gnuplotDataFileName );
        deleteCommand.append( " " );
        deleteCommand.append( gnuplotCommandFileName );
        deleteCommand.append( " " );
        deleteCommand.append( gnuplotTexBaseName );
        deleteCommand.append( ".eps " );
        deleteCommand.append( gnuplotTexBaseName );
        deleteCommand.append( ".tex " );
        deleteCommand.append( fullLatexBaseName );
        deleteCommand.append( ".aux " );
        deleteCommand.append( fullLatexBaseName );
        deleteCommand.append( ".dvi " );
        deleteCommand.append( fullLatexBaseName );
        deleteCommand.append( ".log " );
        deleteCommand.append( fullLatexBaseName );
        deleteCommand.append( ".ps " );
        deleteCommand.append( fullLatexBaseName );
        deleteCommand.append( ".tex" );
        /* "rm LHPC_SpectrumPlotter_gnuplot.dat \
         * LHPC_SlhaPlotter_gnuplot.input \
         * LHPC_SlhaPlotter_gnuplot_TeX.eps \
         * LHPC_SlhaPlotter_gnuplot_TeX.tex \
         * LHPC_SlhaPlotter_LaTeX.aux \
         * LHPC_SlhaPlotter_LaTeX.dvi \
         * LHPC_SlhaPlotter_LaTeX.log \
         * LHPC_SlhaPlotter_LaTeX.ps \
         * LHPC_SlhaPlotter_LaTeX.tex"
         */
        moveCommand.append( " " );
        moveCommand.append( fullLatexBaseName );
        moveCommand.append( ".eps" );
        moveCommand.append( " " );

        fullColumnWidth = ( labelWidth
                            + joinerWidth
                            + flatBitWidth
                            + columnPairOffset
                            + joinerWidth
                            + labelWidth );
      }
      else
      {
        throw std::runtime_error( "Could not parse <SpecificInput>." );
      }
    }

    // This sets up the maximum vertical range, and takes care of labels for
    // lines which would be above the upper limit.
    void SlhaValuePlotter::setUpVerticalRange( double const largestValueFound )
    {
      // If the maximum range is to be automatically calculated, the value for
      // verticalAxisUpperRange was given as zero or a negative number.
      if( verticalAxisUpperRange <= 0.0 )
      {
        verticalAxisUpperRange = ( largestValueFound * automaticScaleFactor );
      }
      else
      {
        for( std::vector< std::list< SlhaValuePlotLine > >::iterator
             whichColumn( columnSet.begin() );
             whichColumn < columnSet.end();
             ++whichColumn )
        {
          for( std::list< SlhaValuePlotLine >::iterator
               whichLine( whichColumn->begin() );
               whichLine != whichColumn->end();
               ++whichLine )
          {
            if( whichLine->getValue() > verticalAxisUpperRange )
            {
              whichLine->relabelForAboveRange( verticalAxisUpperRange );
            }
          }
        }
      }

      // Either way, the label height is set as the fraction times the range.
      labelHeight = ( labelHeightFractionOfRange * verticalAxisUpperRange );
      labelShuffleDistance = ( labelSeparationShuffleFactor * labelHeight );
    }

    // This tries to move all the labels up or down until they are all
    // separated by ( labelHeightFractionOfRange * verticalAxisUpperRange ).
    void SlhaValuePlotter::floatLabels()
    {
      for( std::vector< std::list< SlhaValuePlotLine > >::iterator
           whichColumn( columnSet.begin() );
           whichColumn < columnSet.end();
           ++whichColumn )
      {
        // Each column has 2 sets of labels to float independently, so pointers
        // to the SlhaValuePlotLine objects which have labels on the left are
        // sorted into one vector and on the right to another vector, and each
        // of those vectors is the used to float the labels off each other.
        std::vector< SlhaValuePlotLine* > leftLabels;
        std::vector< SlhaValuePlotLine* > rightLabels;
        for( std::list< SlhaValuePlotLine >::iterator
             columnLine( whichColumn->begin() );
             columnLine != whichColumn->end();
             ++columnLine )
        {
          if( columnLine->hasLabelOnLeftOfColumn() )
          {
            leftLabels.push_back( &(*columnLine) );
          }
          else
          {
            rightLabels.push_back( &(*columnLine) );
          }
        }
        floatLabelsFromOneSide( leftLabels );
        floatLabelsFromOneSide( rightLabels );
      }  // End of loop over columns.
    }


    // This prepares input files for gnuplot to plot the lines.
    void SlhaValuePlotter::prepareTemporaryInputFiles()
    {
      std::string const
      gnuplotDataFileName( "LHPC_SlhaPlotter_gnuplot.dat" );
      std::string const
      gnuplotCommandFileName( "LHPC_SlhaPlotter_gnuplot.input" );
      std::string const
      gnuplotTexBaseName( "LHPC_SlhaPlotter_gnuplot_TeX" );

      std::ofstream gnuplotDataFile( gnuplotDataFileName.c_str() );
      std::ofstream gnuplotCommandFile( gnuplotCommandFileName.c_str() );
      if( !( gnuplotDataFile.is_open()
             &&
             gnuplotCommandFile.is_open() ) )
      {
        throw std::runtime_error( "Could not open "
                                  + gnuplotDataFileName
                                  + " and "
                                  + gnuplotCommandFileName
                                  + " as temporary files!" );
      }

      gnuplotCommandFile
      << "set term epslatex color solid"
      << std::endl << "set output \"" << gnuplotTexBaseName << ".eps\""
      << std::endl << "set format x \"\""
      << std::endl << "set xlabel \"" << horizontalAxisLabel << "\""
      << std::endl << "set ylabel \"" << verticalAxisLabel << "\""
      << std::endl << "unset xtics"
      << std::endl << "set ytics out"
      << std::endl;

      int gnuplotLineIndex( 0 );

      for( int columnIndex( columnSet.size() - 1);
           columnIndex >= 0;
           --columnIndex )
      {
        for( std::list< SlhaValuePlotLine >::iterator
             columnLine( columnSet[ columnIndex ].begin() );
             columnLine != columnSet[ columnIndex ].end();
             ++columnLine )
        {
          columnLine->writeLineData( gnuplotLineIndex,
                                     gnuplotDataFile,
                                     gnuplotCommandFile );
          columnLine->writeLabel( gnuplotCommandFile );
        }
      }

      gnuplotDataFile.close();

      // Finally the plot command is constructed.
      gnuplotCommandFile
      << "plot [0:"
      << ( ( columnSet.size() * fullColumnWidth ) + 2.0 * marginWidth )
      << "] [0:" << verticalAxisUpperRange << "] '" << gnuplotDataFileName
      << "' index 0 notitle with lines ls 1";
      for( int lineIndex( 1 );
           lineIndex < gnuplotLineIndex;
           ++lineIndex )
      {
        gnuplotCommandFile
        << ", '" << gnuplotDataFileName
        << "' index " << lineIndex
        << " notitle with lines ls " << ( lineIndex + 1 );
      }
      gnuplotCommandFile << std::endl;
      gnuplotCommandFile.close();

      std::ofstream latexFile( ( fullLatexBaseName + ".tex" ).c_str() );
      if( !(latexFile.is_open()) )
      {
        throw std::runtime_error( "Could not open "
                                  + fullLatexBaseName
                                  + ".tex as a temporary file!" );
      }
      latexFile
      << "\\documentclass{article}" << std::endl
      << "\\usepackage{graphics}" << std::endl
      << "\\newlength\\labelmover" << std::endl
      << "\\begin{document}" << std::endl
      << "\\pagestyle{empty}" << std::endl
      << "\\begin{center}" << std::endl
      << "\\input{./" << gnuplotTexBaseName << ".tex}" << std::endl
      << "\\end{center}" << std::endl
      << "\\end{document}" << std::endl;
      latexFile.close();
    }

    // This executes the commands to run gnuplot etc. on the prepared input
    // files.
    void SlhaValuePlotter::executeCommands( std::string plotFilename,
                                            bool const shouldCleanUp )
    {
      int systemReturn( 0 );
      systemReturn = system( gnuplotCommand.c_str() );
      // "gnuplot LHPC_SpectrumPlotter_gnuplot.input"
      systemReturn = system( latexCommand.c_str() );
      // "latex LHPC_SpectrumPlotter_LaTeX.tex"
      systemReturn = system( dvipsCommand.c_str() );
      // "dvips LHPC_SpectrumPlotter_LaTeX.dvi"
      systemReturn = system( ps2epsCommand.c_str() );
      // "ps2eps LHPC_SpectrumPlotter_LaTeX.ps"
      if( shouldCleanUp )
      {
        systemReturn = system( deleteCommand.c_str() );
        /* "rm LHPC_SpectrumPlotter_gnuplot.dat \
         * LHPC_SpectrumPlotter_gnuplot.input \
         * LHPC_SpectrumPlotter_gnuplot_TeX.eps \
         * LHPC_SpectrumPlotter_gnuplot_TeX.tex \
         * LHPC_SpectrumPlotter_LaTeX.aux \
         * LHPC_SpectrumPlotter_LaTeX.dvi \
         * LHPC_SpectrumPlotter_LaTeX.log \
         * LHPC_SpectrumPlotter_LaTeX.ps \
         * LHPC_SpectrumPlotter_LaTeX.tex"
         */
      }
      systemReturn = system( (moveCommand + plotFilename).c_str() );
      // "mv LHPC_SpectrumPlotter_LaTeX.eps plotFileName"

      if( systemReturn == -1 )
      {
        std::cout
        << std::endl
        << "Call to \"system\" returned -1. No idea what to do with that.";
        std::cout << std::endl;

      }
    }

    // This adds a new column to columnSet based on the given definition.
    void
    SlhaValuePlotter::addNewColumn( std::string columnDefinitionXml,
                                  LHPC::SlhaSimplisticInterpreter& slhaParser )
    {
      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.loadString( columnDefinitionXml ) )
      {
        columnSet.push_back( std::list< SlhaValuePlotLine >() );
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "ColumnLine" ) )
          {
            addColumnLine( xmlParser.getTrimmedCurrentElementContent(),
                           slhaParser);
          }
        }
      }
      else
      {
        throw std::runtime_error( "Could not parse <ColumnDefinition>." );
      }
    }

    // This adds a new column line to the last column in columnSet based on the
    // given definition.
    void
    SlhaValuePlotter::addColumnLine( std::string lineDefinitionXml,
                                  LHPC::SlhaSimplisticInterpreter& slhaParser )
    {
      double slhaValue( 0.0 );
      std::string labelString( "" );
      bool labelLeftOfColumn( true );
      SlhaValueLineColoring const* lineColoring( NULL );

      BOL::AsciiXmlParser xmlParser;
      if( xmlParser.loadString( lineDefinitionXml ) )
      {
        while( xmlParser.readNextElement() )
        {
          if( xmlParser.currentElementNameMatches( "SlhaValue" ) )
          {
            slhaValue = slhaParser.getDouble(
                                 xmlParser.getTrimmedCurrentElementContent() );
            if( useAbsoluteValues
                &&
                ( slhaValue < 0.0 ) )
            {
              slhaValue = -slhaValue;
            }
            slhaValue *= verticalAxisScaling;
          }
          else if( xmlParser.currentElementNameMatches( "LabelLatex" ) )
          {
            labelString = xmlParser.getTrimmedCurrentElementContent();
          }
          else if( xmlParser.currentElementNameMatches( "LabelSide" ) )
          {
            std::string
            whichSide( xmlParser.getTrimmedCurrentElementContent() );
            BOL::StringParser::transformToLowercase( whichSide );

            labelLeftOfColumn = ( whichSide.compare( "left" ) == 0 );
          }
          else if( xmlParser.currentElementNameMatches( "LineColor" ) )
          {
            std::map< std::string, std::string > const&
            lineAttributes( xmlParser.getCurrentElementAttributes() );
            std::map< std::string, std::string >::const_iterator
            attributeFinder( lineAttributes.find( "DrawType" ) );
            if( attributeFinder != lineAttributes.end() )
            {
              lineColoring = getNewLineColoring( attributeFinder->second,
                                   xmlParser.getTrimmedCurrentElementContent(),
                                                 slhaParser );
            }
            else
            {
              throw std::runtime_error(
                                      "Unknown \"DrawType\" in <LineColor>." );
            }
          }
        }
      }
      double const columnCenterPosition( marginWidth
                                         + ( ( columnSet.size() - 0.5 )
                                             * fullColumnWidth ) );
      columnSet.back().push_back( SlhaValuePlotLine( slhaValue,
                                                     columnCenterPosition,
                                                     flatBitWidth,
                                                     columnPairOffset,
                                                     joinerWidth,
                                                     labelString,
                                                     labelLeftOfColumn,
                                                     lineColoring ) );
    }

    // This tries to move all the labels up or down until they are all
    // separated by ( labelHeightFractionOfRange * verticalAxisUpperRange ).
    void SlhaValuePlotter::floatLabelsFromOneSide(
                                  std::vector< SlhaValuePlotLine* >& labelSet )
    {
      // The labels only need to be checked for overlap if there are at least 2
      // labels.
      if( labelSet.size() > 1 )
      {
        // At least one iteration of the shuffling loop must be performed so
        // that at least the labels are checked to be in position.
        bool notYetFinishedShuffling( true );
        int remainingShuffles( maximumLabelFloatingShuffles );
        // Each column is allowed to try floating the labels apart
        // maximumLabelFloatingShuffles times.

        std::vector< SlhaValuePlotLine* >::iterator
        lowerLine( labelSet.begin() );
        std::vector< SlhaValuePlotLine* >::iterator upperLine( lowerLine );

        while( notYetFinishedShuffling
               &&
               ( remainingShuffles  > 0 ) )
        {
          notYetFinishedShuffling = false;
          // each iteration starts by assuming that it just finds that all
          // the labels are sufficiently separate, with no floating needed.
          --remainingShuffles;
          lowerLine = labelSet.begin();
          // if the lowest label is too low, it is bumped up to its lowest
          // allowed position:
          if( (*lowerLine)->getLabelPosition() < labelHeight )
          {
            notYetFinishedShuffling = true;
            // moving the lowest label means that the column has to be
            // checked again for overlapping labels.
            (*lowerLine)->setLabelPosition( labelHeight );
          }
          upperLine = lowerLine;
          ++upperLine;
          // The lower line iterator is always 1 step behind the upper line
          // iterator.
          while( upperLine < labelSet.end() )
          {
            if( ( (*upperLine)->getLabelPosition()
                  - (*lowerLine)->getLabelPosition() ) < labelHeight )
            {
              // If a pair of overlapping labels is found,
              // notYetFinishedShuffling is set to true since this iteration
              // did not merely confirm that the labels are in order and well
              // separated.
              notYetFinishedShuffling = true;

              // Now the labels are floated away from each other around their
              // average position. Since the conditional does not check for
              // the absolute separation, if labels ever get twisted past
              // each other, the conditional will find them, and since the
              // lines are sorted according to their SLHA value rather than
              // their current label position, the labels get separated in the
              // right direction.
              double const
              labelAverage( 0.5 * ( (*upperLine)->getLabelPosition()
                                    + (*lowerLine)->getLabelPosition() ) );
              (*lowerLine)->setLabelPosition( labelAverage
                                              - labelShuffleDistance );
              (*upperLine)->setLabelPosition( labelAverage
                                              + labelShuffleDistance );
            }
            ++lowerLine;
            ++upperLine;
          }
          // Now upperLine is at columnPointer->end(), so lowerLine is at the
          // highest label, which is ensured to be lower than its allowed
          // highest value:
          if( (*lowerLine)->getLabelPosition() > ( verticalAxisUpperRange
                                                   - labelHeight ) )
          {
            // Moving the highest label means that the column has to be
            // checked again for overlapping labels.
            notYetFinishedShuffling = true;
            (*lowerLine)->setLabelPosition( verticalAxisUpperRange
                                            - labelHeight );
          }
        }  // End of loop over shuffles.
      }  // End of if the column had 2 or more labels on the relevant side.
    }

  } /* namespace SLHA */
} /* namespace LHPC */
