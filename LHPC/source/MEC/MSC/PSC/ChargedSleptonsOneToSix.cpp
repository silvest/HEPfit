/*
 * ChargedSleptonsOneToSix.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "MEC.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    ChargedSleptonsOneToSix::ChargedSleptonsOneToSix( bool const isVerbose,
                                                   bool const flavorConserving,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        spositronOne( PDGIX::spositronOne,
                      -PDGVII::selectronOne,
                      mapAndVectorAndBools,
                      false,
                      "se1c",
                      "${\\tilde{e}}_{1}^{\\ast}$" ),
        antispositronOne( spositronOne,
                          "se1",
                          "${\\tilde{e}}_{1}$" ),
        spositronTwo( PDGIX::spositronTwo,
                      -PDGVII::selectronTwo,
                      mapAndVectorAndBools,
                      false,
                      "se2c",
                      "${\\tilde{e}}_{2}^{\\ast}$" ),
        antispositronTwo( spositronTwo,
                          "se2",
                          "${\\tilde{e}}_{2}$" ),
        spositronThree( PDGIX::spositronThree,
                        -PDGVII::selectronThree,
                        mapAndVectorAndBools,
                        false,
                        "se3c",
                        "${\\tilde{e}}_{3}^{\\ast}$" ),
        antispositronThree( spositronThree,
                            "se3",
                            "${\\tilde{e}}_{3}$" ),
        spositronFour( PDGIX::spositronFour,
                       -PDGVII::selectronFour,
                       mapAndVectorAndBools,
                       false,
                       "se4c",
                       "${\\tilde{e}}_{4}^{\\ast}$" ),
        antispositronFour( spositronFour,
                           "se4",
                           "${\\tilde{e}}_{4}$" ),
        spositronFive( PDGIX::spositronFive,
                       -PDGVII::selectronFive,
                       mapAndVectorAndBools,
                       false,
                       "se5c",
                       "${\\tilde{e}}_{5}^{\\ast}$" ),
        antispositronFive( spositronFive,
                           "se5",
                           "${\\tilde{e}}_{5}$" ),
        spositronSix( PDGIX::spositronSix,
                      -PDGVII::selectronSix,
                      mapAndVectorAndBools,
                      false,
                      "se6c",
                      "${\\tilde{e}}_{6}^{\\ast}$" ),
        antispositronSix( spositronSix,
                          "se6",
                          "${\\tilde{e}}_{6}$" ),
        positiveSleptonPointers( 6,
                                 &spositronOne ),
        negativeSleptonPointers( 6,
                                 &antispositronOne )
    {
      positiveSleptonPointers[ 1 ] = &spositronTwo;
      positiveSleptonPointers[ 2 ] = &spositronThree;
      positiveSleptonPointers[ 3 ] = &spositronFour;
      positiveSleptonPointers[ 4 ] = &spositronFive;
      positiveSleptonPointers[ 5 ] = &spositronSix;
      negativeSleptonPointers[ 1 ] = &antispositronTwo;
      negativeSleptonPointers[ 2 ] = &antispositronThree;
      negativeSleptonPointers[ 3 ] = &antispositronFour;
      negativeSleptonPointers[ 4 ] = &antispositronFive;
      negativeSleptonPointers[ 5 ] = &antispositronSix;
      if( flavorConserving )
      {
        spositronOne.setAsciiName( "seLc" );
        spositronOne.setLatexName( "${\\tilde{e}}_{L}^{\\ast}$" );
        antispositronOne.setAsciiName( "seL" );
        antispositronOne.setLatexName( "${\\tilde{e}}_{L}$" );
        spositronTwo.setAsciiName( "smuLc" );
        spositronTwo.setLatexName( "${\\tilde{\\mu}}_{L}^{\\ast}$" );
        antispositronTwo.setAsciiName( "smuL" );
        antispositronTwo.setLatexName( "${\\tilde{\\mu}}_{L}$" );
        spositronThree.setAsciiName( "sta1c" );
        spositronThree.setLatexName( "${\\tilde{\\tau}}_{1}^{\\ast}$" );
        antispositronThree.setAsciiName( "sta1" );
        antispositronThree.setLatexName( "${\\tilde{\\tau}}_{1}$" );
        spositronFour.setAsciiName( "seRc" );
        spositronFour.setLatexName( "${\\tilde{e}}_{R}^{\\ast}$" );
        antispositronFour.setAsciiName( "seR" );
        antispositronFour.setLatexName( "${\\tilde{e}}_{R}$" );
        spositronFive.setAsciiName( "smuRc" );
        spositronFive.setLatexName( "${\\tilde{\\mu}}_{R}^{\\ast}$" );
        antispositronFive.setAsciiName( "smuR" );
        antispositronFive.setLatexName( "${\\tilde{\\mu}}_{R}$" );
        spositronSix.setAsciiName( "sta2c" );
        spositronSix.setLatexName( "${\\tilde{\\tau}}_{2}^{\\ast}$" );
        antispositronSix.setAsciiName( "sta2" );
        antispositronSix.setLatexName( "${\\tilde{\\tau}}_{2}$" );
      }
    }

    ChargedSleptonsOneToSix::~ChargedSleptonsOneToSix()
    {
      // does nothing.
    }

  }

}
