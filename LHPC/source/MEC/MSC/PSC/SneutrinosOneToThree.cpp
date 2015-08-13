/*
 * SneutrinosOneToThree.cpp
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
    SneutrinosOneToThree::SneutrinosOneToThree( bool const isVerbose,
                                                bool const flavorConserving,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        antisneutrinoOne( PDGIX::antisneutrinoOne,
                          -PDGVII::sneutrinoOne,
                          mapAndVectorAndBools,
                          false,
                          "sv1c",
                          "${\\tilde{{\\nu}}}_{1}^{\\ast}$" ),
        sneutrinoOne( antisneutrinoOne,
                      "sv1",
                      "${\\tilde{{\\nu}}}_{1}$" ),
        antisneutrinoTwo( PDGIX::antisneutrinoTwo,
                          -PDGVII::sneutrinoTwo,
                          mapAndVectorAndBools,
                          false,
                          "sv2c",
                          "${\\tilde{{\\nu}}}_{2}^{\\ast}$" ),
        sneutrinoTwo( antisneutrinoTwo,
                      "sv2",
                      "${\\tilde{{\\nu}}}_{2}$" ),
        antisneutrinoThree( PDGIX::antisneutrinoThree,
                            -PDGVII::sneutrinoThree,
                            mapAndVectorAndBools,
                            false,
                            "sv3c",
                            "${\\tilde{{\\nu}}}_{3}^{\\ast}$" ),
        sneutrinoThree( antisneutrinoThree,
                        "sv3",
                        "${\\tilde{{\\nu}}}_{3}$" ),
        antisneutrinoPointers( 3,
                               &antisneutrinoOne ),
        sneutrinoPointers( 3,
                           &sneutrinoOne )
    {
      antisneutrinoPointers[ 1 ] = &antisneutrinoTwo;
      antisneutrinoPointers[ 2 ] = &antisneutrinoThree;
      sneutrinoPointers[ 1 ] = &sneutrinoTwo;
      sneutrinoPointers[ 2 ] = &sneutrinoThree;
      if( flavorConserving )
      {
        antisneutrinoOne.setAsciiName( "svec" );
        antisneutrinoOne.setLatexName( "${\\tilde{{\\nu}}}_{eL}^{\\ast}$" );
        sneutrinoOne.setAsciiName( "sve" );
        sneutrinoOne.setLatexName( "${\\tilde{{\\nu}}}_{eL}$" );
        antisneutrinoTwo.setAsciiName( "svmc" );
        antisneutrinoTwo.setLatexName(
                                     "${\\tilde{{\\nu}}}_{{\\mu}L}^{\\ast}$" );
        sneutrinoTwo.setAsciiName( "svm" );
        sneutrinoTwo.setLatexName( "${\\tilde{{\\nu}}}_{{\\mu}L}$" );
        antisneutrinoThree.setAsciiName( "svtc" );
        antisneutrinoThree.setLatexName(
                                    "${\\tilde{{\\nu}}}_{{\\tau}L}^{\\ast}$" );
        sneutrinoThree.setAsciiName( "svt" );
        sneutrinoThree.setLatexName( "${\\tilde{{\\nu}}}_{{\\tau}L}$" );
      }
    }

    SneutrinosOneToThree::~SneutrinosOneToThree()
    {
      // does nothing.
    }

  }

}
