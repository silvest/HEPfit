/*
 * NeutralinosOneToFour.cpp
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
    NeutralinosOneToFour::NeutralinosOneToFour( bool const isVerbose,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        neutralinoOne( PDGIX::neutralinoOne,
                       PDGVII::neutralinoOne,
                       mapAndVectorAndBools,
                       true,
                       "N1",
                       "${\\tilde{{\\tilde{{\\chi}}}}}_{1}^{0}$" ),
        neutralinoTwo( PDGIX::neutralinoTwo,
                       PDGVII::neutralinoTwo,
                       mapAndVectorAndBools,
                       true,
                       "N2",
                       "${\\tilde{{\\chi}}}_{2}^{0}$" ),
        neutralinoThree( PDGIX::neutralinoThree,
                         PDGVII::neutralinoThree,
                         mapAndVectorAndBools,
                         true,
                         "N3",
                         "${\\tilde{{\\chi}}}_{3}^{0}$" ),
        neutralinoFour( PDGIX::neutralinoFour,
                        PDGVII::neutralinoFour,
                        mapAndVectorAndBools,
                        true,
                        "N4",
                        "${\\tilde{{\\chi}}}_{4}^{0}$" ),
        neutralinoPointers( 4,
                            &neutralinoOne )
    {
     neutralinoPointers[ 1 ] = &neutralinoTwo;
     neutralinoPointers[ 2 ] = &neutralinoThree;
     neutralinoPointers[ 3 ] = &neutralinoFour;
    }

    NeutralinosOneToFour::~NeutralinosOneToFour()
    {
      // does nothing.
    }

  }

}
