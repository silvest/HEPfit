/*
 * NeutralinosOneToFive.cpp
 *
 *  Created on: Jan 27, 2012
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
    NeutralinosOneToFive::NeutralinosOneToFive( bool const isVerbose,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        NeutralinosOneToFour( isVerbose,
                              defaultFlags ),
        neutralinoFive( PDGIX::neutralinoFive,
                        PDGVII::neutralinoFive,
                        mapAndVectorAndBools,
                        true,
                        "N5",
                        "${\\tilde{{\\tilde{{\\chi}}}}}_{5}^{0}$" )
    {
      neutralinoPointers.push_back( &neutralinoFive );
    }

    NeutralinosOneToFive::~NeutralinosOneToFive()
    {
      // does nothing.
    }

  }

}
