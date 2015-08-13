/*
 * GluinoOneGeneration.cpp
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
    GluinoOneGeneration::GluinoOneGeneration( bool const isVerbose,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        gluinoOne( PDGIX::gluinoFermion,
                   PDGVII::gluinoFermion,
                   mapAndVectorAndBools,
                   true,
                   "go",
                   "${\\tilde{g}}$" )
    {
      // just an initialization list.
    }

    GluinoOneGeneration::~GluinoOneGeneration()
    {
      // does nothing.
    }

  }

}
