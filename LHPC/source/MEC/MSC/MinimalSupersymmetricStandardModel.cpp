/*
 * MinimalSupersymmetricStandardModel.cpp
 *
 *  Created on: Jan 8, 2012
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
    MinimalSupersymmetricStandardModel::MinimalSupersymmetricStandardModel(
                                                          bool const isVerbose,
                                               bool const neutrinosAreMajorana,
                                                   bool const flavorConserving,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        StandardModel( isVerbose,
                       neutrinosAreMajorana,
                       defaultFlags ),
        MssmExtraEwsbSpinZeroBosonSet( isVerbose,
                                       defaultFlags ),
        ChargedSleptonsOneToSix( isVerbose,
                                 flavorConserving,
                                 defaultFlags ),
        SneutrinosOneToThree( isVerbose,
                              flavorConserving,
                              defaultFlags ),
        SquarksOneToSix( isVerbose,
                         flavorConserving,
                         defaultFlags ),
        NeutralinosOneToFour( isVerbose,
                              defaultFlags ),
        CharginosOneToTwo( isVerbose,
                           defaultFlags ),
        GluinoOneGeneration( isVerbose,
                             defaultFlags )
    {
      neutralinoOne.setFlags( &defaultEscapesDetectorBoolVector );
    }

    MinimalSupersymmetricStandardModel::~MinimalSupersymmetricStandardModel()
    {
      // does nothing.
    }

  }

}
