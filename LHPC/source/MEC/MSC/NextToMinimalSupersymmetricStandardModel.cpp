/*
 * NextToMinimalSupersymmetricStandardModel.cpp
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
    NMSSM::NextToMinimalSupersymmetricStandardModel( bool const isVerbose,
                                               bool const neutrinosAreMajorana,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        MinimalSupersymmetricStandardModel( isVerbose,
                                            neutrinosAreMajorana,
                                            defaultFlags ),
        NmssmExtraEwsbSpinZeroBosonSet( isVerbose,
                                        defaultFlags ),
        NeutralinosOneToFive( isVerbose,
                              defaultFlags )
    {
      // just an initialization list.
    }

    NMSSM::~NextToMinimalSupersymmetricStandardModel()
    {
      // does nothing.
    }

  }

}
