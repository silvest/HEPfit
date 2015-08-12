/*
 * NmssmExtraEwsbSpinZeroBosonSet.cpp
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
    NmssmExtraEwsbSpinZeroBosonSet::NmssmExtraEwsbSpinZeroBosonSet(
                                                          bool const isVerbose,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        MssmExtraEwsbSpinZeroBosonSet( isVerbose,
                                       defaultFlags ),
        neutralColorlessScalarThree( PDGIX::neutralColorlessScalarThree,
                                     PDGVII::nmssmHiggsScalarThree,
                                     mapAndVectorAndBools,
                                     true,
                                     "h03",
                                     "$h_{3}^{0}$" ),
        neutralColorlessPseudoscalarTwo(
                                        PDGIX::neutralColorlessPseudoscalarTwo,
                                         PDGVII::nmssmHiggsPseudoscalarTwo,
                                         mapAndVectorAndBools,
                                         true,
                                         "A02",
                                         "$A_{2}^{0}$" )
    {
      neutralScalarsAndPseudoscalarPointers.push_back(
                                                &neutralColorlessScalarThree );
      neutralScalarsAndPseudoscalarPointers.push_back(
                                            &neutralColorlessPseudoscalarTwo );
      ewsbSpinZeroAndOneBosonPointers.push_back(
                                                &neutralColorlessScalarThree );
      ewsbSpinZeroAndOneBosonPointers.push_back(
                                            &neutralColorlessPseudoscalarTwo );
    }

    NmssmExtraEwsbSpinZeroBosonSet::~NmssmExtraEwsbSpinZeroBosonSet()
    {
      // does nothing.
    }

  }

}
