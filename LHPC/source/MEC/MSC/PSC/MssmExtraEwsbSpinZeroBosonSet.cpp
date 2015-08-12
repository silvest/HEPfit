/*
 * MssmExtraEwsbSpinZeroBosonSet.cpp
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
    MssmExtraEwsbSpinZeroBosonSet::MssmExtraEwsbSpinZeroBosonSet(
                                                          bool const isVerbose,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        neutralColorlessScalarTwo( PDGIX::neutralColorlessScalarTwo,
                                   PDGVII::neutralColorlessScalarTwo,
                                   mapAndVectorAndBools,
                                   true,
                                   "h02",
                                   "$h_{2}^{0}$" ),
        neutralColorlessPseudoscalarOne(
                                        PDGIX::neutralColorlessPseudoscalarOne,
                                       PDGVII::neutralColorlessPseudoscalarOne,
                                         mapAndVectorAndBools,
                                         true,
                                         "A01",
                                         "$A_{1}^{0}$" ),
        positiveColorlessSpinZeroBosonOne(
                                      PDGIX::positiveColorlessSpinZeroBosonOne,
                                     PDGVII::positiveColorlessSpinZeroBosonOne,
                                           mapAndVectorAndBools,
                                           false,
                                           "hp1",
                                           "$h_{1}^{+}$" ),
        negativeColorlessSpinZeroBosonOne( positiveColorlessSpinZeroBosonOne,
                                           "hm1",
                                           "$h_{1}^{-}$" ),
        neutralScalarsAndPseudoscalarPointers(),
        chargedColorlessSpinZeroBosonPointers( 2,
                                        &positiveColorlessSpinZeroBosonOne ),
        ewsbSpinZeroAndOneBosonPointers(),
        neutralEwsbSpinZeroAndOneBosonPointers(),
        chargedEwsbSpinZeroAndOneBosonPointers(
                                        chargedColorlessSpinZeroBosonPointers )
    {
      // neutralScalarsAndPseudoscalarPointers should start with the scalar
      // from a StandardModel instance.
      MassEigenstate* vectorFiller( MassEigenstate::findPointerWithCode(
                                              PDGIX::neutralColorlessScalarOne,
                                                                pdgCodeMap ) );
      if( NULL == vectorFiller )
      {
        vectorFiller = MassEigenstate::findPointerWithCode(
                                             PDGVII::neutralColorlessScalarOne,
                                                            pdgCodeMap );
      }
      if( NULL != vectorFiller )
      {
        neutralScalarsAndPseudoscalarPointers.push_back( vectorFiller );
      }
      neutralScalarsAndPseudoscalarPointers.push_back(
                                                  &neutralColorlessScalarTwo );
      neutralScalarsAndPseudoscalarPointers.push_back(
                                            &neutralColorlessPseudoscalarOne );
      chargedColorlessSpinZeroBosonPointers[ 1 ]
      = &negativeColorlessSpinZeroBosonOne;
      chargedEwsbSpinZeroAndOneBosonPointers[ 1 ]
      = chargedColorlessSpinZeroBosonPointers[ 1 ];
      vectorFiller = MassEigenstate::findPointerWithCode( PDGIX::wPlusBosonOne,
                                                          pdgCodeMap );
      if( NULL == vectorFiller )
      {
        vectorFiller
        = MassEigenstate::findPointerWithCode( PDGVII::wPlusBosonOne,
                                               pdgCodeMap );
      }
      if( NULL != vectorFiller )
      {
        chargedEwsbSpinZeroAndOneBosonPointers.push_back( vectorFiller );
        chargedEwsbSpinZeroAndOneBosonPointers.push_back(
                                       &(vectorFiller->getChargeConjugate()) );
      }
      neutralEwsbSpinZeroAndOneBosonPointers
      = neutralScalarsAndPseudoscalarPointers;
      vectorFiller = MassEigenstate::findPointerWithCode( PDGIX::zBosonOne,
                                                          pdgCodeMap );
      if( NULL == vectorFiller )
      {
        vectorFiller
        = MassEigenstate::findPointerWithCode( PDGVII::zBosonOne,
                                               pdgCodeMap );
      }
      if( NULL != vectorFiller )
      {
        neutralEwsbSpinZeroAndOneBosonPointers.push_back( vectorFiller );
      }
      ewsbSpinZeroAndOneBosonPointers = neutralEwsbSpinZeroAndOneBosonPointers;
      ewsbSpinZeroAndOneBosonPointers.insert(
                                         ewsbSpinZeroAndOneBosonPointers.end(),
                                chargedEwsbSpinZeroAndOneBosonPointers.begin(),
                                chargedEwsbSpinZeroAndOneBosonPointers.end() );
    }

    MssmExtraEwsbSpinZeroBosonSet::~MssmExtraEwsbSpinZeroBosonSet()
    {
      // does nothing.
    }

  }

}
