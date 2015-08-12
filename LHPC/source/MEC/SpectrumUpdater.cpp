/*
 * SpectrumUpdater.cpp
 *
 *  Created on: Mar 15, 2012
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
  SpectrumUpdater::SpectrumUpdater() :
    fmassMap( NULL ),
    fmasserrMap( NULL ),
    massMap( NULL ),
    isHoldingDecayFlag( false ),
    decayerCode( 0 ),
    decayWidth( BOL::UsefulStuff::notANumber ),
    decayChannels(),
    branchingRatioAndProducts()
  {
    // just an initialization list.
  }

  SpectrumUpdater::~SpectrumUpdater()
  {
    // does nothing.
  }


  void
  SpectrumUpdater::updateMassEigenstates(
                                MassEigenstateCodeToPointerMap& codeMap ) const
  {
    MassEigenstate* massEigenstateFiller( NULL );
    // 1st the masses from the MASS block are recorded:

    if( NULL != massMap )
    {
      std::map< int, double >::const_iterator
      massMapIterator( massMap->begin() );
      while( massMap->end() != massMapIterator )
      {
        massEigenstateFiller
        = MassEigenstate::findPointerWithCode( massMapIterator->first,
                                               codeMap );
        if( NULL != massEigenstateFiller )
        {
          massEigenstateFiller->recordMass( massMapIterator->second );
          // now the charge-conjugate is considered:
          if( !(massEigenstateFiller->getChargeConjugate(
                                                     ).hasMassBeenRecorded()) )
          {
            massEigenstateFiller->getChargeConjugate().recordMass(
                                                     massMapIterator->second );
          }
        }
        ++massMapIterator;
      }
    }

    // then the masses from the FMASS block are recorded:
    if( NULL != fmassMap )
    {
      std::multimap< int, RunningConstant >::const_iterator
      fmassMapIterator( fmassMap->begin() );
      while( fmassMap->end() != fmassMapIterator )
      {
        massEigenstateFiller
        = MassEigenstate::findPointerWithCode( fmassMapIterator->first,
                                               codeMap );
        if( NULL != massEigenstateFiller )
        {
          massEigenstateFiller->recordMass(
                                           fmassMapIterator->second.getValue(),
                                            0.0,
                                            0.0,
                                          fmassMapIterator->second.getScheme(),
                                         fmassMapIterator->second.getScale() );
          // now the charge-conjugate is considered:
          if( massEigenstateFiller->getChargeConjugate(
                                                ).getAllRecordedMasses().size()
              < massEigenstateFiller->getAllRecordedMasses().size() )
          {
            massEigenstateFiller->getChargeConjugate().recordMass(
                                           fmassMapIterator->second.getValue(),
                                                                   0.0,
                                                                   0.0,
                                          fmassMapIterator->second.getScheme(),
                                         fmassMapIterator->second.getScale() );
          }
        }
        ++fmassMapIterator;
      }
    }

    // then the mass errors from the FMASSERR block are recorded:
    if( NULL != fmasserrMap )
    {
      std::multimap< int, RunningConstantError >::const_iterator
      fmasserrMapIterator( fmasserrMap->begin() );
      while( fmasserrMap->end() != fmasserrMapIterator )
      {
        massEigenstateFiller
        = MassEigenstate::findPointerWithCode( fmasserrMapIterator->first,
                                               codeMap );
        if( NULL != massEigenstateFiller )
        {
          massEigenstateFiller->recordMassError(
                             fmasserrMapIterator->second.getMinusUncertainty(),
                              fmasserrMapIterator->second.getPlusUncertainty(),
                                       fmasserrMapIterator->second.getScheme(),
                                      fmasserrMapIterator->second.getScale() );
          // now the charge-conjugate is considered:
          if( !(massEigenstateFiller->isSelfConjugate()) )
          {
            massEigenstateFiller->getChargeConjugate().recordMassError(
                             fmasserrMapIterator->second.getMinusUncertainty(),
                              fmasserrMapIterator->second.getPlusUncertainty(),
                                       fmasserrMapIterator->second.getScheme(),
                                      fmasserrMapIterator->second.getScale() );
            // recordMassError does nothing if it doesn't find a matching
            // running mass.
          }
        }
        ++fmasserrMapIterator;
      }
    }

    // finally the given decay is recorded:
    if( isHoldingDecayFlag )
    {
      MassEigenstate::MassEigenstatesPairedWithBr decayRecorder;
      massEigenstateFiller
      = MassEigenstate::findPointerWithCode( decayerCode,
                                             codeMap );
      massEigenstateFiller->setDecayWidth( decayWidth );
      size_t numberOfChannels( decayChannels.size() );
      for( size_t whichChannel( 0 );
           numberOfChannels > whichChannel;
           ++whichChannel )
      {
        decayRecorder.clearPointers();
        for( int whichProduct( decayChannels[ whichChannel ].second.size()
                               - 1 );
             0 <= whichProduct;
             --whichProduct )
        {
          decayRecorder.addPointer( MassEigenstate::findPointerWithCode(
                          decayChannels[ whichChannel ].second[ whichProduct ],
                                                                   codeMap ) );
        }
        decayRecorder.setPairedValueAndSortPointers(
                                         decayChannels[ whichChannel ].first );
        massEigenstateFiller->recordDecay( decayRecorder );
      }

      /* next, the mass eigenstate is checked for whether it has a charge
       * conjugate & whether it has already recorded a decay (so that
       * charge-conjugate pairs do not necessarily have the same decay widths
       * & branching ratios: they are initially set to be so when the 1st of a
       * pair is recorded, but if its charge conjugate is found, only the
       * charge conjugate gets its decays overwritten.
       */
      if( !(massEigenstateFiller->isSelfConjugate())
          &&
          !(massEigenstateFiller->getChargeConjugate(
                                                  ).haveDecaysBeenRecorded()) )
      {
        numberOfChannels = massEigenstateFiller->getDecaySet().size();
        for( size_t whichDecay( 0 );
             numberOfChannels > whichDecay;
             ++whichDecay )
        {
          massEigenstateFiller->getChargeConjugate(
                                                ).recordChargeConjugateOfDecay(
                        *(massEigenstateFiller->getDecaySet()[ whichDecay ]) );
        }
        massEigenstateFiller->getChargeConjugate().setDecayWidth( decayWidth );
      }
    }
  }

}
