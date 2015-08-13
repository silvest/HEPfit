/*
 * EitherSignElectronPreselector.cpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#include "LHEF.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      EitherSignElectronPreselector::EitherSignElectronPreselector(
                                                        LhefParser& lhefParser,
                                            double const transverseMomentumCut,
                                             double const pseudorapidityCut ) :
          StandardPreselector( lhefParser,
                               BOL::Vi( PDGIX::positronOne
                                     )( -PDGIX::positronOne
                                     )( PDGVII::positronOne
                                     ).e( -PDGVII::positronOne ),
                               transverseMomentumCut,
                               pseudorapidityCut )
      {
        // just an initialization list.
      }

      EitherSignElectronPreselector::~EitherSignElectronPreselector()
      {
        // does nothing.
      }

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
