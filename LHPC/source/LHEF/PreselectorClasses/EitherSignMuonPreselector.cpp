/*
 * EitherSignMuonPreselector.cpp
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
      EitherSignMuonPreselector::EitherSignMuonPreselector(
                                                        LhefParser& lhefParser,
                                            double const transverseMomentumCut,
                                             double const pseudorapidityCut ) :
          StandardPreselector( lhefParser,
                               BOL::Vi( PDGIX::positronTwo
                                     )( -PDGIX::positronTwo
                                     )( PDGVII::positronTwo
                                     ).e( -PDGVII::positronTwo ),
                               transverseMomentumCut,
                               pseudorapidityCut )
      {
        // just an initialization list.
      }

      EitherSignMuonPreselector::~EitherSignMuonPreselector()
      {
        // does nothing.
      }

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
