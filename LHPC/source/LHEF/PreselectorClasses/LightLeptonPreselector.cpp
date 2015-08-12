/*
 * LightLeptonPreselector.cpp
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
      LightLeptonPreselector::LightLeptonPreselector( LhefParser& lhefParser,
                                            double const transverseMomentumCut,
                                             double const pseudorapidityCut ) :
          StandardPreselector( lhefParser,
                               BOL::Vi( PDGIX::positronOne
                                     )( -PDGIX::positronOne
                                     )( PDGIX::positronTwo
                                     )( -PDGIX::positronTwo
                                     )( PDGVII::positronOne
                                     )( -PDGVII::positronOne
                                     )( PDGVII::positronTwo
                                     ).e( -PDGVII::positronTwo ),
                               transverseMomentumCut,
                               pseudorapidityCut )
      {
        // just an initialization list.
      }

      LightLeptonPreselector::~LightLeptonPreselector()
      {
        // does nothing.
      }

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
