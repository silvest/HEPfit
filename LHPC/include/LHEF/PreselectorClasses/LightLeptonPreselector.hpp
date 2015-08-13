/*
 * LightLeptonPreselector.hpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef LIGHTLEPTONPRESELECTOR_HPP_
#define LIGHTLEPTONPRESELECTOR_HPP_

#include "StandardPreselector.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      // this selects electrons, positrons, muons, & antimuons.
      class LightLeptonPreselector : public StandardPreselector
      {
      public:
        LightLeptonPreselector( LhefParser& lhefParser,
                                double const transverseMomentumCut = 0.0,
                                double const pseudorapidityCut = -1.0 );
        virtual
        ~LightLeptonPreselector();
      };

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
#endif /* LIGHTLEPTONPRESELECTOR_HPP_ */
