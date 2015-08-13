/*
 * EitherSignMuonPreselector.hpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef EITHERSIGNMUONPRESELECTOR_HPP_
#define EITHERSIGNMUONPRESELECTOR_HPP_

#include "StandardPreselector.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      // this selects final-state muons & antimuons.
      class EitherSignMuonPreselector : public StandardPreselector
      {
      public:
        EitherSignMuonPreselector( LhefParser& lhefParser,
                                   double const transverseMomentumCut = 0.0,
                                   double const pseudorapidityCut = -1.0 );
        virtual
        ~EitherSignMuonPreselector();
      };

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
#endif /* EITHERSIGNMUONPRESELECTOR_HPP_ */
