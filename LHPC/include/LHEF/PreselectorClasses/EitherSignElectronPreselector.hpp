/*
 * EitherSignElectronPreselector.hpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef EITHERSIGNELECTRONPRESELECTOR_HPP_
#define EITHERSIGNELECTRONPRESELECTOR_HPP_

#include "StandardPreselector.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      // this selects final-state electrons & positrons.
      class EitherSignElectronPreselector : public StandardPreselector
      {
      public:
        EitherSignElectronPreselector( LhefParser& lhefParser,
                                      double const transverseMomentumCut = 0.0,
                                       double const pseudorapidityCut = -1.0 );
        virtual
        ~EitherSignElectronPreselector();
      };

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
#endif /* EITHERSIGNELECTRONPRESELECTOR_HPP_ */
