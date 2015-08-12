/*
 * JetPreselector.hpp
 *
 *  Created on: Jan 25, 2013
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 */

#ifndef JETPRESELECTOR_HPP_
#define JETPRESELECTOR_HPP_

#include "StandardPreselector.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace PreselectorClass
    {
      // this class selects final-state gluon, quarks, & anti-quarks, except
      // tops or antitops, & bottoms & antibottoms are included if
      // includeBottoms is true.
      class JetPreselector : public StandardPreselector
      {
      public:
        JetPreselector( LhefParser& lhefParser,
                        bool const includeBottoms = true,
                        double const transverseMomentumCut = 0.0,
                        double const pseudorapidityCut = -1.0 );
        virtual
        ~JetPreselector();

      protected:
        static std::vector< int >
        jetCodes( bool const includeBottoms );
      };

    } /* namespace PreselectorClass */
  } /* namespace LHEF */
} /* namespace LHPC */
#endif /* JETPRESELECTOR_HPP_ */
