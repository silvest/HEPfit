/*
 * TransverseMomentumCut.hpp
 *
 *  Created on: Jan 26, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef TRANSVERSEMOMENTUMCUT_HPP_
#define TRANSVERSEMOMENTUMCUT_HPP_

#include "../FilterRule.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace FilterRuleClass
    {
      // this class returns acceptRatherThanReject if it finds that the
      // transverse momentum from lineToCheck is less than momentumCut.
      class TransverseMomentumCut : public FilterRule
      {
      public:
        TransverseMomentumCut( double const momentumCut,
                               bool const acceptRatherThanReject = true );
        virtual
        ~TransverseMomentumCut();

        virtual bool
        operator()( ParticleLine const& lineToCheck ) const;


      protected:
        double momentumCutSquared;
      };



      inline bool
      TransverseMomentumCut::operator()(
                                        ParticleLine const& lineToCheck ) const
      {
        if( momentumCutSquared <= lineToCheck.getTransverseMomentumSquared() )
        {
          return acceptRatherThanReject;
        }
        else
        {
          return !acceptRatherThanReject;
        }
      }

    }
    typedef FilterRuleClass::TransverseMomentumCut FilterOnTransverseMomentum;

  }

}

#endif /* TRANSVERSEMOMENTUMCUT_HPP_ */
