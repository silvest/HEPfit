/*
 * PseudorapidityCut.hpp
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

#ifndef PSEUDORAPIDITYCUT_HPP_
#define PSEUDORAPIDITYCUT_HPP_

#include "../FilterRule.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace FilterRuleClass
    {
      // this class returns acceptRatherThanReject if it finds that
      // lineToCheck.getPseudorapidity() is outside the given range.
      class PseudorapidityCut : public FilterRule
      {
      public:
        PseudorapidityCut( double const minimumPseudorapidity,
                           double const maximumPseudorapidity,
                           bool const acceptRatherThanReject = true );
        PseudorapidityCut( double const absoluteMaximumPseudorapidity,
                           bool const acceptRatherThanReject = true );
        virtual
        ~PseudorapidityCut();

        virtual bool
        operator()( ParticleLine const& lineToCheck ) const;


      protected:
        double minimumPseudorapidity;
        double maximumPseudorapidity;
      };



      inline bool
      PseudorapidityCut::operator()( ParticleLine const& lineToCheck ) const
      {
        double linePseudorapidity( lineToCheck.getPseudorapidity() );
        if( ( minimumPseudorapidity <= linePseudorapidity )
            &&
            ( maximumPseudorapidity >= linePseudorapidity ) )
        {
          return acceptRatherThanReject;
        }
        else
        {
          return !acceptRatherThanReject;
        }
      }

    }
    typedef FilterRuleClass::PseudorapidityCut FilterOnPseudorapidity;

  }

}

#endif /* PSEUDORAPIDITYCUT_HPP_ */
