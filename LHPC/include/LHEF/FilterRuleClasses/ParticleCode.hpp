/*
 * ParticleCode.hpp
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

#ifndef PARTICLECODE_HPP_
#define PARTICLECODE_HPP_

#include <vector>
#include "../FilterRule.hpp"

namespace LHPC
{
  namespace LHEF
  {
    namespace FilterRuleClass
    {
      // this class returns acceptRatherThanReject if it finds that
      // lineToCheck.getParticleCode() matches an int in acceptableValues.
      class ParticleCode : public FilterRule
      {
      public:
        ParticleCode( std::vector< int > const& acceptableValues,
                      bool const acceptRatherThanReject = true );
        ParticleCode( int const acceptableValue,
                      bool const acceptRatherThanReject = true );
        virtual
        ~ParticleCode();

        virtual bool
        operator()( ParticleLine const& lineToCheck ) const;


      protected:
        std::vector< int > acceptableValues;
      };



      inline bool
      ParticleCode::operator()( ParticleLine const& lineToCheck ) const
      {
        bool returnValue( !acceptRatherThanReject );
        for( int valueIndex( acceptableValues.size() - 1 );
             0 <= valueIndex;
             --valueIndex )
        {
          if( lineToCheck.getParticleCode() == acceptableValues[ valueIndex ] )
          {
            returnValue = acceptRatherThanReject;
            valueIndex = -1;
          }
        }
        return returnValue;
      }

    }
    typedef FilterRuleClass::ParticleCode FilterOnParticleCode;

  }

}

#endif /* PARTICLECODE_HPP_ */
