/*
 * FilterRule.hpp
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

#ifndef FILTERRULE_HPP_
#define FILTERRULE_HPP_

#include "ParticleLine.hpp"

namespace LHPC
{
  namespace LHEF
  {
    // this is an abstract base class for filtering ParticleLines based on
    // whether their entries are acceptable.
    class FilterRule
    {
    public:
      FilterRule( bool const acceptRatherThanReject );
      virtual
      ~FilterRule();

      virtual bool
      operator()( ParticleLine const& lineToCheck ) const = 0;

    protected:
      bool const acceptRatherThanReject;
    };

  }

}

#endif /* FILTERRULE_HPP_ */
