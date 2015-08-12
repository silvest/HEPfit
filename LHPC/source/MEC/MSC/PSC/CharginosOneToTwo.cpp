/*
 * CharginosOneToTwo.cpp
 *
 *  Created on: Jan 18, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "MEC.hpp"

namespace LHPC
{
  namespace MassSpectrumClass
  {
    CharginosOneToTwo::CharginosOneToTwo( bool const isVerbose,
                                    std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        positiveCharginoOne( PDGIX::positiveCharginoOne,
                             PDGVII::positiveCharginoOne,
                             mapAndVectorAndBools,
                             false,
                             "C1",
                             "${\\tilde{{\\chi}}}_{1}^{+}$" ),
        negativeCharginoOne( positiveCharginoOne,
                             "C1bar",
                             "${\\tilde{{\\chi}}}_{1}^{-}$" ),
        positiveCharginoTwo( PDGIX::positiveCharginoTwo,
                             PDGVII::positiveCharginoTwo,
                             mapAndVectorAndBools,
                             false,
                             "C2",
                             "${\\tilde{{\\chi}}}_{2}^{+}$" ),
        negativeCharginoTwo( positiveCharginoTwo,
                             "C2bar",
                             "${\\tilde{{\\chi}}}_{2}^{-}$" ),
        positiveCharginoPointers( 2,
                                  &positiveCharginoOne ),
        negativeCharginoPointers( 2,
                                  &negativeCharginoOne )
    {
      positiveCharginoPointers[ 1 ] = &positiveCharginoTwo;
      negativeCharginoPointers[ 1 ] = &negativeCharginoTwo;
    }

    CharginosOneToTwo::~CharginosOneToTwo()
    {
      // does nothing.
    }

  }

}
