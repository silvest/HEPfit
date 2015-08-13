/*
 * InterfaceToClhepLorentzVectorClass.hpp
 *
 *  Created on: Jan 27, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

// including this file gives a way of filling CLHEP::HepLorentzVector instances
// with data from ParticleLine instances.

#ifndef INTERFACETOCLHEPLORENTZVECTORCLASS_HPP_
#define INTERFACETOCLHEPLORENTZVECTORCLASS_HPP_

#ifdef HEP_LORENTZVECTOR_H
#include "ParticleLine.hpp"

namespace LHPC
{
  namespace LHEF
  {
    class InterfaceToClhepLorentzVector
    {
    public:
      static void
      fillFromLine( CLHEP::HepLorentzVector const& vectorToFill,
                    ParticleLine const& lineToConvert );
    };
    typedef InterfaceToClhepLorentzVector LineToVec;




    inline void
    InterfaceToClhepLorentzVector::fillVectorFromParticleLine(
                                   CLHEP::HepLorentzVector const& vectorToFill,
                                            ParticleLine const& lineToConvert )
    {
      vectorToFill.set( lineToConvert.getXMomentum(),
                        lineToConvert.getYMomentum(),
                        lineToConvert.getZMomentum(),
                        lineToConvert.getEnergy() );
    }

  }

}  // end of LHEF namespace

#endif /* HEP_LORENTZVECTOR_H */

#endif /* INTERFACETOCLHEPLORENTZVECTORCLASS_HPP_ */
