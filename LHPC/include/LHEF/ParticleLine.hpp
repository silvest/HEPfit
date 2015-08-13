/*
 * ParticleLine.hpp
 *
 *  Created on: Jan 25, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef PARTICLELINE_HPP_
#define PARTICLELINE_HPP_

#include <string>
#include <vector>
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  namespace LHEF
  {
    /* this is a class to hold the data in a single line of an event in the
     * Les Houches Event File format. instances of this class must be told
     * their place in a VectorlikeArray.
     */
    class ParticleLine
    {
    public:
      static ParticleLine*
      copyParticleLine( ParticleLine const& copySource );
      static ParticleLine&
      sumMomentaAsLine( ParticleLine& lineToBecomeSum,
                        std::list< ParticleLine const* >& linesToSum );
      // this returns a ParticleLine that has NaNs and NULLs for all entries
      // apart from the 4-momentum entries, which are set to the appropriate
      // sums of the entries from the ParticleLines pointed at by the list.
      static double
      transverseMomentumSquaredOfVectorSum(
                                std::list< ParticleLine const* >& linesToSum );
      static double
      transverseMomentumOfVectorSum(
                                std::list< ParticleLine const* >& linesToSum );
      static double
      scalarSumOfTransverseMomentum(
                                std::list< ParticleLine const* >& linesToSum );
      static bool
      isLowToHighTransverseMomentum( ParticleLine const* firstLine,
                                     ParticleLine const* secondLine );
      // this returns true if firstLine has transverse momentum less than or
      // equal to that of secondLine.

      ParticleLine();
      ParticleLine( ParticleLine const& copySource );
      ~ParticleLine();

      int
      getParticleCode() const;
      int
      IDUP() const{ return getParticleCode(); }
      int
      getInitialOrIntermediateOrFinalState() const;
      int
      ISTUP() const{ return getInitialOrIntermediateOrFinalState(); }
      int
      getPrimaryMotherLineNumber() const;
      int
      getSecondaryMotherLineNumber() const;
      int
      MOTHUP( int const whichMother = 1 ) const;
      int
      getColorCode() const;
      int
      getAnticolorCode() const;
      int
      ICOLUP( int const whichColorCode = 1 ) const;
      double
      getXMomentum() const;
      double
      getYMomentum() const;
      double
      getZMomentum() const;
      double
      getEnergy() const;
      double
      getMass() const;
      double
      IPUP( int const whichComponent ) const;
      double
      getProperLifetime() const;
      double
      VTIMTUP() const{ return getProperLifetime(); }
      int
      getTwiceSpin() const;
      int
      SPINUP() const{ return getTwiceSpin(); }
      int
      getOwnLineNumber() const;
      ParticleLine const*
      getPrimaryMotherLinePointer() const;
      ParticleLine const*
      getSecondaryMotherLinePointer() const;
      std::vector< ParticleLine* > const&
      getDaughterLines() const;

      // some basic four-momentum stuff:
      double
      getTransverseMomentum() const;
      double
      getTransverseMomentumSquared() const;
      double
      getThreeMomentumMagnitude() const;
      double
      getPseudorapidity() const;

      // stuff for setting up the ParticleLine:
      bool
      recordLine( int const ownLineNumber,
                  std::string const& lineAsString );
      /* this interprets a string as the data it is meant to correspond to. it
       * returns false if the wrong number of data were given, & in that case
       * it sets all its entries to BOL::UsefulStuff::reallyWrongValue. it
       * returns true if all the entries were filled as expected. this function
       * also clears daughterLines.
       */
      void
      addDaughterLine( ParticleLine* const daughterLinePointer );
      void
      setMotherLines( ParticleLine* const primaryMotherLinePointer,
                      ParticleLine* const secondaryMotherLinePointer );
      // this also calls addDaughterLine on each non-NULL mother line.


    protected:
      // the LHEF data in the order in which they appear in LHE format:
      int particleCode;
      // IDUP, the particle's PDG code.
      int initialOrIntermediateOrFinalState;
      /* ISTUP, "status", which is whether the particle is initial-state (-1),
       * intermediate-state (+2), or final state (+1).
       * (there may be other values, I haven't seen them.)
       */
      int primaryMotherLineNumber;
      /* MOTHUP(1), the number of the line of the particle which produced the
       * particle of this line, if from a decay, or 1 of the particles from an
       * annihilation.
       */
      int secondaryMotherLineNumber;
      // MOTHUP(2), the number of the line of the other particle which produced
      // the particle of this line from an annihilation.
      int colorCode;
      // ICOLUP, the color of the particle if a quark or a gluon or similar (0
      // for a color singlet).
      int anticolorCode;
      /* ICOLUP, the anticolor of the particle if an antiquark or a gluon or
       * similar (0 for a color singlet).
       * I think color & anticolor cover the possibilities, though there might
       * be exotic things that have a different interpretation.
       */
      double xMomentum;
      // PUP(1), the x-component of the particle momentum in GeV.
      double yMomentum;
      // PUP(2), the y-component of the particle momentum in GeV.
      double zMomentum;
      // PUP(3), the z-component of the particle momentum in GeV (the beam is
      // aligned with the z-axis).
      double tMomentum;
      // PUP(4), the energy of the particle in GeV.
      double restMass;
      // PUP(5), the mass of the particle in GeV.
      double properLifetime;
      // VTIMUP, the speed of light divided by the particle's decay width.
      int twiceSpin;
      // SPINUP, 2 * the spin of the particle (negative values indicating
      // left-handed spin) (e.g. a left-handed fermion has -1.0).

      // additional useful data:
      int ownLineNumber;
      // number of the line of the event which this ParticleLine corresponds
      // to (not actually specified by the LHEF format, but handy).
      ParticleLine* primaryMotherLinePointer;
      // pointer to the ParticleLine referred to by primaryMotherLineNumber.
      ParticleLine* secondaryMotherLinePointer;
      // pointer to the ParticleLine referred to by secondaryMotherLineNumber.
      std::vector< ParticleLine* > daughterLines;
      // these are the ParticleLines which refer this ParticleLine with their
      // primaryMotherLineNumbers & secondaryMotherLineNumbers.

      // this is for parsing the data line:
      BOL::VectorlikeArray< std::string > lineParser;
    };





    inline ParticleLine*
    ParticleLine::copyParticleLine( ParticleLine const& copySource )
    {
      return new ParticleLine( copySource );
    }

    inline ParticleLine&
    ParticleLine::sumMomentaAsLine( ParticleLine& lineToBecomeSum,
                                 std::list< ParticleLine const* >& linesToSum )
    {
      lineToBecomeSum.tMomentum = 0.0;
      lineToBecomeSum.xMomentum = 0.0;
      lineToBecomeSum.yMomentum = 0.0;
      lineToBecomeSum.zMomentum = 0.0;
      for( std::list< ParticleLine const* >::iterator
           lineIterator( linesToSum.begin() );
           linesToSum.end() != lineIterator;
           ++lineIterator )
      {
        lineToBecomeSum.tMomentum += (*lineIterator)->getEnergy();
        lineToBecomeSum.xMomentum += (*lineIterator)->getXMomentum();
        lineToBecomeSum.yMomentum += (*lineIterator)->getYMomentum();
        lineToBecomeSum.zMomentum += (*lineIterator)->getZMomentum();
      }
      return lineToBecomeSum;
    }

    inline double
    ParticleLine::transverseMomentumSquaredOfVectorSum(
                                 std::list< ParticleLine const* >& linesToSum )
    {
      double xTotal( 0.0 );
      double yTotal( 0.0 );
      for( std::list< ParticleLine const* >::iterator
           lineIterator( linesToSum.begin() );
           linesToSum.end() != lineIterator;
           ++lineIterator )
      {
        xTotal += (*lineIterator)->getXMomentum();
        yTotal += (*lineIterator)->getYMomentum();
      }
      return ( ( xTotal * xTotal ) + ( yTotal * yTotal ) );
    }

    inline double
    ParticleLine::transverseMomentumOfVectorSum(
                                 std::list< ParticleLine const* >& linesToSum )
    {
      return sqrt( transverseMomentumSquaredOfVectorSum( linesToSum ) );
    }

    inline double
    ParticleLine::scalarSumOfTransverseMomentum(
                                 std::list< ParticleLine const* >& linesToSum )
    {
      double scalarSum( 0.0 );
      for( std::list< ParticleLine const* >::iterator
           lineIterator( linesToSum.begin() );
           linesToSum.end() != lineIterator;
           ++lineIterator )
      {
        scalarSum += (*lineIterator)->getTransverseMomentum();
      }
      return scalarSum;
    }

    inline bool
    ParticleLine::isLowToHighTransverseMomentum( ParticleLine const* firstLine,
                                               ParticleLine const* secondLine )
    // this returns true if firstLine has transverse momentum less than or
    // equal to that of secondLine.
    {
      if( firstLine->getTransverseMomentumSquared()
          <= secondLine->getTransverseMomentumSquared() )
        // actually, the square of the transverse momentum is compared, because
        // it's faster.
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    inline int
    ParticleLine::getParticleCode() const
    {
      return particleCode;
    }

    inline int
    ParticleLine::getInitialOrIntermediateOrFinalState() const
    {
      return initialOrIntermediateOrFinalState;
    }

    inline int
    ParticleLine::getPrimaryMotherLineNumber() const
    {
      return primaryMotherLineNumber;
    }

    inline int
    ParticleLine::getSecondaryMotherLineNumber() const
    {
      return secondaryMotherLineNumber;
    }

    inline int
    ParticleLine::MOTHUP( int const whichMother ) const
    {
      if( 1 == whichMother )
      {
        return getPrimaryMotherLineNumber();
      }
      else if( 2 == whichMother )
      {
        return getSecondaryMotherLineNumber();
      }
      else
      {
        return (int)BOL::UsefulStuff::notANumber;
      }
    }

    inline int
    ParticleLine::getColorCode() const
    {
      return colorCode;
    }

    inline int
    ParticleLine::getAnticolorCode() const
    {
      return anticolorCode;
    }

    inline int
    ParticleLine::ICOLUP( int const whichColorCode ) const
    {
      if( 1 == whichColorCode )
      {
        return getColorCode();
      }
      else if( 2 == whichColorCode )
      {
        return getAnticolorCode();
      }
      else
      {
        return (int)BOL::UsefulStuff::notANumber;
      }
    }

    inline double
    ParticleLine::getXMomentum() const
    {
      return xMomentum;
    }

    inline double
    ParticleLine::getYMomentum() const
    {
      return yMomentum;
    }

    inline double
    ParticleLine::getZMomentum() const
    {
      return zMomentum;
    }

    inline double
    ParticleLine::getEnergy() const
    {
      return tMomentum;
    }

    inline double
    ParticleLine::getMass() const
    {
      return restMass;
    }

    inline double
    ParticleLine::getProperLifetime() const
    {
      return properLifetime;
    }

    inline int
    ParticleLine::getTwiceSpin() const
    {
      return twiceSpin;
    }

    inline int
    ParticleLine::getOwnLineNumber() const
    {
      return ownLineNumber;
    }

    inline ParticleLine const*
    ParticleLine::getPrimaryMotherLinePointer() const
    {
      return primaryMotherLinePointer;
    }

    inline ParticleLine const*
    ParticleLine::getSecondaryMotherLinePointer() const
    {
      return secondaryMotherLinePointer;
    }

    inline std::vector< ParticleLine* > const&
    ParticleLine::getDaughterLines() const
    {
      return daughterLines;
    }

    inline double
    ParticleLine::getTransverseMomentum() const
    {
      return sqrt( getTransverseMomentumSquared() );
    }

    inline double
    ParticleLine::getTransverseMomentumSquared() const
    {
      return ( ( getXMomentum() * getXMomentum() )
               + ( getYMomentum() * getYMomentum() ) );
    }

    inline double
    ParticleLine::getThreeMomentumMagnitude() const
    {
      return sqrt( getTransverseMomentumSquared()
                   + ( getZMomentum() * getZMomentum() ) );
    }

    inline double
    ParticleLine::getPseudorapidity() const
    {
      double threeMomentumMagnitude( getThreeMomentumMagnitude() );
      if( ( threeMomentumMagnitude == getZMomentum() )
          ||
          ( threeMomentumMagnitude == -(getZMomentum()) ) )
      {
        return BOL::UsefulStuff::notANumber;
      }
      else
      {
        return ( 0.5 * log( ( threeMomentumMagnitude + getZMomentum() )
                            / ( threeMomentumMagnitude - getZMomentum() ) ) );
      }
    }

    inline void
    ParticleLine::addDaughterLine( ParticleLine* const daughterLinePointer )
    {
      daughterLines.push_back( daughterLinePointer );
    }

    inline void
    ParticleLine::setMotherLines( ParticleLine* const primaryMotherLinePointer,
                               ParticleLine* const secondaryMotherLinePointer )
    // this also calls addDaughterLine on each non-NULL mother line.
    {
      this->primaryMotherLinePointer = primaryMotherLinePointer;
      if( NULL != primaryMotherLinePointer )
      {
        primaryMotherLinePointer->addDaughterLine( this );
      }
      this->secondaryMotherLinePointer = secondaryMotherLinePointer;
      if( NULL != secondaryMotherLinePointer )
      {
        secondaryMotherLinePointer->addDaughterLine( this );
      }
    }

  }

}

#endif /* PARTICLELINE_HPP_ */
