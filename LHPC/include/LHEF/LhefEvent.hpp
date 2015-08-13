/*
 * LhefEvent.hpp
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

#ifndef LHEFEVENT_HPP_
#define LHEFEVENT_HPP_

#include <iostream>
#include "ParticleLine.hpp"
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  namespace LHEF
  {
    // this class is for holding several ParticleLine instances together with
    // some ints and doubles to describe a complete LHEF-format event.
    class LhefEvent
    {
    public:
      LhefEvent( bool const isVerbose );
      LhefEvent( LhefEvent const& copySource );
      ~LhefEvent();

      int
      getNumberOfParticles() const;
      int
      NUP() const{ return getNumberOfParticles(); }
      int
      getEventId() const;
      int
      IDPRUP() const{ return getEventId(); }
      double
      getEventWeight() const;
      double
      XWGTUP() const{ return getEventWeight(); }
      double
      getEventScale() const;
      double
      SCALUP() const{ return getEventScale(); }
      double
      getAlphaQed() const;
      double
      AQEDUP() const{ return getAlphaQed(); }
      double
      getAlphaQcd() const;
      double
      AQCDUP() const{ return getAlphaQcd(); }
      ParticleLine const&
      getLine( int const whichLine ) const;
      ParticleLine const&
      operator[]( int const whichLine ) const{ return getLine( whichLine ); }
      int
      getEventNumberInFile() const;
      std::string const&
      getAsString() const;
      std::string const
      getAsStringWithTags() const;
      std::string const&
      getOptionalInformation() const;

      // stuff for setting up the LhefEvent:
      bool
      recordEvent( std::string const& eventAsString );
      /* this interprets a string as the data it is meant to correspond to, &
       * ensures that the ParticleLine instances correctly point to each other.
       * if any of the lines had the wrong number of data entries, or if the
       * header number of particles didn't correspond to the number of
       * particle lines, setAsInvalid() is called & then false is returned.
       * true is returned if all the entries were filled as expected.
       */


    protected:
      static std::string const trimmingChars;
      static bool const trueForVerbosity;

      // in the order in which they appear in LHE format:
      // 1st line of the event:
      int numberOfParticles;
      // NUP, the number of particles listed in the event.
      int eventId;
      // IDPRUP, which is not actually defined in hep-ph/0609017...
      double eventWeight;
      // XWGTUP, the weight of the event.
      double eventScale;
      // SCALUP, the energy scale of the event, in GeV.
      double alphaQed;
      // AQEDUP, the value of the electromagnetic alpha at eventScale.
      double alphaQcd;
      // AQCDUP, the value of the color force alpha at eventScale.

      // then numberOfParticles lines, each of which are stored in a
      // ParticleLine:
      BOL::VectorlikeArray< ParticleLine > particleLines;

      int eventNumberInFile;
      // this isn't specified as part of the LHEF format, but it's probably
      // useful.

      // stuff for setting up the LhefEvent:
      BOL::VectorlikeArray< std::string > eventAsLines;
      BOL::VectorlikeArray< std::string > lineAsNumbersAsStrings;
      bool recordingSucceeded;
      std::string headerLine;
      int motherLineNumber;
      ParticleLine* primaryMotherLinePointer;
      ParticleLine* secondaryMotherLinePointer;
      bool const isVerbose;
      std::string eventAsString;
      std::string optionalInformation;

      void
      setAsInvalid();
      // this sets every entry in the header line to
      // BOL::UsefulStuff::notANumber & clears particleLines.
      void
      setUpPointersBetweenParticleLines();
      // this sets up the mother & daughter pointers for each ParticleLine in
      // particleLines.
    };





    inline int
    LhefEvent::getNumberOfParticles() const
    {
      return numberOfParticles;
    }

    inline int
    LhefEvent::getEventId() const
    {
      return eventId;
    }

    inline double
    LhefEvent::getEventWeight() const
    {
      return eventId;
    }

    inline double
    LhefEvent::getEventScale() const
    {
      return eventScale;
    }

    inline double
    LhefEvent::getAlphaQed() const
    {
      return alphaQed;
    }

    inline double
    LhefEvent::getAlphaQcd() const
    {
      return alphaQcd;
    }

    inline ParticleLine const&
    LhefEvent::getLine( int const whichLine ) const
    {
      return particleLines[ whichLine - 1 ];
    }

    inline int
    LhefEvent::getEventNumberInFile() const
    {
      return eventNumberInFile;
    }

    inline std::string const&
    LhefEvent::getAsString() const
    {
      return eventAsString;
    }

    inline std::string const
    LhefEvent::getAsStringWithTags() const
    {
      std::string returnString( "<event>\n" );
      returnString.append( eventAsString );
      returnString.append( "\n</event>" );
      return returnString;
    }

    inline std::string const&
    LhefEvent::getOptionalInformation() const
    {
      return optionalInformation;
    }

    inline void
    LhefEvent::setAsInvalid()
    // this sets every entry in the header line to
    // BOL::UsefulStuff::notANumber & clears particleLines.
    {
      particleLines.clearEntries();
      numberOfParticles = (int)BOL::UsefulStuff::notANumber;
      eventId = (int)BOL::UsefulStuff::notANumber;
      eventWeight = BOL::UsefulStuff::notANumber;
      eventScale = BOL::UsefulStuff::notANumber;
      alphaQed = BOL::UsefulStuff::notANumber;
      alphaQcd = BOL::UsefulStuff::notANumber;

      std::cout
      << std::endl
      << "LHPC::warning! the following string could not be parsed as a valid"
      << " LHEF event:";
      std::cout << std::endl << eventAsString;
      std::cout << std::endl;
    }

  }

}

#endif /* LHEFEVENT_HPP_ */
