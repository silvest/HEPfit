/*
 * LhcoEvent.hpp
 *
 *  Created on: Jun 29, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef LHCOEVENT_HPP_
#define LHCOEVENT_HPP_

#include <list>
#include <vector>
#include "ObjectLine.hpp"


namespace LHPC
{
  namespace LHCO
  {
    class LhcoEvent
    {
    public:
      enum objectType
      {
        photonObject = 0,
        electronObject = 1,
        muonObject = 2,
        tauObject = 3,
        jetObject = 4,
        missingEnergyObject = 6
      };
      LhcoEvent( bool const isVerbose );
      LhcoEvent( LhcoEvent const& copySource );
      ~LhcoEvent();

      int
      getNumberOfObjects() const;
      int
      getEventNumber() const;
      int
      getTriggerWord() const;
      ObjectLine const&
      operator[]( int const whichLine ) const;
      std::list< ObjectLine const* > const&
      getPhotons() const;
      std::list< ObjectLine const* > const&
      getElectrons() const;
      std::list< ObjectLine const* > const&
      getMuons() const;
      std::list< ObjectLine const* > const&
      getTaus() const;
      std::list< ObjectLine const* > const&
      getJets() const;
      ObjectLine const*
      getMissingEnergy() const;
      std::list< ObjectLine const* > const&
      getObjectsOfType( int const whichType ) const;
      std::string const&
      getAsString() const;

      // stuff for setting up the LhcoEvent:
      int
      recordLine( std::string const& lineAsString );
      // this parses the 1st word of lineAsString and returns it, first adding
      // a new ObjectLine if appropriate, noting its pointer in the appropriate
      // std::list. if lineAsString is the start of a new event (as the line
      // begins with 0, which is also what this function will return), it
      // prepares nextEventNumber & nextTriggerWord based on this line.
      void
      prepareForNextEvent();
      // this clears objectLines, copies nextEventNumber & nextTriggerWord
      // into eventNumber & triggerWord, & prepares eventAsString.


    protected:
      static std::string const trimmingChars;
      static bool const trueForVerbosity;
      static std::string const eventCommentLine;
      static int const charactersForEventNumber;
      static int const charactersForTriggerWord;

      int eventNumber;
      // the number of the event, given in line 0.
      int triggerWord;
      // the number encoding information about the trigger, given in line 0.
      int nextEventNumber;
      // eventNumber for the next event.
      int nextTriggerWord;
      // triggerWord for the next event.

      BOL::VectorlikeArray< ObjectLine > objectLines;
      int objectLineNumber;

      ObjectLine const* missingEnergyLinePointer;
      std::vector< std::list< ObjectLine const* > > objectLists;

      bool const isVerbose;
      std::string eventAsString;
      ObjectLine const* interpretingObjectLine;
      // this is for parsing the data line:
      BOL::VectorlikeArray< std::string > lineParser;

      std::list< ObjectLine const* >&
      getObjectList( unsigned int const whichType );
      // the type of object:
      // 0: photon
      // 1: electron
      // 2: muon
      // 3: hadronically-decaying tau lepton
      // 4: jet
      // 5: not defined
      // 6: missing transverse energy
      // anything else that appears is also undefined.
    };





    inline int
    LhcoEvent::getNumberOfObjects() const
    {
      return objectLines.getSize();
    }

    inline int
    LhcoEvent::getEventNumber() const
    {
      return eventNumber;
    }

    inline int
    LhcoEvent::getTriggerWord() const
    {
      return triggerWord;
    }

    inline ObjectLine const&
    LhcoEvent::operator[]( int const whichLine ) const
    {
      return objectLines[ ( whichLine - 1 ) ];
    }

    inline std::list< ObjectLine const* > const&
    LhcoEvent::getPhotons() const
    {
      return objectLists[ (unsigned int)photonObject ];
    }

    inline std::list< ObjectLine const* > const&
    LhcoEvent::getElectrons() const
    {
      return objectLists[ (unsigned int)electronObject ];
    }

    inline std::list< ObjectLine const* > const&
    LhcoEvent::getMuons() const
    {
      return objectLists[ (unsigned int)muonObject ];
    }

    inline std::list< ObjectLine const* > const&
    LhcoEvent::getTaus() const
    {
      return objectLists[ (unsigned int)tauObject ];
    }

    inline std::list< ObjectLine const* > const&
    LhcoEvent::getJets() const
    {
      return objectLists[ (unsigned int)jetObject ];
    }

    inline ObjectLine const*
    LhcoEvent::getMissingEnergy() const
    {
      return missingEnergyLinePointer;
    }

    inline std::list< ObjectLine const* > const&
    LhcoEvent::getObjectsOfType( int const whichType ) const
    {
      return objectLists[ whichType ];
    }

    inline std::string const&
    LhcoEvent::getAsString() const
    {
      return eventAsString;
    }

    inline void
    LhcoEvent::prepareForNextEvent()
    // this clears objectLines, copies nextEventNumber & nextTriggerWord
    // into eventNumber & triggerWord, & prepares eventAsString.
    {
      objectLines.clearEntries();
      for( int whichList( objectLists.size() - 1 );
           0 <= whichList;
           --whichList )
      {
        objectLists[ whichList ].clear();
      }
      missingEnergyLinePointer = NULL;
      eventNumber = nextEventNumber;
      triggerWord = nextTriggerWord;
      eventAsString.assign( eventCommentLine );
      eventAsString.append( "0 " );
      eventAsString.append( BOL::StringParser::intToSpacePaddedString(
                                                               nextEventNumber,
                                                  charactersForEventNumber ) );
      eventAsString.append( " " );
      eventAsString.append( BOL::StringParser::intToSpacePaddedString(
                                                               nextTriggerWord,
                                                  charactersForTriggerWord ) );
    }

    inline std::list< ObjectLine const* >&
    LhcoEvent::getObjectList( unsigned int const whichType )
    // the type of object:
    // 0: photon
    // 1: electron
    // 2: muon
    // 3: hadronically-decaying tau lepton
    // 4: jet
    // 5: not defined
    // 6: missing transverse energy
    // anything else that appears is also undefined.
    {
      if( objectLists.size() <= whichType )
      {
        objectLists.resize( whichType + 1 );
      }
      return objectLists[ whichType ];
    }

  }

}

#endif /* LHCOEVENT_HPP_ */
