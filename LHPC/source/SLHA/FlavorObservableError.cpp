/*
 * FlavorObservableError.cpp
 *
 *  Created on: Apr 1, 2012 (really!)
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "SLHA.hpp"

namespace LHPC
{
  FlavorObservableError::FlavorObservableError() :
      minusUncertainty( BOL::UsefulStuff::notANumber ),
      plusUncertainty( BOL::UsefulStuff::notANumber ),
      evaluationScale( BOL::UsefulStuff::notANumber ),
      daughterParticleCodes()
  {
    // just an initialization list.
  }

  FlavorObservableError::FlavorObservableError(
                                    FlavorObservableError const& copySource ) :
      minusUncertainty( copySource.minusUncertainty ),
      plusUncertainty( copySource.plusUncertainty ),
      evaluationScale( copySource.evaluationScale ),
      daughterParticleCodes( copySource.daughterParticleCodes )
  {
    // just an initialization list.
  }

  FlavorObservableError::~FlavorObservableError()
  {
    // does nothing.
  }


  void
  FlavorObservableError::setFromString( std::string const& valuesString )
  {
    std::string firstRemainder;
    std::string secondRemainder;
    minusUncertainty
    = BOL::StringParser::stringToDouble( BOL::StringParser::firstWordOf(
                                                                  valuesString,
                                                               &firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    plusUncertainty
    = BOL::StringParser::stringToDouble( BOL::StringParser::firstWordOf(
                                                                firstRemainder,
                                                              &secondRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    evaluationScale
    = BOL::StringParser::stringToDouble( BOL::StringParser::firstWordOf(
                                                               secondRemainder,
                                                               &firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    int numberOfDaughterParticles( BOL::StringParser::stringToInt(
                                                BOL::StringParser::firstWordOf(
                                                                firstRemainder,
                                                              &secondRemainder,
                            BOL::StringParser::whitespaceAndNewlineChars ) ) );
    firstRemainder.assign( BOL::StringParser::trimFromFrontAndBack(
                                                               secondRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    daughterParticleCodes.clear();
    while( !(firstRemainder.empty()) )
    {
      daughterParticleCodes.push_back( BOL::StringParser::stringToInt(
                                                BOL::StringParser::firstWordOf(
                                                                firstRemainder,
                                                              &secondRemainder,
                            BOL::StringParser::whitespaceAndNewlineChars ) ) );
      firstRemainder.assign( BOL::StringParser::trimFromFrontAndBack(
                                                               secondRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    }
    if( ( 0 < numberOfDaughterParticles )
        &&
        ( (size_t)numberOfDaughterParticles != daughterParticleCodes.size() ) )
    {
      std::cout
      << std::endl
      << "LHPC::warning! An FOBSERR line declared a different number of"
      << " daughter particles ( " << numberOfDaughterParticles << " ) to the"
      << " actual number of daughter particle codes it had ( "
      << daughterParticleCodes.size()
      << " )! The declared number is being ignored in favor of the number of"
      << " codes read in.";
      std::cout << std::endl;
    }
  }

  std::string
  FlavorObservableError::getAsString() const
  {
    std::string
    returnString( BOL::StringParser::doubleToString( minusUncertainty,
                                                     9,
                                                     3 ) );
    returnString.append( FlavorObservable::spacesBetweenCodes,
                         ' ' );
    returnString.append( BOL::StringParser::doubleToString( plusUncertainty,
                                                            9,
                                                            3 ) );
    returnString.append( FlavorObservable::spacesBetweenCodes,
                         ' ' );
    returnString.append( BOL::StringParser::doubleToString( evaluationScale,
                                                            9,
                                                            3 ) );
    for( std::list< int >::const_iterator
         daughterIterator( daughterParticleCodes.begin() );
         daughterParticleCodes.end() != daughterIterator;
         ++daughterIterator )
    {
      returnString.append( FlavorObservable::spacesBetweenCodes,
                           ' ' );
      returnString.append( BOL::StringParser::intToSpacePaddedString(
                                                             *daughterIterator,
                                       FlavorObservable::minimumDigitsForCodes,
                                                                      "" ) );
    }
    return returnString;
  }

}
