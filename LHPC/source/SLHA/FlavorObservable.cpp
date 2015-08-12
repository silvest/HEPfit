/*
 * FlavorObservable.cpp
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
  int const FlavorObservable::spacesBetweenCodes( 3 );
  int const FlavorObservable::minimumDigitsForCodes( 9 );

  FlavorObservable::FlavorObservable() :
      valueDouble( BOL::UsefulStuff::notANumber ),
      evaluationScale( BOL::UsefulStuff::notANumber ),
      daughterParticleCodes()
  {
    // just an initialization list.
  }

  FlavorObservable::FlavorObservable( FlavorObservable const& copySource ) :
      valueDouble( copySource.valueDouble ),
      evaluationScale( copySource.evaluationScale ),
      daughterParticleCodes( copySource.daughterParticleCodes )
  {
    // just an initialization list.
  }

  FlavorObservable::~FlavorObservable()
  {
    // does nothing.
  }

  void
  FlavorObservable::setFromString( std::string const& valuesString )
  {
    std::string firstRemainder;
    std::string secondRemainder;
    valueDouble
    = BOL::StringParser::stringToDouble( BOL::StringParser::firstWordOf(
                                                                  valuesString,
                                                               &firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    evaluationScale
    = BOL::StringParser::stringToDouble( BOL::StringParser::firstWordOf(
                                                                firstRemainder,
                                                              &secondRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    int numberOfDaughterParticles( BOL::StringParser::stringToInt(
                                                BOL::StringParser::firstWordOf(
                                                               secondRemainder,
                                                               &firstRemainder,
                            BOL::StringParser::whitespaceAndNewlineChars ) ) );
    secondRemainder.assign( BOL::StringParser::trimFromFrontAndBack(
                                                                firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    daughterParticleCodes.clear();
    while( !(secondRemainder.empty()) )
    {
      daughterParticleCodes.push_back( BOL::StringParser::stringToInt(
                                                BOL::StringParser::firstWordOf(
                                                               secondRemainder,
                                                               &firstRemainder,
                            BOL::StringParser::whitespaceAndNewlineChars ) ) );
      secondRemainder.assign( BOL::StringParser::trimFromFrontAndBack(
                                                                firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    }
    if( ( 0 < numberOfDaughterParticles )
        &&
        ( (size_t)numberOfDaughterParticles != daughterParticleCodes.size() ) )
    {
      std::cout
      << std::endl
      << "LHPC::warning! An FOBS line declared a different number of daughter"
      << " particles ( " << numberOfDaughterParticles << " ) to the actual"
      << " number of daughter particle codes it had ( "
      << daughterParticleCodes.size()
      << " )! The declared number is being ignored in favor of the number of"
      << " codes read in.";
      std::cout
      << std::endl
      << "input string: \"" << valuesString << "\"";
      std::cout << std::endl;
    }
  }

  std::string
  FlavorObservable::getAsString() const
  {
    std::string returnString( BOL::StringParser::doubleToString( valueDouble,
                                                                 9,
                                                                 3 ) );
    returnString.append( spacesBetweenCodes,
                         ' ' );
    returnString.append( BOL::StringParser::doubleToString( evaluationScale,
                                                            9,
                                                            3 ) );
    for( std::list< int >::const_iterator
         daughterIterator( daughterParticleCodes.begin() );
         daughterParticleCodes.end() != daughterIterator;
         ++daughterIterator )
    {
      returnString.append( spacesBetweenCodes,
                           ' ' );
      returnString.append( BOL::StringParser::intToSpacePaddedString(
                                                             *daughterIterator,
                                                         minimumDigitsForCodes,
                                                                      "" ) );
    }
    return returnString;
  }

}
