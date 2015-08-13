/*
 * RunningConstantError.hpp
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

#ifndef RUNNINGCONSTANTERROR_HPP_
#define RUNNINGCONSTANTERROR_HPP_

#include <string>
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  // this is a class to hold the information about the mass of a particle in
  // the FLHA format.
  class RunningConstantError
  {
  public:
    RunningConstantError();
    RunningConstantError( RunningConstantError const& copySource );
    ~RunningConstantError();

    double
    getMinusUncertainty() const;
    double
    getPlusUncertainty() const;
    int
    getScheme() const;
    double
    getScale() const;
    void
    setValues( double const minusUncertainty,
               double const plusUncertainty,
               int const schemeType,
               double const evaluationScale );
    void
    setFromString( std::string const& valuesString );
    std::string
    getAsString() const;


  protected:
    double minusUncertainty;
    double plusUncertainty;
    int schemeType;
    double evaluationScale;
  };





  inline double
  RunningConstantError::getMinusUncertainty() const
  {
    return minusUncertainty;
  }

  inline double
  RunningConstantError::getPlusUncertainty() const
  {
    return plusUncertainty;
  }

  inline int
  RunningConstantError::getScheme() const
  {
    return schemeType;
  }

  inline double
  RunningConstantError::getScale() const
  {
    return evaluationScale;
  }

  inline void
  RunningConstantError::setValues( double const minusUncertainty,
                                   double const plusUncertainty,
                                   int const schemeType,
                                   double const evaluationScale )
  {
    this->minusUncertainty = minusUncertainty;
    this->plusUncertainty = plusUncertainty;
    this->schemeType = schemeType;
    this->evaluationScale = evaluationScale;
  }

  inline void
  RunningConstantError::setFromString( std::string const& valuesString )
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
    schemeType
    = BOL::StringParser::stringToInt( BOL::StringParser::firstWordOf(
                                                               secondRemainder,
                                                               &firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    evaluationScale = BOL::StringParser::stringToDouble( firstRemainder );
  }

  inline std::string
  RunningConstantError::getAsString() const
  {
    std::string
    returnString( BOL::StringParser::doubleToString( minusUncertainty,
                                                     9,
                                                     3 ) );
    returnString.append( "   " );
    returnString.append( BOL::StringParser::doubleToString( plusUncertainty,
                                                            9,
                                                            3 ) );
    returnString.append( "   " );
    returnString.append( BOL::StringParser::intToString( schemeType,
                                                         2,
                                                         "",
                                                         "-",
                                                         ' ' ) );
    returnString.append( "   " );
    returnString.append( BOL::StringParser::doubleToString( evaluationScale,
                                                            9,
                                                            3 ) );
    return returnString;
  }

}

#endif /* RUNNINGCONSTANTERROR_HPP_ */
