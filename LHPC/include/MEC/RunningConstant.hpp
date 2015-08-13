/*
 * RunningConstant.hpp
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

#ifndef RUNNINGCONSTANT_HPP_
#define RUNNINGCONSTANT_HPP_

#include <string>
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  // this is a class to hold the information about the mass of a particle in
  // the FLHA format.
  class RunningConstant
  {
  public:
    RunningConstant();
    RunningConstant( RunningConstant const& copySource );
    ~RunningConstant();

    double
    getValue() const;
    int
    getScheme() const;
    double
    getScale() const;
    void
    setValues( double const valueDouble,
               int const schemeType,
               double const evaluationScale );
    void
    setFromString( std::string const& valuesString );
    std::string
    getAsString() const;


  protected:
    double valueDouble;
    int schemeType;
    double evaluationScale;
  };





  inline double
  RunningConstant::getValue() const
  {
    return valueDouble;
  }

  inline int
  RunningConstant::getScheme() const
  {
    return schemeType;
  }

  inline double
  RunningConstant::getScale() const
  {
    return evaluationScale;
  }

  inline void
  RunningConstant::setValues( double const valueDouble,
                              int const schemeType,
                              double const evaluationScale )
  {
    this->valueDouble = valueDouble;
    this->schemeType = schemeType;
    this->evaluationScale = evaluationScale;
  }

  inline void
  RunningConstant::setFromString( std::string const& valuesString )
  {
    std::string firstRemainder;
    std::string secondRemainder;
    valueDouble
    = BOL::StringParser::stringToDouble( BOL::StringParser::firstWordOf(
                                                                  valuesString,
                                                               &firstRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    schemeType
    = BOL::StringParser::stringToInt( BOL::StringParser::firstWordOf(
                                                                firstRemainder,
                                                              &secondRemainder,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
    evaluationScale = BOL::StringParser::stringToDouble( secondRemainder );
  }

  inline std::string
  RunningConstant::getAsString() const
  {
    std::string returnString( BOL::StringParser::doubleToString( valueDouble,
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

#endif /* RUNNINGCONSTANT_HPP_ */
