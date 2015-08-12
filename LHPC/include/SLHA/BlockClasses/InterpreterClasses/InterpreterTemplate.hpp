/*
 * InterpreterTemplate.hpp
 *
 *  Created on: Feb 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef INTERPRETERTEMPLATE_HPP_
#define INTERPRETERTEMPLATE_HPP_

#include "BOLlib/include/BOLlib.hpp"
#include "BlockInterpreter.hpp"
#include "../../../MEC/ExtendedMass.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace InterpreterClass
    {
      // this template class derives from SlhaBlock to provide a base class for
      // blocks with values which are indexed in various ways.
      template< class ValueType >
      class InterpreterTemplate : public BlockInterpreter
      {
      public:
        InterpreterTemplate();
        virtual
        ~InterpreterTemplate();

        virtual void
        setDefaultUnsetValue( ValueType const& defaultUnsetValue );
        virtual void
        setVerbosity( bool const isVerbose );


      protected:
        ValueType defaultUnsetValue;
        bool isVerbose;
        ValueType valueFromString;
        std::string currentWord;
        std::string lineRemainderA;
        std::string lineRemainderB;
        std::string stringFromValue;
        std::string valuePrintingString;

        ValueType const&
        stringToValue( std::string const& stringToConvert );
        // this sets valueFromString according to the interpretation of
        // stringToConvert.
        std::string const&
        valueToString( ValueType const& valueToConvert );
        // this sets stringFromValue according to the interpretation of
        // valueToConvert.
        std::string const&
        valueToPrintingString( ValueType const& valueToPrint );
        // this puts 3 spaces into valuePrintingString, then
        // valueToString( valueToPrint ).
      };





      template< class ValueType >
      inline
      InterpreterTemplate< ValueType >::InterpreterTemplate() :
          BlockInterpreter(),
          defaultUnsetValue(),
          isVerbose( BlockInterpreter::defaultVerbosity ),
          valueFromString(),
          currentWord( "" ),
          lineRemainderA( "" ),
          lineRemainderB( "" ),
          stringFromValue( "no_string_interpretation_given" ),
          valuePrintingString( "   no_string_interpretation_given" )
      {
        // just an initialization list.
      }

      template< class ValueType >
      inline
      InterpreterTemplate< ValueType >::~InterpreterTemplate()
      {
        // does nothing.
      }

      template< class ValueType >
      inline void
      InterpreterTemplate< ValueType >::setDefaultUnsetValue(
                                           ValueType const& defaultUnsetValue )
      {
        this->defaultUnsetValue = defaultUnsetValue;
      }

      template< class ValueType >
      inline void
      InterpreterTemplate< ValueType >::setVerbosity( bool const isVerbose )
      {
        this->isVerbose = isVerbose;
      }

      template< class ValueType >
      inline ValueType const&
      InterpreterTemplate< ValueType >::stringToValue(
                                           std::string const& stringToConvert )
      /* this sets valueFromString according to the interpretation of
       * stringToConvert. this default version leaves stringFromValue as
       * defaultUnsetValue.
       */
      {
        //return defaultUnsetValue;
        valueFromString.setFromString( stringToConvert );
        return valueFromString;
      }

      template<>
      inline int const&
      InterpreterTemplate< int >::stringToValue(
                                           std::string const& stringToConvert )
      // this sets valueFromString according to the interpretation of
      // stringToConvert.
      {
        valueFromString = BOL::StringParser::stringToInt( stringToConvert );
        return valueFromString;
      }

      template<>
      inline double const&
      InterpreterTemplate< double >::stringToValue(
                                           std::string const& stringToConvert )
      // this sets valueFromString according to the interpretation of
      // stringToConvert.
      {
        valueFromString = BOL::StringParser::stringToDouble( stringToConvert );
        return valueFromString;
      }

      template<>
      inline std::pair< double, double > const&
      InterpreterTemplate< std::pair< double, double > >::stringToValue(
                                           std::string const& stringToConvert )
      // this sets valueFromString according to the interpretation of
      // stringToConvert.
      {
        this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                               stringToConvert,
                                                       &(this->lineRemainderA),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
        valueFromString.first
        = BOL::StringParser::stringToDouble( this->currentWord );
        valueFromString.second
        = BOL::StringParser::stringToDouble( this->lineRemainderA );
        return valueFromString;
      }

      template<>
      inline std::string const&
      InterpreterTemplate< std::string >::stringToValue(
                                           std::string const& stringToConvert )
      // this sets valueFromString according to the interpretation of
      // stringToConvert.
      {
        valueFromString.assign( stringToConvert );
        return valueFromString;
      }

      /*template<>
      inline ExtendedMass const&
      InterpreterTemplate< ExtendedMass >::stringToValue(
                                           std::string const& stringToConvert )
      // this sets valueFromString according to the interpretation of
      // stringToConvert.
      {
        this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                               stringToConvert,
                                                       &(this->lineRemainderA),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
        double
        parsedDouble( BOL::StringParser::stringToDouble( currentWord ) );
        this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                          this->lineRemainderA,
                                                       &(this->lineRemainderB),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
        valueFromString.setValues( parsedDouble,
                           BOL::StringParser::stringToInt( this->currentWord ),
                   BOL::StringParser::stringToDouble( this->lineRemainderB ) );
        return valueFromString;
      }*/

      template< class ValueType >
      inline std::string const&
      InterpreterTemplate< ValueType >::valueToString(
                                              ValueType const& valueToConvert )
      /* this sets stringFromValue according to the interpretation of
       * valueToConvert. this default version leaves stringFromValue as a
       * default error string.
       */
      {
        stringFromValue.assign( valueToConvert.getAsString() );
        return stringFromValue;
      }

      template<>
      inline std::string const&
      InterpreterTemplate< int >::valueToString( int const& valueToConvert )
      // this sets stringFromValue according to the interpretation of
      // valueToConvert.
      {
        stringFromValue.assign( BOL::StringParser::intToString( valueToConvert,
                                                                1,
                                                                "",
                                                                "-",
                                                                ' ' ) );
        return stringFromValue;
      }

      template<>
      inline std::string const&
      InterpreterTemplate< double >::valueToString(
                                                 double const& valueToConvert )
      // this sets stringFromValue according to the interpretation of
      // valueToConvert.
      {
        stringFromValue.assign( slhaDoubleMaker.doubleToString(
                                                            valueToConvert ) );
        return stringFromValue;
      }

      template<>
      inline std::string const&
      InterpreterTemplate< std::pair< double, double > >::valueToString(
                            std::pair< double, double > const& valueToConvert )
      // this sets stringFromValue according to the interpretation of
      // valueToConvert.
      {
        stringFromValue.assign( slhaDoubleMaker.doubleToString(
                                                      valueToConvert.first ) );
        stringFromValue.append( "   " );
        stringFromValue.append( slhaDoubleMaker.doubleToString(
                                                     valueToConvert.second ) );
        return stringFromValue;
      }

      template<>
      inline std::string const&
      InterpreterTemplate< std::string >::valueToString(
                                            std::string const& valueToConvert )
      // this sets stringFromValue according to the interpretation of
      // valueToConvert.
      {
        stringFromValue.assign( valueToConvert );
        return stringFromValue;
      }

      /*template<>
      inline std::string const&
      InterpreterTemplate< ExtendedMass >::valueToString(
                                           ExtendedMass const& valueToConvert )
      // this sets stringFromValue according to the interpretation of
      // valueToConvert.
      {
        stringFromValue.assign( slhaDoubleMaker.doubleToString(
                                                  valueToConvert.getValue() ) );
        stringFromValue.append( "   " );
        stringFromValue.append( BOL::StringParser::intToString(
                                                    valueToConvert.getScheme(),
                                                                2,
                                                                "",
                                                                "-",
                                                                ' ' ) );
        stringFromValue.append( "   " );
        stringFromValue.append( slhaDoubleMaker.doubleToString(
                                                 valueToConvert.getScale() ) );
        return stringFromValue;
      }*/

      template< class ValueType >
      inline std::string const&
      InterpreterTemplate< ValueType >::valueToPrintingString(
                                                ValueType const& valueToPrint )
      // this puts 3 spaces into returnString, then
      // valueToString( valueToPrint ).
      {
        valuePrintingString.assign( "   " );
        valuePrintingString.append( this->valueToString( valueToPrint ) );
        return valuePrintingString;
      }

    }

  }

}

#endif /* INTERPRETERTEMPLATE_HPP_ */
