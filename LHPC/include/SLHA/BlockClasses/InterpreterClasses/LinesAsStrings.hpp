/*
 * LinesAsStrings.hpp
 *
 *  Created on: Mar 12, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef LINESASSTRINGS_HPP_
#define LINESASSTRINGS_HPP_

#include <string>
#include "InterpreterTemplate.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace InterpreterClass
    {
      // this class interprets SLHA blocks as just a set of full lines as
      // std::strings.
      class LinesAsStrings : public InterpreterTemplate< std::string >
      {
      public:
        LinesAsStrings();
        virtual
        ~LinesAsStrings();

        std::string
        operator()( int const whichLine ) const;
        // this just returns the equivalent line from currentStringBlock.
        std::string
        operator[]( int const whichLine ) const
        { return (*this)( whichLine ); }
        virtual std::string const&
        getAsString();
        // see base version's description.
        virtual void
        clearEntries();
        // derived classes should clear their interpreted values.


      protected:
        virtual void
        interpretCurrentStringBlock();
      };





      inline std::string
      LinesAsStrings::operator()( int const whichLine ) const
      // this just returns the equivalent line from currentStringBlock.
      {
        std::string returnString( "" );
        if( NULL != this->currentStringBlock )
        {
          returnString.append(
                            (*(this->currentStringBlock))[ whichLine ].first );
          returnString.append(
                           (*(this->currentStringBlock))[ whichLine ].second );
        }
        return returnString;
      }

      inline std::string const&
      LinesAsStrings::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        int numberOfLinesToPrint( 0 );
        if( NULL != this->currentStringBlock )
        {
          numberOfLinesToPrint
          = this->currentStringBlock->getNumberOfBodyLines();
        }
        for( int whichLine( 1 );
             numberOfLinesToPrint >= whichLine;
             ++whichLine )
        {
          this->stringInterpretation.append(
                            (*(this->currentStringBlock))[ whichLine ].first );
          this->stringInterpretation.append(
                           (*(this->currentStringBlock))[ whichLine ].second );
          this->stringInterpretation.append( "\n" );
        }
        return this->stringInterpretation;
      }

      inline void
      LinesAsStrings::clearEntries()
      // this ensures that the entry at soughtIndex exists, filling out with
      // copies of defaultUnsetValue, & returns it.
      {
        // does nothing.
      }

      inline void
      LinesAsStrings::interpretCurrentStringBlock()
      {
        // does nothing.
      }

    }

  }

}


#endif /* LINESASSTRINGS_HPP_ */
