/*
 * BlockInterpreter.hpp
 *
 *  Created on: Mar 11, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef BLOCKINTERPRETER_HPP_
#define BLOCKINTERPRETER_HPP_

#include "BOLlib/include/BOLlib.hpp"
#include "../BaseStringBlock.hpp"

namespace LHPC
{
  namespace SLHA
  {
    // this abstract base class provides a base class for interpreting
    // BaseStringBlock instances.
    class BlockInterpreter
    {
    public:
      static BOL::StringParser const slhaDoubleMaker;
      static BOL::StringParser const slhaIntHelper;
      static BOL::StringParser const particleCodeMaker;
      static bool const defaultVerbosity;

      BlockInterpreter();
      virtual
      ~BlockInterpreter();

      void
      interpretStringBlock(
                       BlockClass::BaseStringBlock const& stringsToInterpret );
      std::string
      getLineWithoutComment( int const whichLine ) const;
      std::string
      getLineWithComment( int const whichLine ) const;
      /* information in comments can be decoded by someone else's code, because
       * putting important information in comments is an entirely unacceptably
       * stupid idea, in my humble opinion.
       */
      double
      getScale() const;
      virtual std::string const&
      getAsString() = 0;
      // derived classes should return their block as a single string of
      // re-formatted interpreted values.
      virtual void
      clearEntries() = 0;
      // derived classes should clear their interpreted values.


    protected:
      BlockClass::BaseStringBlock const* currentStringBlock;
      std::string stringInterpretation;

      virtual void
      interpretCurrentStringBlock() = 0;
    };





    inline void
    BlockInterpreter::interpretStringBlock(
                        BlockClass::BaseStringBlock const& stringsToInterpret )
    {
      clearEntries();
      currentStringBlock = &stringsToInterpret;
      interpretCurrentStringBlock();
    }

    inline std::string
    BlockInterpreter::getLineWithoutComment( int const whichLine ) const
    {
      if( NULL == currentStringBlock )
      {
        std::string lineReturnString;
        lineReturnString.assign( "LHPC::error! this block has not been" );
        lineReturnString.append( " registered with an SlhaParser properly!" );
        return lineReturnString;
      }
      else
      {
        return (*currentStringBlock)[ whichLine ].first;
      }
    }

    inline std::string
    BlockInterpreter::getLineWithComment( int const whichLine ) const
    /* information in comments can be decoded by someone else's code, because
     * putting important information in comments is an entirely unacceptably
     * stupid idea, in my humble opinion.
     */
    {
      std::string lineReturnString;
      if( NULL == currentStringBlock )
      {
        lineReturnString.assign( "LHPC::error! this block has not been" );
        lineReturnString.append( " registered with an SlhaParser properly!" );
      }
      else
      {
        lineReturnString.assign( (*currentStringBlock)[ whichLine ].first );
        lineReturnString.append( (*currentStringBlock)[ whichLine ].second );
      }
      return lineReturnString;
    }

    inline double
    BlockInterpreter::getScale() const
    {
      if( NULL != currentStringBlock )
      {
        return currentStringBlock->getScale();
      }
      else
      {
        return 0.0;
      }
    }

  }

}

#endif /* BLOCKINTERPRETER_HPP_ */
