/*
 * DenseSinglyIndexed.hpp
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

#ifndef DENSESINGLYINDEXED_HPP_
#define DENSESINGLYINDEXED_HPP_

#include "IndexedInterpreter.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace InterpreterClass
    {
      // this template class interprets SLHA blocks that have a single int
      // index with a single ValueClass value.
      template< class ValueClass >
      class DenseSinglyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        DenseSinglyIndexed();
        virtual
        ~DenseSinglyIndexed();

        ValueClass&
        operator()( int const soughtIndex );
        /* this returns the ValueClass indexed by soughtIndex.. if there is no
         * element at soughtIndex, valueVector is extended with copies of
         * defaultUnsetValue until there is an element at soughtIndex.
         */
        ValueClass const&
        operator()( int const soughtIndex ) const;
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at soughtIndex if there isn't an
         * entry there already.
         */
        ValueClass&
        operator[]( int const soughtIndex )
        { return (*this)( soughtIndex ); }
        ValueClass const&
        operator[]( int const soughtIndex ) const
        { return (*this)( soughtIndex ); }
        bool
        hasEntry( int const soughtIndex ) const;
        // this returns true if there is an entry at soughtIndex.
        virtual std::string const&
        getAsString();
        // see base version's description.
        virtual void
        clearEntries();
        // derived classes should clear their interpreted values.


      protected:
        std::vector< ValueClass > valueVector;
        size_t recordingIndex;

        ValueClass&
        findOrMakeEntry( int const soughtIndex );
        // this ensures that the entry at soughtIndex exists, filling out with
        // copies of defaultUnsetValue, & returns it.
        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      DenseSinglyIndexed< ValueClass >::DenseSinglyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueVector(),
          recordingIndex( 0 )
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      DenseSinglyIndexed< ValueClass >::~DenseSinglyIndexed()
      {
        // does nothing.
      }


      template< class ValueClass >
      inline ValueClass&
      DenseSinglyIndexed< ValueClass >::operator()( int const soughtIndex )
      /* this returns the ValueClass indexed by soughtIndex. if there is no
       * element at soughtIndex, valueVector is extended with copies of
       * defaultUnsetValue until there is an element at soughtIndex.
       */
      {
        return findOrMakeEntry( soughtIndex );
      }

      template< class ValueClass >
      inline ValueClass const&
      DenseSinglyIndexed< ValueClass >::operator()(
                                                  int const soughtIndex ) const
      /* const version of above, though it returns defaultUnsetValue rather
       * than copying in a new element at soughtIndex if there isn't an
       * entry there already.
       */
      {
        if( hasEntry( soughtIndex ) )
        {
          return valueVector[ soughtIndex - 1 ];
        }
        else
        {
          return this->defaultUnsetValue;
        }
      }

      template< class ValueClass >
      inline bool
      DenseSinglyIndexed< ValueClass >::hasEntry( int const soughtIndex ) const
      {
        if( ( 0 < soughtIndex )
            &&
            ( (size_t)soughtIndex <= valueVector.size() ) )
        {
          return true;
        }
        else
        {
          return false;
        }
      }

      template< class ValueClass >
      inline std::string const&
      DenseSinglyIndexed< ValueClass >::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        for( size_t soughtIndex( 1 );
             valueVector.size() >= soughtIndex;
             ++soughtIndex )
        {
          this->indexPrintingVector[ 0 ] = soughtIndex;
          this->stringInterpretation.append( this->indicesToPrintingString() );
          // SLHA indices are in the sane starts-at-one format, while C++ code
          // uses the silly starts-at-zero format.
          this->stringInterpretation.append( this->valueToPrintingString(
                                            findOrMakeEntry( soughtIndex ) ) );
          this->stringInterpretation.append( "\n" );
        }
        return this->stringInterpretation;
      }

      template< class ValueClass >
      inline void
      DenseSinglyIndexed< ValueClass >::clearEntries()
      // this ensures that the entry at soughtIndex exists, filling out with
      // copies of defaultUnsetValue, & returns it.
      {
        valueVector.clear();
      }

      template< class ValueClass >
      inline ValueClass&
      DenseSinglyIndexed< ValueClass >::findOrMakeEntry(
                                                        int const soughtIndex )
      // this ensures that the entry at soughtIndex exists, filling out with
      // copies of defaultUnsetValue, & returns it.
      {
        if( valueVector.size() < (size_t)soughtIndex )
        {
          valueVector.resize( soughtIndex,
                              this->defaultUnsetValue );
        }
        return valueVector[ soughtIndex - 1 ];
      }

      template< class ValueClass >
      inline void
      DenseSinglyIndexed< ValueClass >::interpretCurrentStringBlock()
      {
        for( int whichLine( this->currentStringBlock->getNumberOfBodyLines() );
             0 < whichLine;
             --whichLine )
        {
          this->currentWord.assign( BOL::StringParser::firstWordOf(
                              (*(this->currentStringBlock))[ whichLine ].first,
                                                       &(this->lineRemainderA),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
          if( !(this->currentWord.empty()) )
          {
            recordingIndex
            = BOL::StringParser::stringToInt( this->currentWord );
            this->currentWord.assign( BOL::StringParser::trimFromFrontAndBack(
                                                          this->lineRemainderA,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            if( ( 0 < recordingIndex )
                &&
                !(this->currentWord.empty()) )
            {
              findOrMakeEntry( recordingIndex )
              = this->stringToValue( this->currentWord );
            }
            else if( this->isVerbose )
            {
              std::cout
              << std::endl
              << "LHPC::SLHA::error! expected to find 1 positive index then a"
              << " value, instead found \""
              << (*(this->currentStringBlock))[ whichLine ].first << "\"";
              std::cout << std::endl;
            }
          }
        }
      }

    }

  }

}

#endif /* DENSESINGLYINDEXED_HPP_ */
