/*
 * DenseDoublyIndexed.hpp
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

#ifndef DENSEDOUBLYINDEXED_HPP_
#define DENSEDOUBLYINDEXED_HPP_

#include "IndexedInterpreter.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace InterpreterClass
    {
      /* this template class derives from InterpreterTemplate to interpret
       * SLHA blocks that have a pair of int indices with a single ValueClass
       * value.
       */
      template< class ValueClass >
      class DenseDoublyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        DenseDoublyIndexed();
        virtual
        ~DenseDoublyIndexed();

        ValueClass&
        operator()( int const firstIndex,
                    int const secondIndex );
        /* this returns the ValueClass indexed by firstIndex & secondIndex. if
         * there is no element at the sought indices, valueMatrix is extended
         * with copies of defaultUnsetValue until there is an element at the
         * sought indices.
         */
        ValueClass const&
        operator()( int const firstIndex,
                    int const secondIndex ) const;
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at the sought indices if there isn't
         * an entry there already.
         */
        ValueClass&
        operator()( std::pair< int, int > const& indexPair )
        { return (*this)( indexPair.first,
                          indexPair.second ); }
        ValueClass const&
        operator()( std::pair< int, int > const& indexPair ) const
        { return (*this)( indexPair.first,
                          indexPair.second ); }
        bool
        hasEntry( int const firstIndex,
                  int const secondIndex ) const;
        // this returns true if there is an entry at the sought indices.
        bool
        hasEntry( std::pair< int, int > const& indexPair ) const
        { return hasEntry( indexPair.first,
                           indexPair.second ); }
        virtual std::string const&
        getAsString();
        // see base version's description.
        virtual void
        clearEntries();
        // derived classes should clear their interpreted values.


      protected:
        BOL::VectorlikeArray< std::vector< ValueClass > > valueMatrix;
        int firstRecordingIndex;
        size_t secondRecordingIndex;
        int largestFirstIndex;
        size_t largestSecondIndex;

        ValueClass&
        findOrMakeEntry( int firstIndex,
                         int const secondIndex );
        // this ensures that the entry at the sought indices exists, filling
        // out with copies of defaultUnsetValue, & returns it.
        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      DenseDoublyIndexed< ValueClass >::DenseDoublyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueMatrix(),
          firstRecordingIndex( 0 ),
          secondRecordingIndex( 0 ),
          largestFirstIndex( 0 ),
          largestSecondIndex( 0 )
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      DenseDoublyIndexed< ValueClass >::~DenseDoublyIndexed()
      {
        // does nothing.
      }


      template< class ValueClass >
      inline ValueClass&
      DenseDoublyIndexed< ValueClass >::operator()( int const firstIndex,
                                                    int const secondIndex )
      /* this returns the ValueClass indexed by firstIndex & secondIndex. if
       * there is no element at the sought indices, valueMatrix is extended
       * with copies of defaultUnsetValue until there is an element at the
       * sought indices.
       */
      {
        return findOrMakeEntry( firstIndex,
                                secondIndex );
      }

      template< class ValueClass >
      inline ValueClass const&
      DenseDoublyIndexed< ValueClass >::operator()( int const firstIndex,
                                                  int const secondIndex ) const
      /* const version of above, though it returns defaultUnsetValue rather
       * than copying in a new element at the sought indices if there isn't an
       * entry there already.
       */
      {
        if( hasEntry( firstIndex,
                      secondIndex ) )
        {
          return valueMatrix[ firstIndex - 1 ][ secondIndex - 1 ];
        }
        else
        {
          return this->defaultUnsetValue;
        }
      }


      template< class ValueClass >
      inline bool
      DenseDoublyIndexed< ValueClass >::hasEntry( int const firstIndex,
                                                  int const secondIndex ) const
      {
        if( ( 0 < firstIndex )
            &&
            ( 0 < secondIndex )
            &&
            ( firstIndex <= valueMatrix.getSize() )
            &&
            ( (size_t)secondIndex <= valueMatrix[ firstIndex - 1 ].size() ) )
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
      DenseDoublyIndexed< ValueClass >::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        for( int firstIndex( 1 );
             largestFirstIndex >= firstIndex;
             ++firstIndex )
        {
          for( size_t secondIndex( 1 );
               largestSecondIndex >= secondIndex;
               ++secondIndex )
          {
            this->indexPrintingVector[ 0 ] = firstIndex;
            this->indexPrintingVector[ 1 ] = secondIndex;
            this->stringInterpretation.append(
                                             this->indicesToPrintingString() );
            // SLHA indices are in the sane starts-at-one format, while C++
            // code uses the silly starts-at-zero format.
            this->stringInterpretation.append( this->valueToPrintingString(
                                                   findOrMakeEntry( firstIndex,
                                                             secondIndex ) ) );
            this->stringInterpretation.append( "\n" );
          }
        }
        return this->stringInterpretation;
      }

      template< class ValueClass >
      inline void
      DenseDoublyIndexed< ValueClass >::clearEntries()
      // this ensures that the entry at soughtIndex exists, filling out with
      // copies of defaultUnsetValue, & returns it.
      {
        valueMatrix.clearEntries();
      }

      template< class ValueClass >
      inline ValueClass&
      DenseDoublyIndexed< ValueClass >::findOrMakeEntry( int firstIndex,
                                                        int const secondIndex )
      // this ensures that the entry at the sought indices exists, filling
      // out with copies of defaultUnsetValue, & returns it.
      {
        if( firstIndex > largestFirstIndex )
        {
          largestFirstIndex = firstIndex;
        }
        if( (size_t)secondIndex > largestSecondIndex )
        {
          largestSecondIndex = secondIndex;
        }
        while( valueMatrix.getSize() < firstIndex )
        {
          valueMatrix.newEnd().clear();
          // empty std::vectors are added as necessary.
        }
        if( valueMatrix[ (--firstIndex) ].size() < (size_t)secondIndex )
          // the conditional does the job of converting from sane starts-at-one
          // format to silly starts-at-zero format.
        {
          valueMatrix[ firstIndex ].resize( secondIndex,
                                            this->defaultUnsetValue );
        }
        return valueMatrix[ firstIndex ][ secondIndex - 1 ];
      }

      template< class ValueClass >
      inline void
      DenseDoublyIndexed< ValueClass >::interpretCurrentStringBlock()
      {
        valueMatrix.clearEntries();
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
            firstRecordingIndex
            = BOL::StringParser::stringToInt( this->currentWord );
            this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                          this->lineRemainderA,
                                                       &(this->lineRemainderB),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            secondRecordingIndex
            = BOL::StringParser::stringToInt( this->currentWord );
            this->currentWord.assign( BOL::StringParser::trimFromFrontAndBack(
                                                          this->lineRemainderB,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            if( ( 0 < firstRecordingIndex )
                &&
                ( 0 < secondRecordingIndex )
                &&
                !(this->currentWord.empty()) )
            {
              findOrMakeEntry( firstRecordingIndex,
                               secondRecordingIndex )
              = this->stringToValue( this->currentWord );
            }
            else if( this->isVerbose )
            {
              std::cout
              << std::endl
              << "LHPC::SLHA::error! expected to find 2 positive indices then"
              << " a value, instead found \""
              << (*(this->currentStringBlock))[ whichLine ].first << "\"";
              std::cout << std::endl;
            }
          }
        }
      }

    }

  }

}

#endif /* DENSEDOUBLYINDEXED_HPP_ */
