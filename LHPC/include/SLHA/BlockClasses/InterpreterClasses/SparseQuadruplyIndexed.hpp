/*
 * SparseQuadruplyIndexed.hpp
 *
 *  Created on: Mar 19, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SPARSEQUADRUPLYINDEXED_HPP_
#define SPARSEQUADRUPLYINDEXED_HPP_

#include <map>
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
      class SparseQuadruplyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        typedef typename std::pair< int, int > IntPair;
        typedef typename std::pair< IntPair, IntPair > IntQuadruple;
        SparseQuadruplyIndexed();
        virtual
        ~SparseQuadruplyIndexed();

        ValueClass&
        operator()( IntQuadruple const& indexQuadruple );
        /* this returns the ValueClass mapped to by indexQuadruple.first.first,
         * indexQuadruple.first.second, indexQuadruple.second.first &
         * indexPair.second.second. if there is no element at the sought
         * indices, a new one is made & copied from defaultUnsetValue.
         */
        ValueClass const&
        operator()( IntQuadruple const& indexQuadruple ) const;
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at the sought indices if there isn't
         * an entry there already.
         */
        ValueClass&
        operator()( int const firstIndex,
                    int const secondIndex,
                    int const thirdIndex,
                    int const fourthIndex )
        { return (*this)( std::make_pair( std::make_pair( firstIndex,
                                                          secondIndex ),
                                          std::make_pair( thirdIndex,
                                                          fourthIndex ) ) ); }
        ValueClass const&
        operator()( int const firstIndex,
                    int const secondIndex,
                    int const thirdIndex,
                    int const fourthIndex ) const
        { return (*this)( std::make_pair( std::make_pair( firstIndex,
                                                          secondIndex ),
                                          std::make_pair( thirdIndex,
                                                          fourthIndex ) ) ); }
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at the sought indices if there isn't
         * an entry there already.
         */
        bool
        hasEntry( IntQuadruple const& indexQuadruple ) const;
        // this returns true if there is an entry at the sought indices.
        bool
        hasEntry( int const firstIndex,
                  int const secondIndex,
                  int const thirdIndex,
                  int const fourthIndex ) const
        { return (*this)( hasEntry( std::make_pair( std::make_pair( firstIndex,
                                                                 secondIndex ),
                                                    std::make_pair( thirdIndex,
                                                        fourthIndex ) ) ) ) ; }
        virtual std::string const&
        getAsString();
        // see base version's description.
        virtual void
        clearEntries();
        // derived classes should clear their interpreted values.


      protected:
        typedef typename
        std::map< IntQuadruple, ValueClass >::const_iterator
        mapIterator;

        std::map< IntQuadruple, ValueClass > valueMap;
        IntQuadruple mapKey;
        std::pair< IntQuadruple, ValueClass > valueRecorder;

        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      SparseQuadruplyIndexed< ValueClass >::SparseQuadruplyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueMap(),
          mapKey( std::pair< int, int >( 0,
                                         0 ),
                  std::pair< int, int >( 0,
                                         0 ) ),
          valueRecorder()
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      SparseQuadruplyIndexed< ValueClass >::~SparseQuadruplyIndexed()
      {
        // does nothing.
      }


      template< class ValueClass >
      inline ValueClass&
      SparseQuadruplyIndexed< ValueClass >::operator()(
                                           IntQuadruple const& indexQuadruple )
      /* this returns the ValueClass mapped to by indexPair.first &
       * indexPair.second. if there is no element at the sought indices, a
       * new one is made & copied from defaultUnsetValue.
       */
      {
        if( 0 >= valueMap.count( indexQuadruple ) )
        {
          valueMap[ indexQuadruple ] = this->defaultUnsetValue;
        }
        return valueMap[ indexQuadruple ];
      }

      template< class ValueClass >
      inline ValueClass const&
      SparseQuadruplyIndexed< ValueClass >::operator()(
                                     IntQuadruple const& indexQuadruple ) const
      /* const version of above, though it returns defaultUnsetValue rather
       * than copying in a new element at the sought indices if there isn't
       * an entry there already.
       */
      {
        mapIterator valueFinder( valueMap.find( indexQuadruple ) );
        if( valueMap.end() != valueFinder )
        {
          return valueFinder->second;
        }
        else
        {
          return this->defaultUnsetValue;
        }
      }

      template< class ValueClass >
      inline bool
      SparseQuadruplyIndexed< ValueClass >::hasEntry(
                                     IntQuadruple const& indexQuadruple ) const
      // this returns true if there is an entry at soughtIndex.
      {
        return ( 0 < valueMap.count( indexQuadruple ) );
      }

      template< class ValueClass >
      inline std::string const&
      SparseQuadruplyIndexed< ValueClass >::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        mapIterator valueFinder( valueMap.begin() );
        while( valueFinder != valueMap.end() )
        {
          this->indexPrintingVector[ 0 ] = valueFinder->first.first.first;
          this->indexPrintingVector[ 1 ] = valueFinder->first.first.second;
          this->indexPrintingVector[ 2 ] = valueFinder->first.second.first;
          this->indexPrintingVector[ 3 ] = valueFinder->first.second.second;
          this->stringInterpretation.append( this->indicesToPrintingString() );
          this->stringInterpretation.append( this->valueToPrintingString(
                                                       valueFinder->second ) );
          this->stringInterpretation.append( "\n" );
          ++valueFinder;
        }
        return this->stringInterpretation;
      }

      template< class ValueClass >
      inline void
      SparseQuadruplyIndexed< ValueClass >::clearEntries()
      {
        valueMap.clear();
      }

      template< class ValueClass >
      inline void
      SparseQuadruplyIndexed< ValueClass >::interpretCurrentStringBlock()
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
            valueRecorder.first.first.first
            = BOL::StringParser::stringToInt( this->currentWord );
            this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                          this->lineRemainderA,
                                                       &(this->lineRemainderB),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            if( !(this->currentWord.empty()) )
            {
              valueRecorder.first.first.second
              = BOL::StringParser::stringToInt( this->currentWord );
              this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                        this->lineRemainderB,
                                                     &(this->lineRemainderA),
                            BOL::StringParser::whitespaceAndNewlineChars ) );
              if( !(this->currentWord.empty()) )
              {
                valueRecorder.first.second.first
                = BOL::StringParser::stringToInt( this->currentWord );
                this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                          this->lineRemainderA,
                                                       &(this->lineRemainderB),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
                if( !(this->currentWord.empty()) )
                {
                  valueRecorder.first.second.second
                  = BOL::StringParser::stringToInt( this->currentWord );
                  valueRecorder.second = this->stringToValue(
                                       BOL::StringParser::trimFromFrontAndBack(
                                                          this->lineRemainderB,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
                  valueMap.insert( valueRecorder );
                }
              }
            }
          }
        }
      }

    }

  }

}

#endif /* SPARSEQUADRUPLYINDEXED_HPP_ */
