/*
 * SparseManyIndexed.hpp
 *
 *  Created on: Aug 28, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SPARSEMANYINDEXED_HPP_
#define SPARSEMANYINDEXED_HPP_

#include <map>
#include "IndexedInterpreter.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace InterpreterClass
    {
      // this template class interprets SLHA blocks that have a specified
      // number of int indices with a single ValueClass value.
      template< class ValueClass >
      class SparseManyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        SparseManyIndexed();
        virtual
        ~SparseManyIndexed();

        ValueClass&
        operator()( std::vector< int > const& indexVector );
        /* this returns the ValueClass mapped to by indexVector. if there is no
         * element at the sought indices, a new one is made & copied from
         * defaultUnsetValue.
         */
        ValueClass const&
        operator()( std::vector< int > const& indexVector ) const;
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at the sought indices if there isn't
         * an entry there already.
         */
        ValueClass&
        operator()( std::string const& indicesAsString );
        // this takes a list of indices in string form, separated by
        // commas or semicolons, with optional whitespace, or just whitespace.
        ValueClass const&
        operator()( std::string const& indicesAsString ) const;
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at the sought indices if there isn't
         * an entry there already.
         */
        bool
        hasEntry( std::vector< int > const& indexVector ) const;
        // this returns true if there is an entry at the sought indices.
        bool
        hasEntry( std::string const& indicesAsString ) const;
        // this returns true if there is an entry at the vector parsed from
        // indicesAsString see operator() for the format).
        virtual std::string const&
        getAsString();
        // see base version's description.
        virtual void
        clearEntries();
        // derived classes should clear their interpreted values.


      protected:
        typedef typename
        std::map< std::vector< int >, ValueClass >::const_iterator
        mapIterator;

        std::map< std::vector< int >, ValueClass > valueMap;
        std::vector< int > mapKey;
        std::pair< std::vector< int >, ValueClass > valueRecorder;

        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      SparseManyIndexed< ValueClass >::SparseManyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueMap(),
          mapKey(),
          valueRecorder()
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      SparseManyIndexed< ValueClass >::~SparseManyIndexed()
      {
        // does nothing.
      }


      template< class ValueClass >
      inline ValueClass&
      SparseManyIndexed< ValueClass >::operator()(
                                        std::vector< int > const& indexVector )
      /* this returns the ValueClass mapped to by indexVector. if there is no
       * element at the sought indices, a new one is made & copied from
       * defaultUnsetValue.
       */
      {
        if( 0 >= valueMap.count( indexVector ) )
        {
          valueMap[ indexVector ] = this->defaultUnsetValue;
        }
        return valueMap[ indexVector ];
      }

      template< class ValueClass >
      inline ValueClass const&
      SparseManyIndexed< ValueClass >::operator()(
                                  std::vector< int > const& indexVector ) const
      /* const version of above, though it returns defaultUnsetValue rather
       * than copying in a new element at the sought indices if there isn't
       * an entry there already.
       */
      {
        mapIterator valueFinder( valueMap.find( indexVector ) );
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
      inline ValueClass&
      SparseManyIndexed< ValueClass >::operator()(
                                           std::string const& indicesAsString )
      // this takes a list of indices in string form, separated by
      // commas or semicolons, with optional whitespace, or just whitespace.
      {
        return (*this)( BOL::StringParser::stringToIntVector(
                                                           indicesAsString ) );
      }


      template< class ValueClass >
      inline ValueClass const&
      SparseManyIndexed< ValueClass >::operator()(
                                     std::string const& indicesAsString ) const
      /* const version of above, though it returns defaultUnsetValue rather
       * than copying in a new element at the sought indices if there isn't
       * an entry there already.
       */
      {
        return (*this)( BOL::StringParser::stringToIntVector(
                                                           indicesAsString ) );
      }

      template< class ValueClass >
      inline bool
      SparseManyIndexed< ValueClass >::hasEntry(
                                  std::vector< int > const& indexVector ) const
      // this returns true if there is an entry at the sought indices.
      {
        return ( 0 < valueMap.count( indexVector ) );
      }

      template< class ValueClass >
      inline bool
      SparseManyIndexed< ValueClass >::hasEntry(
                                     std::string const& indicesAsString ) const
      // this returns true if there is an entry at the vector parsed from
      // indicesAsString see operator() for the format).
      {
        return hasEntry( BOL::StringParser::stringToIntVector(
                                                           indicesAsString ) );
      }

      template< class ValueClass >
      inline std::string const&
      SparseManyIndexed< ValueClass >::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        mapIterator valueFinder( valueMap.begin() );
        while( valueFinder != valueMap.end() )
        {
          this->indexPrintingVector = valueFinder->first;
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
      SparseManyIndexed< ValueClass >::clearEntries()
      {
        valueMap.clear();
      }

      template< class ValueClass >
      inline void
      SparseManyIndexed< ValueClass >::interpretCurrentStringBlock()
      {
        for( int whichLine( this->currentStringBlock->getNumberOfBodyLines() );
             0 < whichLine;
             --whichLine )
        {
          valueRecorder.first.clear();
          this->lineRemainderB.assign(
                            (*(this->currentStringBlock))[ whichLine ].first );
          for( int indexCount( this->indexDigitsVector.size() );
               0 < indexCount;
               --indexCount )
          {
            this->lineRemainderA.assign( this->lineRemainderB );
            this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                          this->lineRemainderA,
                                                       &(this->lineRemainderB),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            if( this->currentWord.empty() )
            {
              break;
            }
            valueRecorder.first.push_back(
                         BOL::StringParser::stringToInt( this->currentWord ) );
          }
          valueRecorder.second
          = this->stringToValue( BOL::StringParser::trimFromFrontAndBack(
                                                          this->lineRemainderB,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
          valueMap.insert( valueRecorder );
        }
      }

    }

  }

}



#endif /* SPARSEMANYINDEXED_HPP_ */
