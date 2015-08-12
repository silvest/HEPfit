/*
 * MultipleDoublyIndexed.hpp
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

#ifndef MULTIPLEDOUBLYINDEXED_HPP_
#define MULTIPLEDOUBLYINDEXED_HPP_

#include <map>
#include <list>
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
      class MultipleDoublyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        MultipleDoublyIndexed();
        virtual
        ~MultipleDoublyIndexed();

        std::list< ValueClass* >
        operator()( std::pair< int, int > const& indexPair );
        /* this returns a std::list of pointers to the ValueClass instances
         * mapped to by indexPair.first & indexPair.second. if there is no
         * element at the sought indices, an empty std::list is returned.
         */
        std::list< ValueClass const* >
        operator()( std::pair< int, int > const& indexPair ) const;
        // const version of above.
        std::list< ValueClass* >
        operator()( int const firstIndex,
                    int const secondIndex )
        { return (*this)( std::make_pair( firstIndex,
                                          secondIndex ) ); }
        std::list< ValueClass const* >
        operator()( int const firstIndex,
                    int const secondIndex ) const
        { return (*this)( std::make_pair( firstIndex,
                                          secondIndex ) ); }
        std::multimap< std::pair< int, int >, ValueClass >&
        getValueMap();
        std::multimap< std::pair< int, int >, ValueClass > const&
        getValueMap() const;
        bool
        hasEntry( std::pair< int, int > const& indexPair ) const;
        // this returns true if there is an entry at the sought indices.
        bool
        hasEntry( int const firstIndex,
                  int const secondIndex ) const
        { return hasEntry( std::make_pair( firstIndex,
                                           secondIndex ) ); }
        virtual std::string const&
        getAsString();
        // see base version's description.
        virtual void
        clearEntries();
        // derived classes should clear their interpreted values.


      protected:
        typedef typename
        std::multimap< std::pair< int, int >, ValueClass >::const_iterator
        mapIterator;

        std::multimap< std::pair< int, int >, ValueClass > valueMap;
        std::pair< std::pair< int, int >, ValueClass > valueRecorder;

        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      MultipleDoublyIndexed< ValueClass >::MultipleDoublyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueMap(),
          valueRecorder()
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      MultipleDoublyIndexed< ValueClass >::~MultipleDoublyIndexed()
      {
        // does nothing.
      }

      template< class ValueClass >
      inline std::list< ValueClass* >
      MultipleDoublyIndexed< ValueClass >::operator()(
                                       std::pair< int, int > const& indexPair )
      /* this returns a std::list of pointers to the ValueClass instances
       * mapped to by indexPair.first & indexPair.second. if there is no
       * element at the sought indices, an empty std::list is returned.
       */
      {
        std::list< ValueClass* > returnList;
        std::pair< mapIterator, mapIterator >
        rangeIterators( valueMap.equal_range( indexPair ) );
        while( rangeIterators.first != rangeIterators.second )
        {
          returnList.push_back( &(rangeIterators.first->second) );
          ++(rangeIterators.first);
        }
        return returnList;
      }

      template< class ValueClass >
      inline std::list< ValueClass const* >
      MultipleDoublyIndexed< ValueClass >::operator()(
                                 std::pair< int, int > const& indexPair ) const
      // const version of above.
      {
        std::list< ValueClass const* > returnList;
        std::pair< mapIterator, mapIterator >
        rangeIterators( valueMap.equal_range( indexPair ) );
        while( rangeIterators.first != rangeIterators.second )
        {
          returnList.push_back( &(rangeIterators.first->second) );
          ++(rangeIterators.first);
        }
        return returnList;
      }

      template< class ValueClass >
      inline std::multimap< std::pair< int, int >, ValueClass >&
      MultipleDoublyIndexed< ValueClass >::getValueMap()
      {
        return valueMap;
      }

      template< class ValueClass >
      inline std::multimap< std::pair< int, int >, ValueClass > const&
      MultipleDoublyIndexed< ValueClass >::getValueMap() const
      {
        return valueMap;
      }

      template< class ValueClass >
      inline bool
      MultipleDoublyIndexed< ValueClass >::hasEntry(
                                 std::pair< int, int > const& indexPair ) const
      // this returns true if there is an entry at soughtIndex.
      {
        return ( 0 < valueMap.count( indexPair ) );
      }

      template< class ValueClass >
      inline std::string const&
      MultipleDoublyIndexed< ValueClass >::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        mapIterator valueFinder( valueMap.begin() );
        while( valueFinder != valueMap.end() )
        {
          this->indexPrintingVector[ 0 ] = valueFinder->first.first;
          this->indexPrintingVector[ 1 ] = valueFinder->first.second;
          this->stringInterpretation.append( this->indicesToPrintingString() );
          this->stringInterpretation.append( this->valueToPrintingString(
                                                       valueFinder->second ) );
          /* negative particle codes can be avoided in any block, I think.
           * well, the exact format is only specified for the MASS block, &
           * in that case, negative particle codes can be avoided, since the
           * charge-conjugate just has the same mass. if the charge-conjugate
           * pair has a mass splitting, I think they are defined with 2
           * separate codes. if negative codes end up here, it's not my
           * fault...
           */
          this->stringInterpretation.append( "\n" );
          ++valueFinder;
        }
        return this->stringInterpretation;
      }

      template< class ValueClass >
      inline void
      MultipleDoublyIndexed< ValueClass >::clearEntries()
      {
        valueMap.clear();
      }

      template< class ValueClass >
      inline void
      MultipleDoublyIndexed< ValueClass >::interpretCurrentStringBlock()
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
            valueRecorder.first.first
            = BOL::StringParser::stringToInt( this->currentWord );
            this->currentWord.assign( BOL::StringParser::firstWordOf(
                                                          this->lineRemainderA,
                                                       &(this->lineRemainderB),
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            if( !(this->currentWord.empty()) )
            {
              valueRecorder.first.second
              = BOL::StringParser::stringToInt( this->currentWord );
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

  }

}

#endif /* MULTIPLEDOUBLYINDEXED_HPP_ */
