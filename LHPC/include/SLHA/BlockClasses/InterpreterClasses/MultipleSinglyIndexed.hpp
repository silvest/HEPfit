/*
 * MultipleSinglyIndexed.hpp
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

#ifndef MULTIPLESINGLYINDEXED_HPP_
#define MULTIPLESINGLYINDEXED_HPP_

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
      class MultipleSinglyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        MultipleSinglyIndexed();
        virtual
        ~MultipleSinglyIndexed();

        std::list< ValueClass* >
        operator()( int const soughtIndex );
        /* this returns a std::list of pointers to the ValueClass instances
         * mapped to by soughtIndex. if there is no element at soughtIndex, an
         * empty std::list is returned.
         */
        std::list< ValueClass const* >
        operator()( int const soughtIndex ) const;
        // const version of above.
        std::list< ValueClass* >
        operator[]( int const soughtIndex )
        { return (*this)( soughtIndex ); }
        std::list< ValueClass* > const
        operator[]( int const soughtIndex ) const
        { return (*this)( soughtIndex ); }
        std::multimap< int, ValueClass >&
        getValueMap();
        std::multimap< int, ValueClass > const&
        getValueMap() const;
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
        typedef typename
        std::multimap< int, ValueClass >::const_iterator mapIterator;

        std::multimap< int, ValueClass > valueMap;
        std::pair< int, ValueClass > valueRecorder;

        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      MultipleSinglyIndexed< ValueClass >::MultipleSinglyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueMap(),
          valueRecorder()
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      MultipleSinglyIndexed< ValueClass >::~MultipleSinglyIndexed()
      {
        // does nothing.
      }

      template< class ValueClass >
      inline std::list< ValueClass* >
      MultipleSinglyIndexed< ValueClass >::operator()( int const soughtIndex )
      /* this returns a std::list of pointers to the ValueClass instances
       * mapped to by soughtIndex. if there is no element at soughtIndex, an
       * empty std::list is returned.
       */
      {
        std::list< ValueClass* > returnList;
        std::pair< mapIterator, mapIterator >
        rangeIterators( valueMap.equal_range( soughtIndex ) );
        while( rangeIterators.first != rangeIterators.second )
        {
          returnList.push_back( &(rangeIterators.first->second) );
          ++(rangeIterators.first);
        }
        return returnList;
      }

      template< class ValueClass >
      inline std::list< ValueClass const* >
      MultipleSinglyIndexed< ValueClass >::operator()(
                                                  int const soughtIndex ) const
      // const version of above.
      {
        std::list< ValueClass const* > returnList;
        std::pair< mapIterator, mapIterator >
        rangeIterators( valueMap.equal_range( soughtIndex ) );
        while( rangeIterators.first != rangeIterators.second )
        {
          returnList.push_back( &(rangeIterators.first->second) );
          ++(rangeIterators.first);
        }
        return returnList;
      }

      template< class ValueClass >
      inline std::multimap< int, ValueClass >&
      MultipleSinglyIndexed< ValueClass >::getValueMap()
      {
        return valueMap;
      }

      template< class ValueClass >
      inline std::multimap< int, ValueClass > const&
      MultipleSinglyIndexed< ValueClass >::getValueMap() const
      {
        return valueMap;
      }

      template< class ValueClass >
      inline bool
      MultipleSinglyIndexed< ValueClass >::hasEntry(
                                                  int const soughtIndex ) const
      // this returns true if there is an entry at soughtIndex.
      {
        return ( 0 < valueMap.count( soughtIndex ) );
      }

      template< class ValueClass >
      inline std::string const&
      MultipleSinglyIndexed< ValueClass >::getAsString()
      // see base version's description.
      {
        this->stringInterpretation.clear();
        mapIterator valueFinder( valueMap.begin() );
        while( valueFinder != valueMap.end() )
        {
          this->indexPrintingVector[ 0 ] = valueFinder->first;
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
      MultipleSinglyIndexed< ValueClass >::clearEntries()
      {
        valueMap.clear();
      }

      template< class ValueClass >
      inline void
      MultipleSinglyIndexed< ValueClass >::interpretCurrentStringBlock()
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
            valueRecorder.first
            = BOL::StringParser::stringToInt( this->currentWord );
            valueRecorder.second
            = this->stringToValue( BOL::StringParser::trimFromFrontAndBack(
                                                          this->lineRemainderA,
                              BOL::StringParser::whitespaceAndNewlineChars ) );
            valueMap.insert( valueRecorder );
          }
        }
      }

    }

  }

}

#endif /* MULTIPLESINGLYINDEXED_HPP_ */
