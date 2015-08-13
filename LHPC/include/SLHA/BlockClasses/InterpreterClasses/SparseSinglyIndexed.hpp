/*
 * SparseSinglyIndexed.hpp
 *
 *  Created on: Feb 7, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SPARSESINGLYINDEXED_HPP_
#define SPARSESINGLYINDEXED_HPP_

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
      class SparseSinglyIndexed : public IndexedInterpreter< ValueClass >
      {
      public:
        SparseSinglyIndexed();
        virtual
        ~SparseSinglyIndexed();

        ValueClass&
        operator()( int const soughtIndex );
        /* this returns the ValueClass mapped to by soughtIndex. if there is no
         * element at soughtIndex, a new one is made & copied from
         * defaultUnsetValue.
         */
        ValueClass const&
        operator()( int const soughtIndex ) const;
        /* const version of above, though it returns defaultUnsetValue rather
         * than copying in a new element at soughtIndex if there isn't an entry
         * there already.
         */
        ValueClass&
        operator[]( int const soughtIndex )
        { return (*this)( soughtIndex ); }
        ValueClass const&
        operator[]( int const soughtIndex ) const
        { return (*this)( soughtIndex ); }
        std::map< int, ValueClass >&
        getValueMap();
        std::map< int, ValueClass > const&
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
        std::map< int, ValueClass >::const_iterator mapIterator;

        std::map< int, ValueClass > valueMap;
        std::pair< int, ValueClass > valueRecorder;

        virtual void
        interpretCurrentStringBlock();
      };





      template< class ValueClass >
      inline
      SparseSinglyIndexed< ValueClass >::SparseSinglyIndexed() :
          IndexedInterpreter< ValueClass >(),
          valueMap(),
          valueRecorder()
      {
        // just an initialization list.
      }

      template< class ValueClass >
      inline
      SparseSinglyIndexed< ValueClass >::~SparseSinglyIndexed()
      {
        // does nothing.
      }


      template< class ValueClass >
      inline ValueClass&
      SparseSinglyIndexed< ValueClass >::operator()( int const soughtIndex )
      /* this returns the ValueClass mapped to by soughtIndex for the data
       * with lowest energy scale. if there is no element at soughtIndex, a new
       * one is made & copied from defaultUnsetValue.
       */
      {
        if( 0 >= valueMap.count( soughtIndex ) )
        {
          valueMap[ soughtIndex ] = this->defaultUnsetValue;
        }
        return valueMap[ soughtIndex ];
      }

      template< class ValueClass >
      inline ValueClass const&
      SparseSinglyIndexed< ValueClass >::operator()(
                                                  int const soughtIndex ) const
      /* const version of above, though it returns defaultUnsetValue rather
       * than copying in a new element at soughtIndex if there isn't an entry
       * there already.
       */
      {
        mapIterator valueFinder( valueMap.find( soughtIndex ) );
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
      inline std::map< int, ValueClass >&
      SparseSinglyIndexed< ValueClass >::getValueMap()
      {
        return valueMap;
      }

      template< class ValueClass >
      inline std::map< int, ValueClass > const&
      SparseSinglyIndexed< ValueClass >::getValueMap() const
      {
        return valueMap;
      }

      template< class ValueClass >
      inline bool
      SparseSinglyIndexed< ValueClass >::hasEntry(
                                                  int const soughtIndex ) const
      // this returns true if there is an entry at soughtIndex.
      {
        return ( 0 < valueMap.count( soughtIndex ) );
      }

      template< class ValueClass >
      inline std::string const&
      SparseSinglyIndexed< ValueClass >::getAsString()
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
      SparseSinglyIndexed< ValueClass >::clearEntries()
      {
        valueMap.clear();
      }

      template< class ValueClass >
      inline void
      SparseSinglyIndexed< ValueClass >::interpretCurrentStringBlock()
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

#endif /* SPARSESINGLYINDEXED_HPP_ */
