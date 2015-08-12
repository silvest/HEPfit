/*
 * DoublyIndexedMultipleEntriesBlock.hpp
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

#ifndef DOUBLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_
#define DOUBLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_

#include "../../MEC/ExtendedMass.hpp"
#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/MultipleDoublyIndexed.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this template class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted as having a single int
     * index which does not have to have entries for each value (nor even
     * necessarily positive index values), allowing for multiple entries with
     * the same index.
     */
    template< class ValueClass >
    class DoublyIndexedMultipleEntriesBlock : public IndexedBlockTemplate<
                                                                    ValueClass,
                        InterpreterClass::MultipleDoublyIndexed< ValueClass > >
    {
    public:
      DoublyIndexedMultipleEntriesBlock( std::string const& blockName,
                                         ValueClass const& defaultUnsetValue,
                                         bool const isVerbose = false,
                                         int const firstIndexDigits = 5,
                                         int const secondIndexDigits = 5 );
      virtual
      ~DoublyIndexedMultipleEntriesBlock();

      std::list< ValueClass* >
      operator()( std::pair< int, int > const& indexPair );
      // this returns operator() of the lowest-scale interpreter.
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
      bool
      hasEntry( std::pair< int, int > const& indexPair ) const;
      // this returns hasEntry( indexPair ) of the lowest-scale interpreter.
      bool
      hasEntry( int const firstIndex,
                int const secondIndex ) const
      { return hasEntry( std::make_pair( firstIndex,
                                         secondIndex ) ); }
    };





    template< class ValueClass >
    inline
    DoublyIndexedMultipleEntriesBlock< ValueClass
                                          >::DoublyIndexedMultipleEntriesBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                    int const firstIndexDigits,
                                                int const secondIndexDigits ) :
        IndexedBlockTemplate< ValueClass,
                       InterpreterClass::MultipleDoublyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                           BOL::Vi( firstIndexDigits ).e( secondIndexDigits ) )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    DoublyIndexedMultipleEntriesBlock< ValueClass
                                        >::~DoublyIndexedMultipleEntriesBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline std::list< ValueClass* >
    DoublyIndexedMultipleEntriesBlock< ValueClass >::operator()(
                                       std::pair< int, int > const& indexPair )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexPair );
    }

    template< class ValueClass >
    inline std::list< ValueClass const* >
    DoublyIndexedMultipleEntriesBlock< ValueClass >::operator()(
                                 std::pair< int, int > const& indexPair ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexPair );
    }

    template< class ValueClass >
    inline bool
    DoublyIndexedMultipleEntriesBlock< ValueClass >::hasEntry(
                                 std::pair< int, int > const& indexPair ) const
    // this returns hasEntry( soughtIndex ) of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ].hasEntry( indexPair );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace

#endif /* DOUBLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_ */
