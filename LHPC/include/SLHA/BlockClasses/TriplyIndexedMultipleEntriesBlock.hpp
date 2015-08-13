/*
 * TriplyIndexedMultipleEntriesBlock.hpp
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

#ifndef TRIPLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_
#define TRIPLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_

#include "../../MEC/ExtendedMass.hpp"
#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/MultipleTriplyIndexed.hpp"

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
    class TriplyIndexedMultipleEntriesBlock : public IndexedBlockTemplate<
                                                                    ValueClass,
                        InterpreterClass::MultipleTriplyIndexed< ValueClass > >
    {
    public:
      TriplyIndexedMultipleEntriesBlock( std::string const& blockName,
                                         ValueClass const& defaultUnsetValue,
                                         bool const isVerbose = false,
                                         int const firstIndexDigits = 5,
                                         int const secondIndexDigits = 5,
                                         int const thirdIndexDigits = 5 );
      virtual
      ~TriplyIndexedMultipleEntriesBlock();

      std::list< ValueClass* >
      operator()( std::pair< std::pair< int, int >, int > const& indexTriple );
      // this returns operator() of the lowest-scale interpreter.
      std::list< ValueClass const* >
      operator()(
            std::pair< std::pair< int, int >, int > const& indexTriple ) const;
      // const version of above.
      std::list< ValueClass* >
      operator()( int const firstIndex,
                  int const secondIndex,
                  int const thirdIndex )
      { return (*this)( std::make_pair( std::make_pair( firstIndex,
                                                        secondIndex ),
                                        thirdIndex ) ); }
      std::list< ValueClass const* >
      operator()( int const firstIndex,
                  int const secondIndex,
                  int const thirdIndex ) const
      { return (*this)( std::make_pair( std::make_pair( firstIndex,
                                                        secondIndex ),
                                        thirdIndex ) ); }
      bool
      hasEntry(
            std::pair< std::pair< int, int >, int > const& indexTriple ) const;
      // this returns hasEntry( indexTriple ) of the lowest-scale interpreter.
      bool
      hasEntry( int const firstIndex,
                int const secondIndex,
                int const thirdIndex ) const
      { return hasEntry( std::make_pair( std::make_pair( firstIndex,
                                                         secondIndex ),
                                         thirdIndex ) ); }
    };





    template< class ValueClass >
    inline
    TriplyIndexedMultipleEntriesBlock< ValueClass
                                          >::TriplyIndexedMultipleEntriesBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                    int const firstIndexDigits,
                                                   int const secondIndexDigits,
                                                 int const thirdIndexDigits ) :
        IndexedBlockTemplate< ValueClass,
                       InterpreterClass::MultipleTriplyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                                                      BOL::Vi( firstIndexDigits
                                                           )( secondIndexDigits
                                                      ).e( thirdIndexDigits ) )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    TriplyIndexedMultipleEntriesBlock< ValueClass
                                        >::~TriplyIndexedMultipleEntriesBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline std::list< ValueClass* >
    TriplyIndexedMultipleEntriesBlock< ValueClass >::operator()(
                   std::pair< std::pair< int, int >, int > const& indexTriple )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexTriple );
    }

    template< class ValueClass >
    inline std::list< ValueClass const* >
    TriplyIndexedMultipleEntriesBlock< ValueClass >::operator()(
             std::pair< std::pair< int, int >, int > const& indexTriple ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexTriple );
    }

    template< class ValueClass >
    inline bool
    TriplyIndexedMultipleEntriesBlock< ValueClass >::hasEntry(
             std::pair< std::pair< int, int >, int > const& indexTriple ) const
    // this returns hasEntry( soughtIndex ) of the lowest-scale interpreter.
    {
      return
      this->dataBlocks[ this->lowestScaleIndex ].hasEntry( indexTriple );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace

#endif /* TRIPLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_ */
