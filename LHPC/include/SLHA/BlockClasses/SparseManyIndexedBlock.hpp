/*
 * SparseManyIndexedBlock.hpp
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

#ifndef SPARSEMANYINDEXEDBLOCK_HPP_
#define SPARSEMANYINDEXEDBLOCK_HPP_


#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/SparseManyIndexed.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this template class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted as having a set of int
     * indices which does not have to have entries for each value (nor even
     * necessarily positive index values).
     */
    template< class ValueClass >
    class SparseManyIndexedBlock : public IndexedBlockTemplate< ValueClass,
                            InterpreterClass::SparseManyIndexed< ValueClass > >
    {
    public:
      SparseManyIndexedBlock( std::string const& blockName,
                              int const numberOfIndices,
                              ValueClass const& defaultUnsetValue,
                              bool const isVerbose = false,
                              int const indexDigits = 9 );
      virtual
      ~SparseManyIndexedBlock();

      ValueClass&
      operator()( std::vector< int > const& indexVector );
      // this returns operator() of the lowest-scale interpreter.
      ValueClass const&
      operator()( std::vector< int > const& indexVector ) const;
      // const version of above.
      ValueClass&
      operator()( std::string const& indicesAsString );
      // this returns operator() of the lowest-scale interpreter.
      ValueClass const&
      operator()( std::string const& indicesAsString ) const;
      // const version of above.
      bool
      hasEntry( std::vector< int > const& indexVector ) const;
      // this returns hasEntry( indexVector ) of the lowest-scale interpreter.
      bool
      hasEntry( std::string const& indicesAsString ) const;
      // this returns hasEntry( indicesAsString ) of the lowest-scale
      // interpreter.
    };





    template< class ValueClass >
    inline
    SparseManyIndexedBlock< ValueClass >::SparseManyIndexedBlock(
                                                  std::string const& blockName,
                                                     int const numberOfIndices,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                      int const indexDigits ) :
        IndexedBlockTemplate< ValueClass,
                           InterpreterClass::SparseManyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                                           std::vector< int >( numberOfIndices,
                                                               indexDigits ) )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    SparseManyIndexedBlock< ValueClass >::~SparseManyIndexedBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline ValueClass&
    SparseManyIndexedBlock< ValueClass >::operator()(
                                        std::vector< int > const& indexVector )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexVector );
    }

    template< class ValueClass >
    inline ValueClass const&
    SparseManyIndexedBlock< ValueClass >::operator()(
                                  std::vector< int > const& indexVector ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexVector );
    }

    template< class ValueClass >
    inline ValueClass&
    SparseManyIndexedBlock< ValueClass >::operator()(
                                           std::string const& indicesAsString )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indicesAsString );
    }

    template< class ValueClass >
    inline ValueClass const&
    SparseManyIndexedBlock< ValueClass >::operator()(
                                     std::string const& indicesAsString ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indicesAsString );
    }

    template< class ValueClass >
    inline bool
    SparseManyIndexedBlock< ValueClass >::hasEntry(
                                  std::vector< int > const& indexVector ) const
    // derived classes over-ride this to interpret their data as a
    // std::string.
    {
      return
      this->dataBlocks[ this->lowestScaleIndex ].hasEntry( indexVector );
    }

    template< class ValueClass >
    inline bool
    SparseManyIndexedBlock< ValueClass >::hasEntry(
                                     std::string const& indicesAsString ) const
    // derived classes over-ride this to interpret their data as a
    // std::string.
    {
      return
      this->dataBlocks[ this->lowestScaleIndex ].hasEntry( indicesAsString );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace



#endif /* SPARSEMANYINDEXEDBLOCK_HPP_ */
