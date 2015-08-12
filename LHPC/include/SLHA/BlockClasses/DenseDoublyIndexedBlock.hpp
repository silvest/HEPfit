/*
 * DenseDoublyIndexedBlock.hpp
 *
 *  Created on: Mar 12, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef DENSEDOUBLYINDEXEDBLOCK_HPP_
#define DENSEDOUBLYINDEXEDBLOCK_HPP_

#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/DenseDoublyIndexed.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this template class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted as having a pair of int
     * indices which have to have entries for each value, each index beginning
     * at 1.
     */
    template< class ValueClass >
    class DenseDoublyIndexedBlock : public IndexedBlockTemplate< ValueClass,
                           InterpreterClass::DenseDoublyIndexed< ValueClass > >
    {
    public:
      DenseDoublyIndexedBlock( std::string const& blockName,
                               ValueClass const& defaultUnsetValue,
                               bool const isVerbose = false,
                               int const firstIndexDigits = 2,
                               int const secondIndexDigits = 2 );
      virtual
      ~DenseDoublyIndexedBlock();

      ValueClass&
      operator()( int const firstIndex,
                  int const secondIndex );
      // this returns operator() of the lowest-scale interpreter.
      ValueClass const&
      operator()( int const firstIndex,
                  int const secondIndex ) const;
      // const version of above.
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
      // this returns hasEntry( firstIndex, secondIndex ) of the
      // lowest-scale interpreter.
      bool
      hasEntry( std::pair< int, int > const& indexPair ) const
      { return hasEntry( indexPair.first,
                         indexPair.second ); }
    };





    template< class ValueClass >
    inline
    DenseDoublyIndexedBlock< ValueClass >::DenseDoublyIndexedBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                    int const firstIndexDigits,
                                                int const secondIndexDigits ) :
        IndexedBlockTemplate< ValueClass,
                          InterpreterClass::DenseDoublyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                           BOL::Vi( firstIndexDigits ).e( secondIndexDigits ) )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    DenseDoublyIndexedBlock< ValueClass >::~DenseDoublyIndexedBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline ValueClass&
    DenseDoublyIndexedBlock< ValueClass >::operator()( int const firstIndex,
                                                       int const secondIndex )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( firstIndex,
                                                         secondIndex );
    }

    template< class ValueClass >
    inline ValueClass const&
    DenseDoublyIndexedBlock< ValueClass >::operator()( int const firstIndex,
                                                  int const secondIndex ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( firstIndex,
                                                         secondIndex );
    }

    template< class ValueClass >
    inline bool
    DenseDoublyIndexedBlock< ValueClass >::hasEntry( int const firstIndex,
                                                  int const secondIndex ) const
    // derived classes over-ride this to interpret their data as a
    // std::string.
    {
      return this->dataBlocks[ this->lowestScaleIndex ].hasEntry( firstIndex,
                                                                 secondIndex );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace


#endif /* DENSEDOUBLYINDEXEDBLOCK_HPP_ */
