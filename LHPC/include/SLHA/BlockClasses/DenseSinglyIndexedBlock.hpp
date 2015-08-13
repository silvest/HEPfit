/*
 * DenseSinglyIndexedBlock.hpp
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

#ifndef DENSESINGLYINDEXEDBLOCK_HPP_
#define DENSESINGLYINDEXEDBLOCK_HPP_

#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/DenseSinglyIndexed.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this template class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted as having a single int
     * index which has to have entries for each value, beginning with index 1.
     */
    template< class ValueClass >
    class DenseSinglyIndexedBlock : public IndexedBlockTemplate< ValueClass,
                           InterpreterClass::DenseSinglyIndexed< ValueClass > >
    {
    public:
      DenseSinglyIndexedBlock( std::string const& blockName,
                               ValueClass const& defaultUnsetValue,
                               bool const isVerbose = false,
                               int const indexDigits = 5 );
      virtual
      ~DenseSinglyIndexedBlock();

      ValueClass&
      operator()( int const soughtIndex );
      // this returns operator() of the lowest-scale interpreter.
      ValueClass const&
      operator()( int const soughtIndex ) const;
      // const version of above.
      bool
      hasEntry( int const soughtIndex ) const;
      // this returns hasEntry( soughtIndex ) of the lowest-scale interpreter.
    };





    template< class ValueClass >
    inline
    DenseSinglyIndexedBlock< ValueClass >::DenseSinglyIndexedBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                      int const indexDigits ) :
        IndexedBlockTemplate< ValueClass,
                          InterpreterClass::DenseSinglyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                                                         std::vector< int >( 1,
                                                                indexDigits ) )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    DenseSinglyIndexedBlock< ValueClass >::~DenseSinglyIndexedBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline ValueClass&
    DenseSinglyIndexedBlock< ValueClass >::operator()( int const soughtIndex )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( soughtIndex );
    }

    template< class ValueClass >
    inline ValueClass const&
    DenseSinglyIndexedBlock< ValueClass >::operator()(
                                                  int const soughtIndex ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( soughtIndex );
    }

    template< class ValueClass >
    inline bool
    DenseSinglyIndexedBlock< ValueClass >::hasEntry(
                                                  int const soughtIndex ) const
    // derived classes over-ride this to interpret their data as a
    // std::string.
    {
      return
      this->dataBlocks[ this->lowestScaleIndex ].hasEntry( soughtIndex );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace


#endif /* DENSESINGLYINDEXEDBLOCK_HPP_ */
