/*
 * SparseQuadruplyIndexedBlock.hpp
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

#ifndef SPARSEQUADRUPLYINDEXEDBLOCK_HPP_
#define SPARSEQUADRUPLYINDEXEDBLOCK_HPP_

#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/SparseQuadruplyIndexed.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this template class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted as having a pair of int
     * indices which does not have to have entries for each value (nor even
     * necessarily positive index values).
     */
    template< class ValueClass >
    class SparseQuadruplyIndexedBlock : public IndexedBlockTemplate<
                                                                    ValueClass,
                       InterpreterClass::SparseQuadruplyIndexed< ValueClass > >
    {
    public:
      SparseQuadruplyIndexedBlock( std::string const& blockName,
                                   ValueClass const& defaultUnsetValue,
                                   bool const isVerbose = false,
                                   int const firstIndexDigits = 9,
                                   int const secondIndexDigits = 9,
                                   int const thirdIndexDigits = 2,
                                   int const fourthIndexDigits = 2 );
      virtual
      ~SparseQuadruplyIndexedBlock();

      ValueClass&
      operator()( std::pair< std::pair< int, int >,
                             std::pair< int, int > > const& indexQuadruple );
      // this returns operator() of the lowest-scale interpreter.
      ValueClass const&
      operator()( std::pair< std::pair< int, int >,
                         std::pair< int, int > > const& indexQuadruple ) const;
      // const version of above.
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
      bool
      hasEntry( std::pair< std::pair< int, int >,
                         std::pair< int, int > > const& indexQuadruple ) const;
      // this returns hasEntry( indexPair ) of the lowest-scale interpreter.
      bool
      hasEntry( int const firstIndex,
                int const secondIndex,
                int const thirdIndex,
                int const fourthIndex ) const
      { return hasEntry( std::make_pair( std::make_pair( firstIndex,
                                                         secondIndex ),
                                         std::make_pair( thirdIndex,
                                                         fourthIndex ) ) ); }
    };





    template< class ValueClass >
    inline
    SparseQuadruplyIndexedBlock< ValueClass >::SparseQuadruplyIndexedBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                    int const firstIndexDigits,
                                                   int const secondIndexDigits,
                                                    int const thirdIndexDigits,
                                                int const fourthIndexDigits ) :
        IndexedBlockTemplate< ValueClass,
                      InterpreterClass::SparseQuadruplyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                                                   BOL::Vi( firstIndexDigits )(
                                                           secondIndexDigits )(
                                                          thirdIndexDigits ).e(
                                                          fourthIndexDigits ) )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    SparseQuadruplyIndexedBlock< ValueClass >::~SparseQuadruplyIndexedBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline ValueClass&
    SparseQuadruplyIndexedBlock< ValueClass >::operator()(
                                              std::pair< std::pair< int, int >,
                                std::pair< int, int > > const& indexQuadruple )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexQuadruple );
    }

    template< class ValueClass >
    inline ValueClass const&
    SparseQuadruplyIndexedBlock< ValueClass >::operator()(
                                              std::pair< std::pair< int, int >,
                          std::pair< int, int > > const& indexQuadruple ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( indexQuadruple );
    }

    template< class ValueClass >
    inline bool
    SparseQuadruplyIndexedBlock< ValueClass >::hasEntry(
                                              std::pair< std::pair< int, int >,
                          std::pair< int, int > > const& indexQuadruple ) const
    // derived classes over-ride this to interpret their data as a
    // std::string.
    {
      return
      this->dataBlocks[ this->lowestScaleIndex ].hasEntry( indexQuadruple );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace


#endif /* SPARSEQUADRUPLYINDEXEDBLOCK_HPP_ */
