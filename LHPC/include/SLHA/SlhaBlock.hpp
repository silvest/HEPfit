/*
 * SlhaBlock.hpp
 *
 *  Created on: Feb 1, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SLHABLOCK_HPP_
#define SLHABLOCK_HPP_

#include <string>
#include "BOLlib/include/BOLlib.hpp"
#include "BlockClasses/BaseSlhaBlock.hpp"
#include "BlockClasses/BaseStringBlock.hpp"
#include "BlockClasses/InterpreterClasses/BlockInterpreter.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* issues for specific blocks:
     * MODSEL: annoying mix of ints & doubles (everything but Q_max is an int).
     *       - set it as a string block & interpret with the general filters?
     *       - set it as a double block, but write a separate routine to
     *         write it in a way SPheno can understand?
     *       - set it as a double block, but get Werner to fix SPheno so it's
     *         not so silly?
     * SPINFO: multiple lines with the same index.
     *       - leave it as a SlhaBlock, because it really has no useful
     *         information
     * SPHENOINPUT: multiple lines with the same index.
     *            - leave it as a SlhaBlock, because it really has no
     *              useful information
     * SPHENOCROSSSECTIONS: just stupid.
     *                    - leave it as a SlhaBlock, to be interpreted by
     *                      the general filters.
     */


    /* this is a class to hold a block from a file in the SLHA format, only
     * interpreting it at the basic level: as set of strings. if multiple
     * blocks with the same name but different scales are recorded, the copy
     * with the lowest scale is defaulted to if no scale is given when seeking
     * information.
     * classes which interpret the strings further derive from this class.
     */
    template< class ValueClass, class BlockParser >
    class SlhaBlock : public BaseSlhaBlock
    {
    public:
      SlhaBlock( std::string const& blockName,
                 ValueClass const& defaultUnsetValue,
                 bool const isVerbose );
      virtual
      ~SlhaBlock();

      BlockParser&
      operator[]( int whichScaleIndex );
      /* the interpreters are indexed in the order in which they were
       * read (or created). the index starts at 1 rather than 0, as well. if 0
       * is given as the argument, the index corresponding to the copy with the
       * lowest scale is used.
       */
      BlockParser const&
      operator[]( int whichScaleIndex ) const;
      // const version of above.
      int
      getNumberOfCopiesWithDifferentScale() const;
      bool
      hasRecordedScale( double const soughtScale,
                        int& indexForLowerScale,
                        int& indexForUpperScale,
                        double& fractionFromLowerScale );
      /* this looks for the pair of blocks with energy scales closest to
       * soughtScale.
       * if there are no recorded copies of this block, none of the references
       * given are changed & false is returned.
       * if there is only one copy of the block, both indexForLowerScale &
       * indexForUpperScale are set to 1, fractionFromLowerScale is set to NaN,
       * & true is returned.
       * if there are 2 or more copies of the block, indexForLowerScale &
       * indexForUpperScale are set as described below, fractionFromLowerScale
       * is set to be
       * ( ( soughtScale - [ scale of copy with lower scale ] )
       *   / [ difference of copy scales ] ), & true is returned.
       * fractionFromLowerScale will thus be between 0.0 & 1.0 if there are
       * copies with scales above & below soughtScale, but may be negative if
       * soughtScale is lower than the lowest scale of the copies, or greater
       * than 1.0 if soughtScale is higher than the highest scale of the
       * copies.
       * since the copies of the block with different scales are not
       * necessarily recorded in order of scale, indexForLowerScale will not
       * necessarily be smaller in value than indexForUpperScale.
       * indexForLowerScale & indexForUpperScale are set as follows:
       * if soughtScale is lower than the lowest scale of the copies,
       * indexForLowerScale is set to the index of the copy with lowest scale,
       * & indexForUpperScale is set to the index of the copy with the next
       * lowest scale.
       * if soughtScale is higher than the highest scale of the copies,
       * indexForUpperScale is set to the index of the copy with highest scale,
       * & indexForLowerScale is set to the index of the copy with the next
       * highest scale.
       * otherwise, indexForLowerScale is set to the index of the copy with
       * highest scale which is still lower than soughtScale, &
       * indexForUpperScale is set to the index of the copy with lowest scale
       * which is still higher than soughtScale.
       */
      virtual std::string const&
      interpretAsString( bool onlyShowScalesGreaterThanZero = true );
      /* derived classes should override this to form their strings directly
       * from their data (because the block may be being used as a way of
       * writing an input file in the SLHA format, or maybe because the
       * original formatting was incorrect, as it is in most spectrum
       * generators...). by default it just joins together the strings from
       * blocksAsSingleStrings. each different-scale copy is concatenated. if
       * onlyShowScalesGreaterThanZero is true, blocks with scales of 0.0 or
       * less just do not have the "Q=" etc. printed.
       */
      virtual void
      respondToObservedSignal();
      /* by default, PushedToObserver instances don't respond without a pushed
       * value, but this is over-ridden to clear all the data that this block
       * has interpreted or had assigned.
       */
      virtual void
      respondToPush( BlockClass::BaseStringBlock const& pushedValue );
      // this adds a new BlockParser to dataBlocks & tells it to interpret
      // pushedValue. it also sorts out scale indices.


    protected:
      bool const isVerbose;
      ValueClass const defaultUnsetValue;
      BOL::VectorlikeArray< BlockParser > dataBlocks;
      int lowestScaleIndex;
      bool hasHadPushSinceLastReset;
      std::list< std::pair< int, double > > scaleOrderedIndices;
      std::list< std::pair< int, double > >::iterator scaleIndexIterator;
      std::string stringInterpretation;

      virtual void
      prepareNewDataBlock();
      // derived classes can insert extra arguments to any new block
      // interpreters by over-riding this.
    };





    template< class ValueClass, class BlockParser >
    inline
    SlhaBlock< ValueClass, BlockParser >::SlhaBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                     bool const isVerbose ) :
        BaseSlhaBlock( blockName ),
        isVerbose( isVerbose ),
        defaultUnsetValue( defaultUnsetValue ),
        dataBlocks( 1 ),
        lowestScaleIndex( 0 ),
        hasHadPushSinceLastReset( true ),
        scaleOrderedIndices(),
        scaleIndexIterator(),
        stringInterpretation( "" )
    {
      dataBlocks.getFront().setDefaultUnsetValue( defaultUnsetValue );
      dataBlocks.getFront().setVerbosity( isVerbose );
    }

    template< class ValueClass, class BlockParser >
    inline
    SlhaBlock< ValueClass, BlockParser >::~SlhaBlock()
    {
      // does nothing.
    }


    template< class ValueClass, class BlockParser >
    inline BlockParser&
    SlhaBlock< ValueClass, BlockParser >::operator[]( int whichScaleIndex )
    /* the interpreters are indexed in the order in which they were
     * read (or created). the index starts at 1 rather than 0, as well. if 0
     * is given as the argument, the index corresponding to the copy with the
     * lowest scale is used.
     */
    {
      if( 0 == whichScaleIndex )
      {
        return dataBlocks[ lowestScaleIndex ];
      }
      else
      {
        return dataBlocks[ (--whichScaleIndex) ];
      }
    }

    template< class ValueClass, class BlockParser >
    inline BlockParser const&
    SlhaBlock< ValueClass, BlockParser >::operator[](
                                                    int whichScaleIndex ) const
    // const version of above.
    {
      if( 0 == whichScaleIndex )
      {
        return dataBlocks[ lowestScaleIndex ];
      }
      else
      {
        return dataBlocks[ (--whichScaleIndex) ];
      }
    }

    template< class ValueClass, class BlockParser >
    inline int
    SlhaBlock< ValueClass, BlockParser >::getNumberOfCopiesWithDifferentScale(
                                                                        ) const
    {
      if( hasHadPushSinceLastReset )
      {
        return dataBlocks.getSize();
      }
      else
      {
        return 0;
      }
    }

    template< class ValueClass, class BlockParser >
    inline bool
    SlhaBlock< ValueClass, BlockParser >::hasRecordedScale(
                                                      double const soughtScale,
                                                       int& indexForLowerScale,
                                                       int& indexForUpperScale,
                                               double& fractionFromLowerScale )
    /* this looks for the pair of blocks with energy scales closest to
     * soughtScale.
     * if there are no recorded copies of this block, none of the references
     * given are changed & false is returned.
     * if there is only one copy of the block, both indexForLowerScale &
     * indexForUpperScale are set to 1, fractionFromLowerScale is set to NaN,
     * & true is returned.
     * if there are 2 or more copies of the block, indexForLowerScale &
     * indexForUpperScale are set as described below, fractionFromLowerScale
     * is set to be
     * ( ( soughtScale - [ scale of copy with lower scale ] )
     *   / [ difference of copy scales ] ), & true is returned.
     * fractionFromLowerScale will thus be between 0.0 & 1.0 if there are
     * copies with scales above & below soughtScale, but may be negative if
     * soughtScale is lower than the lowest scale of the copies, or greater
     * than 1.0 if soughtScale is higher than the highest scale of the
     * copies.
     * since the copies of the block with different scales are not
     * necessarily recorded in order of scale, indexForLowerScale will not
     * necessarily be smaller in value than indexForUpperScale.
     * indexForLowerScale & indexForUpperScale are set as follows:
     * if soughtScale is lower than the lowest scale of the copies,
     * indexForLowerScale is set to the index of the copy with lowest scale,
     * & indexForUpperScale is set to the index of the copy with the next
     * lowest scale.
     * if soughtScale is higher than the highest scale of the copies,
     * indexForUpperScale is set to the index of the copy with highest scale,
     * & indexForLowerScale is set to the index of the copy with the next
     * highest scale.
     * otherwise, indexForLowerScale is set to the index of the copy with
     * highest scale which is still lower than soughtScale, &
     * indexForUpperScale is set to the index of the copy with lowest scale
     * which is still higher than soughtScale.
     */
    {
      return BlockClass::BaseStringBlock::findScaleIndices( soughtScale,
                                                           scaleOrderedIndices,
                                                            indexForLowerScale,
                                                            indexForUpperScale,
                                                      fractionFromLowerScale );
    }

    template< class ValueClass, class BlockParser >
    inline std::string const&
    SlhaBlock< ValueClass, BlockParser >::interpretAsString(
                                           bool onlyShowScalesGreaterThanZero )
    /* derived classes should override getStringForScale to form their
     * strings directly from their data (because the block may be being used as
     * a way of writing an input file in the SLHA format, or maybe because the
     * original formatting was incorrect, as it is in most spectrum
     * generators...). by default it just joins together the strings from
     * blocksAsSingleStrings. each different-scale copy is concatenated. if
     * onlyShowScalesGreaterThanZero is true, blocks with scales of 0.0 or
     * less just do not have the "Q=" etc. printed.
     */
    {
      stringInterpretation.clear();
      for( int scaleIndex( 0 );
           dataBlocks.getSize() > scaleIndex;
           ++scaleIndex )
      {
        stringInterpretation.append(
                          BlockClass::BaseStringBlock::blockIdentifierString );
        stringInterpretation.append( " " );
        stringInterpretation.append( blockName );
        if( onlyShowScalesGreaterThanZero
            ||
            ( 0.0 > dataBlocks[ scaleIndex ].getScale() ) )
        {
          stringInterpretation.append( " Q= " );
          stringInterpretation.append(
                              BlockInterpreter::slhaDoubleMaker.doubleToString(
                                       dataBlocks[ scaleIndex ].getScale() ) );
        }
        stringInterpretation.append( "\n" );
        stringInterpretation.append( dataBlocks[ scaleIndex ].getAsString() );
      }
      return stringInterpretation;
    }

    template< class ValueClass, class BlockParser >
    inline void
    SlhaBlock< ValueClass, BlockParser >::respondToObservedSignal()
    /* by default, PushedToObserver instances don't respond without a pushed
     * value, but this is over-ridden to clear all the data that this block
     * has interpreted or had assigned.
     */
    {
      dataBlocks.setSize( 1 ).getBack().clearEntries();
      scaleOrderedIndices.clear();
      hasHadPushSinceLastReset = false;
    }

    template< class ValueClass, class BlockParser >
    inline void
    SlhaBlock< ValueClass, BlockParser >::respondToPush(
                               BlockClass::BaseStringBlock const& pushedValue )
    // this adds a new BlockParser to dataBlocks & tells it to interpret
    // pushedValue. it also sorts out scale indices.
    {
      if( hasHadPushSinceLastReset )
      {
        dataBlocks.newEnd().setDefaultUnsetValue( defaultUnsetValue );
        dataBlocks.getBack().setVerbosity( isVerbose );
        prepareNewDataBlock();
        dataBlocks.getBack().interpretStringBlock( pushedValue );
        if( pushedValue.getScale()
            < dataBlocks[ lowestScaleIndex ].getScale() )
          /* since hasHadPushSinceLastReset is true, lowestScaleIndex is a
           * valid index for dataBlocks, so the comparison is valid, & if true,
           * lowestScaleIndex is set correctly.
           */
        {
          lowestScaleIndex = dataBlocks.getLastIndex();
        }
      }
      else
      {
        dataBlocks[ 0 ].interpretStringBlock( pushedValue );
      }
      hasHadPushSinceLastReset = true;
      scaleIndexIterator = scaleOrderedIndices.begin();
      while( ( scaleIndexIterator != scaleOrderedIndices.end() )
             &&
             ( scaleIndexIterator->second < pushedValue.getScale() ) )
      {
        ++scaleIndexIterator;
      }
      // now scaleIndexIterator should either be at the index with scale just
      // above blockScale, or at the end of the list.
      scaleOrderedIndices.insert( scaleIndexIterator,
                           std::pair< int, double >( dataBlocks.getLastIndex(),
                                                    pushedValue.getScale() ) );
    }

    template< class ValueClass, class BlockParser >
    inline void
    SlhaBlock< ValueClass, BlockParser >::prepareNewDataBlock()
    // derived classes can insert extra arguments to any new block interpreters
    // by over-riding this.
    {
      // does nothing.
    }

  }

}

#endif /* SLHABLOCK_HPP_ */
