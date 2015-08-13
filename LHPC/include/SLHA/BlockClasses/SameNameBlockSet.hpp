/*
 * SameNameBlockSet.hpp
 *
 *  Created on: Mar 4, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SAMENAMEBLOCKSET_HPP_
#define SAMENAMEBLOCKSET_HPP_

#include <string>
#include <list>
#include "BOLlib/include/BOLlib.hpp"
#include "BaseStringBlock.hpp"

namespace LHPC
{
  namespace SLHA
  {
    typedef
    BOL::PushingObserved< BlockClass::BaseStringBlock > StringBlockPusher;

    /* instances of this class hold together all copies of a block which have
     * the same name (if there are 2 or more copies, it is assumed that they
     * have different energy scale values).
     */
    class SameNameBlockSet : public StringBlockPusher
    {
    public:
      SameNameBlockSet( std::string const& blockName );
      ~SameNameBlockSet();

      std::string const&
      getName() const;
      // this returns the name in uppercase.
      bool
      nameMatches( std::string const& nameToCompare ) const;
      // this returns true if nameToCompare matches blockNameInUppercase
      // ignoring case.
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
      BlockClass::BaseStringBlock &
      operator[]( int whichScaleIndex );
      /* the BaseStringBlocks are indexed in the order in which they were
       * read. the index starts at 1 rather than 0, as well. if 0 is given as
       * the argument, the index corresponding to the copy with the lowest
       * scale is used.
       */
      BlockClass::BaseStringBlock const&
      operator[]( int whichScaleIndex ) const;
      // const version of above.
      int
      getNumberOfCopiesWithDifferentScale() const;
      int
      getLowestScaleIndex() const;
      void
      clearEntries();
      // this clears all the data that this block set has recorded.
      void
      recordHeader( std::string const& headerString,
                    std::string const& commentString,
                    double const blockScale );
      // this prepares a new BaseBlockAsString for the impending block being
      // read as strings.
      void
      recordBodyLine( std::string const& dataString,
                      std::string const& commentString );
      // this sends the lines to the last BaseBlockAsString prepared by
      // recordHeader( ... ).
      void
      finishRecordingLines();
      // this pushes currentStringBlock to all the observers so that they
      // interpret it.


    protected:
      std::string blockNameInUppercase;
      BOL::VectorlikeArray< BlockClass::BaseStringBlock > stringBlocks;
      std::list< std::pair< int, double > > scaleOrderedIndices;
      std::list< std::pair< int, double > >::iterator scaleIndexIterator;
      int lowestScaleIndex;
      BlockClass::BaseStringBlock* currentStringBlock;
    };





    inline std::string const&
    SameNameBlockSet::getName() const
    // this returns the name in uppercase.
    {
      return blockNameInUppercase;
    }

    inline bool
    SameNameBlockSet::nameMatches( std::string const& nameToCompare ) const
    // this returns true if nameToCompare matches blockNameInUppercase
    // ignoring case.
    {
      return BOL::StringParser::stringsMatchIgnoringCase( blockNameInUppercase,
                                                          nameToCompare );
    }

    inline bool
    SameNameBlockSet::hasRecordedScale( double const soughtScale,
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

    inline BlockClass::BaseStringBlock &
    SameNameBlockSet::operator[]( int whichScaleIndex )
    /* the BaseStringBlocks are indexed in the order in which they were
     * read. the index starts at 1 rather than 0, as well. if 0 is given as
     * the argument, the index corresponding to the copy with the lowest
     * scale is used.
     */
    {
      if( 0 == whichScaleIndex )
      {
        return stringBlocks[ lowestScaleIndex ];
      }
      else
      {
        return stringBlocks[ (--whichScaleIndex) ];
      }
    }

    inline BlockClass::BaseStringBlock const&
    SameNameBlockSet::operator[]( int whichScaleIndex ) const
    // const version of above.
    {
      if( 0 == whichScaleIndex )
      {
        return stringBlocks[ lowestScaleIndex ];
      }
      else
      {
        return stringBlocks[ (--whichScaleIndex) ];
      }
    }

    inline int
    SameNameBlockSet::getNumberOfCopiesWithDifferentScale() const
    {
      return stringBlocks.getSize();
    }

    inline int
    SameNameBlockSet::getLowestScaleIndex() const
    {
      return lowestScaleIndex;
    }

    inline void
    SameNameBlockSet::clearEntries()
    // this clears all the data that this block set has recorded.
    {
      stringBlocks.clearEntries();
      scaleOrderedIndices.clear();
      lowestScaleIndex = -1;
      updateObservers();
      // respondToObservedSignal() has been over-ridden for the observers so
      // that this clears their entries.
    }

    inline void
    SameNameBlockSet::recordBodyLine( std::string const& dataString,
                                      std::string const& commentString )
    // this sends the lines to the last BaseBlockAsString prepared by
    // recordHeader( ... ).
    {
      currentStringBlock->recordBodyLine( dataString,
                                          commentString );
    }

    inline void
    SameNameBlockSet::finishRecordingLines()
    // this pushes currentStringBlock to all the observers so that they
    // interpret it.
    {
      updateObservers( *currentStringBlock );
    }

  }

}

#endif /* SAMENAMEBLOCKSET_HPP_ */
