/*
 * BaseStringBlock.hpp
 *
 *  Created on: Mar 3, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef BASEBLOCKASSTRINGS_HPP_
#define BASEBLOCKASSTRINGS_HPP_

#include <string>
#include <list>
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace BlockClass
    {
      // this class holds a SLHA block as a set of pairs of strings (data with
      // comment).
      class BaseStringBlock
      {
      public:
        static std::string const blockIdentifierString;
        static std::string const decayIdentifierString;

        static bool
        findScaleIndices( double const soughtScale,
              std::list< std::pair< int, double > > const& scaleOrderedIndices,
                          int& indexForLowerScale,
                          int& indexForUpperScale,
                          double& fractionFromLowerScale );
        /* this looks for the pair of indices with energy scales closest to
         * soughtScale from scaleOrderedIndices (which is assumed to already be
         * ordered by scale).
         * if scaleOrderedIndices is empty, none of the references given are
         * changed & false is returned.
         * if there is only one copy of the block, both indexForLowerScale &
         * indexForUpperScale are set to 1, fractionFromLowerScale is set to
         * NaN, & true is returned.
         * if there are 2 or more copies of the block, indexForLowerScale &
         * indexForUpperScale are set as described below,
         * fractionFromLowerScale is set to be
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
         * indexForLowerScale is set to the index of the copy with lowest
         * scale, & indexForUpperScale is set to the index of the copy with the
         * next lowest scale.
         * if soughtScale is higher than the highest scale of the copies,
         * indexForUpperScale is set to the index of the copy with highest
         * scale, & indexForLowerScale is set to the index of the copy with the
         * next highest scale.
         * otherwise, indexForLowerScale is set to the index of the copy with
         * highest scale which is still lower than soughtScale, &
         * indexForUpperScale is set to the index of the copy with lowest scale
         * which is still higher than soughtScale.
         */

        BaseStringBlock();
        ~BaseStringBlock();

        BaseStringBlock*
        recordHeader( std::string const& headerString,
                      std::string const& commentString,
                      double const blockScale );
        // this sets the block to be just a header, so subsequent
        // recordBodyLine(...) calls write the block anew.
        void
        recordBodyLine( std::string const& dataString,
                        std::string const& commentString );
        /* this records dataString & commentString in stringPairArray &
         * then copies dataString into comparisonString, trims it of
         * whitespace, & calls interpretBodyLine() if comparisonString is not
         * then empty.
         */
        int
        getNumberOfBodyLines() const;
        // this returns the number of body lines, so the size of
        // stringPairArray minus 1.
        std::pair< std::string, std::string >&
        operator[]( int const whichLine );
        /* the std::pair< std::string, std::string > at index 0 is the block
         * name (with optional scale & anything else that appeared before the
         * '#') paired with its comment, & the rest of the
         * std::pair< std::string, std::string >s are the data lines paired
         * with their comments as recorded.
         */
        std::pair< std::string, std::string > const&
        operator[]( int const whichLine ) const;
        // const version of above.
        double
        getScale() const;


      protected:
        std::string blockAsStringWithHeader;
        BOL::VectorlikeArray< std::pair< std::string, std::string > >
        stringPairArray;
        /* the std::pair< std::string, std::string > at index 0 is the block
         * name (with optional scale & anything else that appeared before the
         * '#') paired with its comment, & the rest of the
         * std::pair< std::string, std::string >s are the data lines paired
         * with their comments as recorded.
         */
        double blockScale;
      };





      inline BaseStringBlock*
      BaseStringBlock::recordHeader( std::string const& headerString,
                                     std::string const& commentString,
                                     double const blockScale )
      {
        stringPairArray.setSize( 1 );
        stringPairArray.getFront().first.assign( headerString );
        stringPairArray.getFront().second.assign( commentString );
        this->blockScale = blockScale;
        return this;
      }

      inline void
      BaseStringBlock::recordBodyLine( std::string const& dataString,
                                       std::string const& commentString )
      /* this records dataString & commentString in stringPairArray & then
       * copies dataString into comparisonString, trims it of whitespace, &
       * calls interpretBodyLine() if comparisonString is not then empty.
       */
      {
        stringPairArray.newEnd().first.assign( dataString );
        stringPairArray.getBack().second.assign( commentString );
      }

      inline std::pair< std::string, std::string >&
      BaseStringBlock::operator[]( int const whichLine )
      /* the std::pair< std::string, std::string > at index 0 is the block
       * name (with optional scale & anything else that appeared before the
       * '#') paired with its comment, & the rest of the
       * std::pair< std::string, std::string >s are the data lines paired
       * with their comments as recorded.
       */
      {
        return stringPairArray[ whichLine ];
      }

      inline std::pair< std::string, std::string > const&
      BaseStringBlock::operator[]( int const whichLine ) const
      // const version of above.
      {
        return stringPairArray[ whichLine ];
      }

      inline int
      BaseStringBlock::getNumberOfBodyLines() const
      // this returns the number of body lines, so the size of
      // stringPairArray minus 1.
      {
        return stringPairArray.getLastIndex();
      }

      inline double
      BaseStringBlock::getScale() const
      {
        return blockScale;
      }

    }

  }

}

#endif /* BASEBLOCKASSTRINGS_HPP_ */
