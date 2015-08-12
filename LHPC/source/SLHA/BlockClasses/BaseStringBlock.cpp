/*
 * BaseStringBlock.cpp
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

#include <list>
#include "SLHA.hpp"
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace BlockClass
    {
      std::string const BaseStringBlock::blockIdentifierString( "BLOCK" );
      std::string const BaseStringBlock::decayIdentifierString( "DECAY" );


      BaseStringBlock::BaseStringBlock() :
          blockAsStringWithHeader( "" ),
          stringPairArray(),
          blockScale( BOL::UsefulStuff::notANumber )
      {
        // just an initialization list.
      }

      BaseStringBlock::~BaseStringBlock()
      {
        // does nothing.
      }


      bool
      BaseStringBlock::findScaleIndices( double const soughtScale,
              std::list< std::pair< int, double > > const& scaleOrderedIndices,
                                         int& indexForLowerScale,
                                         int& indexForUpperScale,
                                         double& fractionFromLowerScale )
      /* this looks for the pair of indices with energy scales closest to
       * soughtScale from scaleOrderedIndices (which is assumed to already be
       * ordered by scale).
       * if scaleOrderedIndices is empty, none of the references given are
       * changed & false is returned.
       * if there is only one copy of the block, both indexForLowerScale &
       * indexForUpperScale are set to 1, fractionFromLowerScale is set to
       * NaN, & true is returned.
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
        if( scaleOrderedIndices.empty() )
        {
          return false;
        }
        else
        {
          if( 1 == scaleOrderedIndices.size() )
          {
            indexForLowerScale = 1;
            indexForUpperScale = 1;
            fractionFromLowerScale = BOL::UsefulStuff::notANumber;
          }
          else
            // otherwise, there are 2 or more copies.
          {
            std::list< std::pair< int, double > >::const_iterator
            scaleIterator( scaleOrderedIndices.begin() );
            indexForLowerScale = scaleIterator->first;
            // indexForLowerScale starts at the lowest scale.
            double lowerScale( scaleIterator->second );
            indexForUpperScale = (++scaleIterator)->first;
            // indexForUpperScale starts at the next lowest scale.
            double upperScale( scaleIterator->second );
            while( ( scaleOrderedIndices.end() != scaleIterator )
                   &&
                   ( soughtScale >= scaleIterator->second ) )
              // until scaleIterator goes past soughtScale or hits the end of
              // the list...
            {
              indexForLowerScale = indexForUpperScale;
              indexForUpperScale = (++scaleIterator)->first;
              // both indices are moved up one copy.
              lowerScale = upperScale;
              upperScale = scaleIterator->second;
            }
            /* now if soughtScale was lower than the lowest scale, both indices
             * are still at the lowest pair of scales; if soughtScale was
             * higher than the highest scale, scaleIterator got to the end of
             * the list & the indices now are at the highest pair of scales;
             * otherwise indexForLowerScale is at the highest scale lower than
             * soughtScale & indexForUpperScale is at the lowest scale higher
             * than soughtScale.
             */
            if( 0.0 == ( upperScale - lowerScale ) )
            {
              fractionFromLowerScale = BOL::UsefulStuff::notANumber;
            }
            else
            {
              fractionFromLowerScale
              = ( ( soughtScale - lowerScale ) / ( upperScale - lowerScale ) );
            }
          }
          return true;
        }
      }

    }

  }

}
