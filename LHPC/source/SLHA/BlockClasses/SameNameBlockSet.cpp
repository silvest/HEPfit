/*
 * SameNameBlockSet.cpp
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

#include "SLHA.hpp"

namespace LHPC
{
  namespace SLHA
  {
    SameNameBlockSet::SameNameBlockSet( std::string const& blockName ) :
        StringBlockPusher(),
        blockNameInUppercase( BOL::StringParser::trimFromFrontAndBack(
                                                                     blockName,
                              BOL::StringParser::whitespaceAndNewlineChars ) ),
        stringBlocks(),
        scaleOrderedIndices(),
        scaleIndexIterator(),
        lowestScaleIndex( -1 ),
        currentStringBlock( NULL )
    {
      BOL::StringParser::transformToUppercase( blockNameInUppercase );
    }

    SameNameBlockSet::~SameNameBlockSet()
    {
      // does nothing.
    }


    void
    SameNameBlockSet::recordHeader( std::string const& headerString,
                                    std::string const& commentString,
                                    double const blockScale )
    // this prepares a new BaseBlockAsString for the impending block being
    // read as strings.
    {
      currentStringBlock = stringBlocks.newEnd().recordHeader( headerString,
                                                               commentString,
                                                               blockScale );
      if( 1 == stringBlocks.getSize() )
        // if this is the only copy so far...
      {
        lowestScaleIndex = 0;
      }
      else if( blockScale < stringBlocks[ lowestScaleIndex ].getScale() )
        /* otherwise, defaultDataBlockIndex is a valid index for stringBlocks,
         * so the comparison is valid, & if true, defaultDataBlockIndex is set
         * correctly.
         */
      {
        lowestScaleIndex = stringBlocks.getLastIndex();
      }
      scaleIndexIterator = scaleOrderedIndices.begin();
      while( ( scaleIndexIterator != scaleOrderedIndices.end() )
             &&
             ( scaleIndexIterator->second < blockScale ) )
      {
        ++scaleIndexIterator;
      }
      // now scaleIndexIterator should either be at the index with scale just
      // above blockScale, or at the end of the list.
      scaleOrderedIndices.insert( scaleIndexIterator,
                         std::pair< int, double >( stringBlocks.getLastIndex(),
                                                   blockScale ) );
    }

  }

}
