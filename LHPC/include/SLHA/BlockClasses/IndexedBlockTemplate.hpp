/*
 * IndexedBlockTemplate.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef INDEXEDBLOCKTEMPLATE_HPP_
#define INDEXEDBLOCKTEMPLATE_HPP_

#include "../SlhaBlock.hpp"

namespace LHPC
{
  namespace SLHA
  {
    // this abstract base template class covers the common functions for a
    // block with at least 1 index.
    template< class ValueClass, class IndexedParser >
    class IndexedBlockTemplate : public SlhaBlock< ValueClass, IndexedParser >
    {
    public:
      IndexedBlockTemplate( std::string const& blockName,
                            ValueClass const& defaultUnsetValue,
                            bool const isVerbose,
                            std::vector< int > const& indexDigitsVector );
      virtual
      ~IndexedBlockTemplate();


    protected:
      std::vector< int > const indexDigitsVector;

      virtual void
      prepareNewDataBlock();
      // derived classes can insert extra arguments to any new block
      // interpreters by over-riding this.
    };





    template< class ValueClass, class IndexedParser >
    inline
    IndexedBlockTemplate< ValueClass, IndexedParser >::IndexedBlockTemplate(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                std::vector< int > const& indexDigitsVector ) :
        SlhaBlock< ValueClass, IndexedParser >( blockName,
                                                defaultUnsetValue,
                                                isVerbose ),
        indexDigitsVector( indexDigitsVector )
    {
      this->dataBlocks.getFront().setIndexDigits( indexDigitsVector );
    }

    template< class ValueClass, class IndexedParser >
    inline
    IndexedBlockTemplate< ValueClass, IndexedParser >::~IndexedBlockTemplate()
    {
      // does nothing.
    }

    template< class ValueClass, class IndexedParser >
    inline void
    IndexedBlockTemplate< ValueClass, IndexedParser >::prepareNewDataBlock()
    {
      this->dataBlocks.getBack().setIndexDigits( indexDigitsVector );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace

#endif /* INDEXEDBLOCKTEMPLATE_HPP_ */
