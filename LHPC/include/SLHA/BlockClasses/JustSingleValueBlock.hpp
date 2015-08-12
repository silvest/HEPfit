/*
 * JustSingleValueBlock.hpp
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

#ifndef JUSTSINGLEVALUEBLOCK_HPP_
#define JUSTSINGLEVALUEBLOCK_HPP_

#include "../SlhaBlock.hpp"
#include "InterpreterClasses/JustSingleValue.hpp"

namespace LHPC
{
  namespace SLHA
  {
    // this template class interprets all the blocks with the same name, though
    // differing scale values, which are interpreted as having a single value.
    template< class ValueClass >
    class JustSingleValueBlock : public SlhaBlock< ValueClass,
                              InterpreterClass::JustSingleValue< ValueClass > >
    {
    public:
      JustSingleValueBlock( std::string const& blockName,
                            ValueClass const& defaultUnsetValue,
                            bool const isVerbose = false );
      virtual
      ~JustSingleValueBlock();

      ValueClass&
      operator()();
      // this returns operator() of the lowest-scale interpreter.
      ValueClass const&
      operator()() const;
      // const version of above.
      bool
      hasEntry() const;
      // this returns hasEntry() of the lowest-scale interpreter.
    };





    template< class ValueClass >
    inline
    JustSingleValueBlock< ValueClass >::JustSingleValueBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                       bool const isVerbose ) :
        SlhaBlock< ValueClass,
                   InterpreterClass::JustSingleValue< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                    isVerbose )
    {
      // just an initialization list.
    }

    template< class ValueClass >
    inline
    JustSingleValueBlock< ValueClass >::~JustSingleValueBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline ValueClass&
    JustSingleValueBlock< ValueClass >::operator()()
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]();
    }

    template< class ValueClass >
    inline ValueClass const&
    JustSingleValueBlock< ValueClass >::operator()() const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]();
    }

    template< class ValueClass >
    inline bool
    JustSingleValueBlock< ValueClass >::hasEntry() const
    // this returns hasEntry() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ].hasEntry();
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace

#endif /* JUSTSINGLEVALUEBLOCK_HPP_ */
