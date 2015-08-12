/*
 * LinesAsStringsBlock.hpp
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

#ifndef LINESASSTRINGSBLOCK_HPP_
#define LINESASSTRINGSBLOCK_HPP_

#include <string>
#include "../SlhaBlock.hpp"
#include "InterpreterClasses/LinesAsStrings.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted just as a set of lines as
     * std::strings.
     */
    class LinesAsStringsBlock : public SlhaBlock< std::string,
                                             InterpreterClass::LinesAsStrings >
    {
    public:
      LinesAsStringsBlock( std::string const& blockName,
                           bool const isVerbose = false );
      virtual
      ~LinesAsStringsBlock();

      std::string
      operator()( int const whichLine ) const;
      // this returns operator() of the lowest-scale interpreter.
    };





    inline std::string
    LinesAsStringsBlock::operator()( int const whichLine ) const
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( whichLine );
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace


#endif /* LINESASSTRINGSBLOCK_HPP_ */
