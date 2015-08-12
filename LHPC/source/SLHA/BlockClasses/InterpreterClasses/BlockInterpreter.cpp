/*
 * BlockInterpreter.cpp
 *
 *  Created on: Mar 11, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "SLHA.hpp"
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  namespace SLHA
  {
    BOL::StringParser const BlockInterpreter::slhaDoubleMaker( 9,
                                                               ' ',
                                                               9,
                                                               3 );
    BOL::StringParser const BlockInterpreter::slhaIntHelper( 1,
                                                             ' ',
                                                             1,
                                                             1,
                                                             "" );
    BOL::StringParser const BlockInterpreter::particleCodeMaker( 9,
                                                                 ' ',
                                                                 1,
                                                                 1,
                                                                 "" );
    bool const BlockInterpreter::defaultVerbosity( false );


    BlockInterpreter::BlockInterpreter() :
        currentStringBlock( NULL ),
        stringInterpretation( "" )
    {
      // just an initialization list.
    }

    BlockInterpreter::~BlockInterpreter()
    {
      // does nothing.
    }

  }

}
