/*
 * LinesAsStrings.cpp
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

#include "SLHA.hpp"

namespace LHPC
{
  namespace SLHA
  {
    namespace InterpreterClass
    {
      LinesAsStrings::LinesAsStrings() :
          InterpreterTemplate< std::string >()
      {
        this->defaultUnsetValue.assign( "" );
      }

      LinesAsStrings::~LinesAsStrings()
      {
        // does nothing.
      }

    }  // end of InterpreterClass namespace

  }  // end of SLHA namespace

}  // end of LHPC namespace
