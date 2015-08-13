/*
 * BasicObserved.cpp
 *
 *  Created on: Mar 4, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#include "BasicObserved.hpp"

namespace BOL
{
  BasicObserved::BasicObserved() :
      observerList(),
      observerIterator()
  {
    // just an initialization list.
  }

  BasicObserved::~BasicObserved()
  {
    removeAllObservers();
  }

}
