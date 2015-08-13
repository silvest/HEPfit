/*
 * PushedToObserver.hpp
 *
 *  Created on: Mar 4, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef PUSHEDTOOBSERVER_HPP_
#define PUSHEDTOOBSERVER_HPP_

#include <list>
#include "BasicObserver.hpp"

namespace BOL
{
  // this is a template class for observers which get PushedClass references
  // pushed to them with their observed subject's updates.
  template< typename PushedClass >
  class PushedToObserver : public BasicObserver
  {
  public:
    PushedToObserver();
    virtual
    ~PushedToObserver();

    virtual void
    respondToObservedSignal();
    // by default, PushedToObserver instances don't respond without a pushed
    // value, but this can be over-ridden if it would be useful.
    virtual void
    respondToPush( PushedClass const& pushedValue ) = 0;
  };





  template< typename PushedClass >
  inline
  PushedToObserver< PushedClass >::PushedToObserver() :
      BasicObserver()
  {
    // just an initialization list.
  }

  template< typename PushedClass >
  inline
  PushedToObserver< PushedClass >::~PushedToObserver()
  {
    // does nothing.
  }


  template< typename PushedClass >
  inline void
  PushedToObserver< PushedClass >::respondToObservedSignal()
  // by default, PushedToObserver instances don't respond without a pushed
  // value, but this can be over-ridden if it would be useful.
  {
    // does nothing.
  }

}

#endif /* PUSHEDTOOBSERVER_HPP_ */
