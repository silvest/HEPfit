/*
 * PushingObserved.hpp
 *
 *  Created on: Mar 4, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef PUSHINGOBSERVED_HPP_
#define PUSHINGOBSERVED_HPP_

#include <list>
#include "PushedToObserver.hpp"

namespace BOL
{
  /* this template class holds a list of PushedToObserver< PushedClass >
   * pointers & calls respondToObservedSignal( pushedValue ) on them all with
   * updateObservers( pushedValue ).
   */
  template< typename PushedClass >
  class PushingObserved
  {
  public:
    PushingObserved();
    virtual
    ~PushingObserved();

    virtual void
    updateObservers( PushedClass const& pushedValue );
    /* this goes through observerList & either removes the observer if it has
     * flagged its bool as false, or calls its respondToPush( pushedValue )
     * otherwise.
     */
    virtual void
    updateObservers();
    /* this goes through observerList & either removes the observer if it has
     * flagged its bool as false, or calls its respondToObservedSignal()
     * otherwise.
     */
    virtual void
    registerObserver( PushedToObserver< PushedClass >* const joiningObserver );
    virtual void
    removeObserver( PushedToObserver< PushedClass >* const leavingObserver );
    /* this goes through observerList & removes any observers which have
     * changed their flags to false, & also removes leavingObserver if found
     * (1st asking it to discard its bool pointer).
     */
    void
    removeAllObservers();
    // this clears observerList after asking all its observers to discard their
    // bool pointers for this PushingObserved instance.


  protected:
    typedef typename
    std::pair< PushedToObserver< PushedClass >*, bool > observerWithBool;
    typedef typename
    std::list< observerWithBool >::iterator observerWithBoolListIterator;
    std::list< observerWithBool > observerList;
    observerWithBoolListIterator observerIterator;
  };





  template< typename PushedClass >
  inline
  PushingObserved< PushedClass >::PushingObserved() :
      observerList(),
      observerIterator()
  {
    // just an initialization list.
  }

  template< typename PushedClass >
  inline
  PushingObserved< PushedClass >::~PushingObserved()
  {
    removeAllObservers();
  }


  template< typename PushedClass >
  inline void
  PushingObserved< PushedClass >::updateObservers(
                                               PushedClass const& pushedValue )
  /* this goes through observerList & either removes the observer if it has
   * flagged its bool as false, or calls its respondToPush( pushedValue )
   * otherwise.
   */
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( observerIterator->second )
      {
        observerIterator->first->respondToPush( pushedValue );
        ++observerIterator;
      }
      else
      {
        observerIterator = observerList.erase( observerIterator );
      }
    }
  }

  template< typename PushedClass >
  inline void
  PushingObserved< PushedClass >::updateObservers()
  /* this goes through observerList & either removes the observer if it has
   * flagged its bool as false, or calls its respondToObservedSignal()
   * otherwise.
   */
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( observerIterator->second )
      {
        observerIterator->first->respondToObservedSignal();
        ++observerIterator;
      }
      else
      {
        observerIterator = observerList.erase( observerIterator );
      }
    }
  }

  template< typename PushedClass >
  inline void
  PushingObserved< PushedClass >::registerObserver(
                       PushedToObserver< PushedClass >* const joiningObserver )
  {
    observerList.push_back( observerWithBool( joiningObserver,
                                              true ) );
    joiningObserver->acceptFlagFromObserved( &(observerList.back().second) );
    /* new observers are given a pointer to their associated bools, which
     * indicate that their observers should be fine for
     * respondToPush( ... ) calls if the bool is true, but if the bool is
     * false, this PushingObserved instance knows to remove the observer with
     * its next pass through the list.
     */
  }

  template< typename PushedClass >
  inline void
  PushingObserved< PushedClass >::removeObserver(
                       PushedToObserver< PushedClass >* const leavingObserver )
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( !(observerIterator->second) )
      {
        observerIterator = observerList.erase( observerIterator );
      }
      else if( leavingObserver == observerIterator->first )
      {
        leavingObserver->discardFlagFromObserved(
                                                 &(observerIterator->second) );
        observerIterator = observerList.erase( observerIterator );
      }
      else
      {
        ++observerIterator;
      }
    }
  }

  template< typename PushedClass >
  inline void
  PushingObserved< PushedClass >::removeAllObservers()
  // this clears observerList after asking all its observers to discard their
  // bool pointers for this BasicObserved instance.
  {
    observerIterator = observerList.begin();
    while( observerList.end() != observerIterator )
    {
      if( observerIterator->second )
      {
        observerIterator->first->discardFlagFromObserved(
                                                 &(observerIterator->second) );
      }
      ++observerIterator;
    }
    observerList.clear();
  }

}

#endif /* PUSHINGOBSERVED_HPP_ */
