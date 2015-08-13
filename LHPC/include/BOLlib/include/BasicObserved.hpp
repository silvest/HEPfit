/*
 * BasicObserved.hpp
 *
 *  Created on: Mar 4, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef BASICOBSERVED_HPP_
#define BASICOBSERVED_HPP_

#include <list>
#include "BasicObserver.hpp"

namespace BOL
{
  // this class holds a list of BasicObserver pointers & calls
  // respondToObservedSignal() on them all with updateObservers().
  class BasicObserved
  {
  public:
    BasicObserved();
    virtual
    ~BasicObserved();

    void
    updateObservers();
    /* this goes through observerList & either removes the observer if it has
     * flagged its bool as false, or calls its respondToObservedSignal()
     * otherwise.
     */
    void
    registerObserver( BasicObserver* const joiningObserver );
    void
    removeObserver( BasicObserver* const leavingObserver );
    /* this goes through observerList & removes any observers which have changed
     * their flags to false, & also removes leavingObserver if found (1st asking
     * it to discard its bool pointer).
     */
    void
    removeAllObservers();
    // this clears observerList after asking all its observers to discard their
    // bool pointers for this BasicObserved instance.


  protected:
    std::list< std::pair< BasicObserver*, bool > > observerList;
    std::list< std::pair< BasicObserver*, bool > >::iterator observerIterator;
  };





  inline void
  BasicObserved::updateObservers()
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

  inline void
  BasicObserved::registerObserver( BasicObserver* const joiningObserver )
  {
    observerList.push_back( std::pair< BasicObserver*, bool >( joiningObserver,
                                                               true ) );
    joiningObserver->acceptFlagFromObserved( &(observerList.back().second) );
    /* new observers are given a pointer to their associated bools, which
     * indicate that their observers should be fine for
     * respondToObservedSignal() calls if the bool is true, but if the bool is
     * false, this BasicObserved instance knows to remove the observer with
     * its next pass through the list.
     */
  }

  inline void
  BasicObserved::removeObserver( BasicObserver* const leavingObserver )
  /* this goes through observerList & removes any observers which have changed
   * their flags to false, & also removes leavingObserver if found (1st asking
   * it to discard its bool pointer).
   */
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

  inline void
  BasicObserved::removeAllObservers()
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

#endif /* BASICOBSERVED_HPP_ */
