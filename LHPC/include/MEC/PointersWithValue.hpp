/*
 * PointersWithValue.hpp
 *
 *  Created on: Jan 8, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef POINTERSWITHVALUE_HPP_
#define POINTERSWITHVALUE_HPP_

#include <list>
#include <algorithm>
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  /* this is a class to hold pointers (intended to be MassEigenstate pointers)
   * paired with values (intended to be branching ratios as double), with
   * functions for matching.
   */
  template< class PointerClass, class ValueClass >
  class PointersWithValue
  {
  public:
    PointersWithValue();
    PointersWithValue(
             PointersWithValue< PointerClass, ValueClass > const& copySource );
    ~PointersWithValue();

    void
    becomeCopyOf(
             PointersWithValue< PointerClass, ValueClass > const& copySource );
    void
    operator=(
             PointersWithValue< PointerClass, ValueClass > const& copySource ){
                                                  becomeCopyOf( copySource ); }
    PointersWithValue< PointerClass, ValueClass >*
    copyAndReturnPointer(
             PointersWithValue< PointerClass, ValueClass > const& copySource );
    std::list< PointerClass >&
    getPointerList();
    std::list< PointerClass > const&
    getPointerList() const;
    int
    getListSize() const;
    void
    addPointer( PointerClass const pointerToAdd );
    ValueClass
    getPairedValue() const;
    void
    clearPointers();
    void
    setPairedValueAndSortPointers( ValueClass const pairedValue );
    bool
    matchesExactly( PointerClass const firstPointer,
                    PointerClass const secondPointer ) const;
    bool
    matchesExactly( std::list< PointerClass > const& sortedList ) const;
    bool
    matchesExactly( int const firstCode,
                    int const secondCode,
                    bool (*comparisonFunction)( PointerClass const,
                                                int const ) ) const;
    bool
    containsAllAsSubset( std::list< PointerClass > const& sortedList ) const;
    bool
    containsAnyAsSubset( std::list< PointerClass > const& sortedList ) const;


  protected:
    std::list< PointerClass > pointerList;
    ValueClass pairedValue;
  };





  template< class PointerClass, class ValueClass >
  inline
  PointersWithValue< PointerClass, ValueClass >::PointersWithValue() :
      pointerList(),
      pairedValue( BOL::UsefulStuff::notANumber )
  {
    // just an initialization list.
  }

  template< class PointerClass, class ValueClass >
  inline
  PointersWithValue< PointerClass, ValueClass >::PointersWithValue(
            PointersWithValue< PointerClass, ValueClass > const& copySource ) :
      pointerList( copySource.pointerList ),
      pairedValue( copySource.pairedValue )
  {
    // just an initialization list.
  }

  template< class PointerClass, class ValueClass >
  inline
  PointersWithValue< PointerClass, ValueClass >::~PointersWithValue()
  {
    // does nothing.
  }


  template< class PointerClass, class ValueClass >
  inline void
  PointersWithValue< PointerClass, ValueClass >::becomeCopyOf(
             PointersWithValue< PointerClass, ValueClass >  const& copySource )
  {
    this->pointerList = copySource.pointerList;
    this->pairedValue = copySource.pairedValue;
  }

  template< class PointerClass, class ValueClass >
  inline PointersWithValue< PointerClass, ValueClass >*
  PointersWithValue< PointerClass, ValueClass >::copyAndReturnPointer(
           PointersWithValue< PointerClass, ValueClass > const& copySource )
  {
    becomeCopyOf( copySource );
    return this;
  }

  template< class PointerClass, class ValueClass >
  inline std::list< PointerClass >&
  PointersWithValue< PointerClass, ValueClass >::getPointerList()
  {
    return pointerList;
  }

  template< class PointerClass, class ValueClass >
  inline std::list< PointerClass > const&
  PointersWithValue< PointerClass, ValueClass >::getPointerList() const
  {
    return pointerList;
  }

  template< class PointerClass, class ValueClass >
  inline int
  PointersWithValue< PointerClass, ValueClass >::getListSize() const
  {
    return (int)(pointerList.size());
  }

  template< class PointerClass, class ValueClass >
  inline void
  PointersWithValue< PointerClass, ValueClass >::addPointer(
                                              PointerClass const pointerToAdd )
  {
    pointerList.push_back( pointerToAdd );
  }

  template< class PointerClass, class ValueClass >
  inline ValueClass
  PointersWithValue< PointerClass, ValueClass >::getPairedValue() const
  {
    return pairedValue;
  }

  template< class PointerClass, class ValueClass >
  inline void
  PointersWithValue< PointerClass, ValueClass >::clearPointers()
  {
    pointerList.clear();
  }

  template< class PointerClass, class ValueClass >
  inline void
  PointersWithValue< PointerClass, ValueClass >::setPairedValueAndSortPointers(
                                                 ValueClass const pairedValue )
  {
    this->pairedValue = pairedValue;
    pointerList.sort();
  }

  template< class PointerClass, class ValueClass >
  inline bool
  PointersWithValue< PointerClass, ValueClass >::matchesExactly(
                                               PointerClass const firstPointer,
                                       PointerClass const secondPointer ) const
  {
    if( ( 2 == pointerList.size() )
        &&
        ( ( ( firstPointer == pointerList.front() )
            &&
            ( secondPointer == pointerList.back() ) )
          ||
          ( ( secondPointer == pointerList.front() )
            &&
            ( firstPointer == pointerList.back() ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  template< class PointerClass, class ValueClass >
  inline bool
  PointersWithValue< PointerClass, ValueClass >::matchesExactly(
                            std::list< PointerClass > const& sortedList ) const
  {
    if( sortedList.size() == pointerList.size() )
    {
      return std::equal( sortedList.begin(),
                         sortedList.end(),
                         pointerList.begin() );
    }
    else
    {
      return false;
    }
  }

  template< class PointerClass, class ValueClass >
  inline bool
  PointersWithValue< PointerClass, ValueClass >::matchesExactly(
                                                           int const firstCode,
                                                          int const secondCode,
                                bool (*comparisonFunction)( PointerClass const,
                                                            int const ) ) const
  {
    if( ( 2 == pointerList.size() )
        &&
        ( ( (*comparisonFunction)( pointerList.front(),
                                   firstCode )
            &&
            (*comparisonFunction)( pointerList.back(),
                                   secondCode ) )
          ||
          ( (*comparisonFunction)( pointerList.back(),
                                   firstCode )
            &&
            (*comparisonFunction)( pointerList.front(),
                                   secondCode ) ) ) )
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  template< class PointerClass, class ValueClass >
  inline bool
  PointersWithValue< PointerClass, ValueClass >::containsAllAsSubset(
                            std::list< PointerClass > const& sortedList ) const
  {
    bool returnBool( false );
    if( pointerList.size() >= sortedList.size() )
    {
      int allowedMismatches( pointerList.size() - sortedList.size() );
      // the difference in list sizes is the most comparisons that can fail.
      typename std::list< PointerClass >::const_iterator
      fromThisSet( pointerList.begin() );
      typename std::list< PointerClass >::const_iterator
      toBeCompared( sortedList.begin() );
      while( ( allowedMismatches >= 0 )
             &&
             ( fromThisSet != pointerList.end() ) )
      {
        if( *fromThisSet < *toBeCompared )
        {
          --allowedMismatches;
          // note that a comparison failed.
        }
        else
        {
          if( *fromThisSet == *toBeCompared )
          {
            // move on to the next sought pointer, as this one was just found:
            ++toBeCompared;
          }
          else
          {
            allowedMismatches = -1;
            /* otherwise, since the 2 lists are ordered & sortedList is smaller
             * than pointerList, if *thisSetPointer > *comparisonPointer then
             * there's no way that all of sortedList can be in pointerList.
             */
          }
        }
        ++fromThisSet;
      }
      if( allowedMismatches >= 0 )
      {
        returnBool = true;
      }
    }
    return returnBool;
  }

  template< class PointerClass, class ValueClass >
  inline bool
  PointersWithValue< PointerClass, ValueClass >::containsAnyAsSubset(
                            std::list< PointerClass > const& sortedList ) const
  {
    bool returnBool( false );
    typename std::list< PointerClass >::const_iterator
    fromThisSet( pointerList.begin() );
    typename std::list< PointerClass >::const_iterator
    toBeCompared( sortedList.begin() );
    /* both lists are sorted, so we walk each iterator along until it catches
     * up with the other: either it overtakes, & we switch to iterating the
     * other, or they are equal, in which case we should stop looking, and
     * return true.
     */
    while( !returnBool
           &&
           ( fromThisSet != pointerList.end() )
           &&
           ( toBeCompared != sortedList.end() ) )
    {
      if( *fromThisSet < *toBeCompared )
      {
        ++fromThisSet;
      }
      else if( *fromThisSet == *toBeCompared )
      {
        returnBool = true;
      }
      else
      {
        ++toBeCompared;
      }
    }
    return returnBool;
  }

}

#endif /* POINTERSWITHVALUE_HPP_ */
