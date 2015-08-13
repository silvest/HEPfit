/*
 * VectorlikeArray.hpp
 *
 *  Created on: Jan 4, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef VECTORLIKEARRAY_HPP_
#define VECTORLIKEARRAY_HPP_

#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <list>

namespace BOL
{
  /* this is a template class to do the job of a container similar to a
   * std::vector, but by acting like an array with a number of elements
   * attached, which can be re-sized. a std::vector is usually going to be a
   * better choice, but there might be circumstances where a VectorlikeArray is
   * preferable.
   */
  template< typename StoredClass >
  class VectorlikeArray
  {
  public:
    VectorlikeArray( int const initialSize = 0 );
    VectorlikeArray( VectorlikeArray< StoredClass > const& copySource,
        StoredClass* (*storedClassCopyFunction)( StoredClass const& ) = NULL );
    /* the copy constructor requires a function that returns a pointer to a
     * StoredClass instance that was made with new, since the VectorlikeArray
     * instance will be responsible for freeing its memory with delete. (the
     * function should of course make the new StoredClass instance as a copy
     * of the StoredClass instance that is passed as its argument.) if the
     * function pointer for copying given is NULL, the StoredClass copies are
     * attempted to be initialized with the default StoredClass copy
     * constructor (i.e. with new StoredClass( StoredClass const& ) with the
     * appropriate argument).
     */
    ~VectorlikeArray();

    StoredClass&
    operator[]( int const soughtIndex );
    // this returns a reference to the requested element. it does not check
    // that the element is in the range!
    StoredClass const&
    operator[]( int const soughtIndex ) const;
    // this returns a reference to the requested element. it does not check
    // that the element is in the range!
    bool
    isEmpty() const;
    // this returns false if the number of elements which are deemed current is
    // greater than zero, true otherwise.
    int
    getSize() const;
    // this returns the number of elements which are deemed current.
    int
    getLastIndex() const;
    // this returns the number for accessing the last of elements which are
    // deemed current.
    VectorlikeArray< StoredClass >&
    setSize( int const newSize );
    // this sets the size of the currently-available container, ensuring that
    // enough elements exist.
    VectorlikeArray< StoredClass >&
    increaseSize( int const numberToAdd = 1 );
    // this increases the size of the currently-available container, ensuring
    // that enough elements exist.
    StoredClass&
    newEnd();
    // this increases the size of the currently-available container by 1,
    // ensuring that enough elements exist, & returns the last element.
    VectorlikeArray< StoredClass >&
    clearEntries();
    // this is just setSize( 0 ).
    StoredClass*
    getPointer( int const soughtIndex );
    // if an element out of range is asked for, this returns NULL.
    StoredClass const*
    getPointer( int const soughtIndex ) const;
    // if an element out of range is asked for, this returns NULL.
    StoredClass&
    safeAt( int const soughtIndex );
    // if an element out of range is asked for, this throws an exception.
    StoredClass const&
    safeAt( int const soughtIndex ) const;
    // if an element out of range is asked for, this throws an exception.
    StoredClass&
    getFront();
    // this returns the 1st element deemed current, throwing an out_of_range
    // exception if the container is empty.
    StoredClass const&
    getFront() const;
    // this returns the 1st element deemed current, throwing an out_of_range
    // exception if the container is empty.
    StoredClass&
    getBack();
    // this returns the last element deemed current, throwing an out_of_range
    // exception if the container is empty.
    StoredClass const&
    getBack() const;
    // this returns the last element deemed current, throwing an out_of_range
    // exception if the container is empty.
    std::vector< StoredClass* >&
    getAsVector( std::vector< StoredClass* >& vectorToFill ) const;
    // this fills the given vector with pointers to the current elements,
    // returning a reference to the vector.
    std::list< StoredClass* >&
    getAsList( std::list< StoredClass* >& listToFill ) const;
    // this fills the given list with pointers to the current elements,
    // returning a reference to the list.
    VectorlikeArray< StoredClass >&
    memoryFreeingResize( int const newSize = 0 );
    // this frees up as much memory as possible held in this VectorlikeArray by
    // leaving it as big as specified by newSize.
    VectorlikeArray< StoredClass >&
    becomeCopyOf( VectorlikeArray< StoredClass > const& copySource,
                  void (*firstArgumentBecomesCopyOfSecond)( StoredClass&,
                                                 StoredClass const& ) = NULL );
    /* this sets this instance of VectorlikeArray to be a copy of copySource.
     * if the function pointer that is meant to point to a function that sets a
     * StoredClass instance to be a copy of another is given as NULL, the
     * StoredClass copies are attempted to be copied with the default
     * StoredClass operator= (i.e. with newCopy = copiedInstance
     * appropriately). it returns a reference to this instance.
     */


  protected:
    StoredClass** pointerArray;
    // this is the pointer to the array of pointers to the held StoredClass
    // instances.
    int allocatedSizeOfArray;
    // this is the number of pointers in the array in memory.
    int indexOfEndConstructed;
    // this is the index of the last pointer in the array which points to a
    // constructed StoredClass instance.
    int currentSizeOfArray;
    // this is the number of pointers as the user should see it.
    int currentEndIndex;
    // this is the index of the last pointer in the array as the user should
    // see it.

    StoredClass*
    unsafeGetPointer( int const soughtIndex );
    // this does not check bounds!
    StoredClass const*
    unsafeGetPointer( int const soughtIndex )
    const;
    // this does not check bounds!
    StoredClass*
    safeGetPointer( int const soughtIndex );
    // if an element out of range is asked for, this returns NULL.
    StoredClass const*
    safeGetPointer( int const soughtIndex ) const;
    // if an element out of range is asked for, this returns NULL.
    StoredClass&
    getReference( int const soughtIndex );
    // if an element out of range is asked for, this throws an exception.
    StoredClass const&
    getReference( int const soughtIndex )
    const;
    // if an element out of range is asked for, this throws an exception.
  };



  template< typename StoredClass >
  inline
  VectorlikeArray< StoredClass >::VectorlikeArray( int const initialSize ) :
      pointerArray( new StoredClass*[ initialSize ] ),
      allocatedSizeOfArray( initialSize ),
      indexOfEndConstructed( -1 ),
      currentSizeOfArray( initialSize ),
      currentEndIndex( initialSize - 1 )
  {
    while( currentEndIndex > indexOfEndConstructed )
    {
      ++indexOfEndConstructed;
      pointerArray[ indexOfEndConstructed ] = new StoredClass;
    }
  }

  /* the copy constructor requires a function that returns a pointer to a
   * StoredClass instance that was made with new, since the VectorlikeArray
   * instance will be responsible for freeing its memory with delete. (the
   * function should of course make the new StoredClass instance as a copy
   * of the StoredClass instance that is passed as its argument.) if the
   * function pointer for copying given is NULL, the StoredClass copies are
   * attempted to be initialized with the default StoredClass copy
   * constructor (i.e. with new StoredClass( StoredClass const& ) with the
   * appropriate argument).
   */
  template< typename StoredClass >
  inline
  VectorlikeArray< StoredClass >::VectorlikeArray(
                              VectorlikeArray< StoredClass > const& copySource,
              StoredClass* (*storedClassCopyFunction)( StoredClass const& ) ) :
      pointerArray( new StoredClass*[ copySource.currentSizeOfArray ] ),
      allocatedSizeOfArray( copySource.currentSizeOfArray ),
      indexOfEndConstructed( -1 ),
      currentSizeOfArray( copySource.currentSizeOfArray ),
      currentEndIndex( copySource.currentEndIndex )
  {
    if( NULL == storedClassCopyFunction )
    {
      while( currentEndIndex > indexOfEndConstructed )
      {
        ++indexOfEndConstructed;
        pointerArray[ indexOfEndConstructed ]
        = new StoredClass( copySource[ indexOfEndConstructed ] );
      }
    }
    else
    {
      while( currentEndIndex > indexOfEndConstructed )
      {
        ++indexOfEndConstructed;
        pointerArray[ indexOfEndConstructed ]
        = (*storedClassCopyFunction)( copySource[ indexOfEndConstructed ] );
      }
    }
  }

  template< typename StoredClass >
  inline
  VectorlikeArray< StoredClass >::~VectorlikeArray()
  {
    while( 0 <= indexOfEndConstructed )
    {
      delete pointerArray[ indexOfEndConstructed ];
      --indexOfEndConstructed;
    }
    delete[] pointerArray;
  }


  template< typename StoredClass >
  inline StoredClass&
  VectorlikeArray< StoredClass >::operator[]( int const soughtIndex )
  // this returns a reference to the requested element. it does not check that
  // the element is in the range!
  {
    return *(unsafeGetPointer( soughtIndex ));
  }

  template< typename StoredClass >
  inline StoredClass const&
  VectorlikeArray< StoredClass >::operator[]( int const soughtIndex ) const
  // this returns a reference to the requested element. it does not check that
  // the element is in the range!
  {
    return *(unsafeGetPointer( soughtIndex ));
  }

  template< typename StoredClass >
  inline bool
  VectorlikeArray< StoredClass >::isEmpty() const
  // this returns false if the number of elements which are deemed current is
  // greater than zero, true otherwise.
  {
    if( 0 < currentSizeOfArray )
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  template< typename StoredClass >
  inline int
  VectorlikeArray< StoredClass >::getSize()
  const
  // this returns the number of elements which are deemed current.
  {
    return currentSizeOfArray;
  }

  template< typename StoredClass >
  inline int
  VectorlikeArray< StoredClass >::getLastIndex() const
  // this returns the number for getPointer to access the last of elements
  // which are deemed current.
  {
    return currentEndIndex;
  }

  template< typename StoredClass >
  inline VectorlikeArray< StoredClass >&
  VectorlikeArray< StoredClass >::setSize( int const newSize )
  // this sets the size of the currently-available container, ensuring that
  // enough elements exist.
  {
    // first the conceptual end of the array is moved:
    currentSizeOfArray = newSize;
    currentEndIndex = ( newSize - 1 );

    // then the array of pointers is ensured to be large enough:
    if( currentSizeOfArray > allocatedSizeOfArray )
    {
      if( !(0 < allocatedSizeOfArray ) )
      {
        allocatedSizeOfArray = currentSizeOfArray;
      }
      else
      {
        while( currentSizeOfArray > allocatedSizeOfArray )
        {
          allocatedSizeOfArray += allocatedSizeOfArray;
          // the array size is doubled if it gets increased, following the
          // logic of the std::vector class.
        }
      }
      // now allocatedSizeOfArray is large enough.

      StoredClass**
      newPointerArray = new StoredClass*[ allocatedSizeOfArray ];
      // now newCurrentPointersArray is a new, larger array where the pointers
      // can be put.

      for( int copyIndex( indexOfEndConstructed );
           0 <= copyIndex;
           --copyIndex )
      {
        newPointerArray[ copyIndex ] = pointerArray[ copyIndex ];
      }
      // now newCurrentPointersArray has copied all the contents of
      // pointerArray, so it is OK to delete[] pointerArray:
      delete[] pointerArray;
      pointerArray = newPointerArray;
    }

    // next it is ensured that enough pointers point to constructed instances:
    while( currentEndIndex > indexOfEndConstructed )
    {
      ++indexOfEndConstructed;
      pointerArray[ indexOfEndConstructed ] = new StoredClass;
    }

    return *this;
  }

  template< typename StoredClass >
  inline VectorlikeArray< StoredClass >&
  VectorlikeArray< StoredClass >::increaseSize( int const numberToAdd )
  // this increases the size of the currently-available container, ensuring
  // that enough elements exist.
  {
    return setSize( numberToAdd + currentSizeOfArray );
  }

  template< typename StoredClass >
  inline StoredClass&
  VectorlikeArray< StoredClass >::newEnd()
  // this increases the size of the currently-available container by 1,
  // ensuring that enough elements exist, & returns the last element.
  {
    increaseSize( 1 );
    return *(unsafeGetPointer( currentEndIndex ));
  }

  template< typename StoredClass >
  inline VectorlikeArray< StoredClass >&
  VectorlikeArray< StoredClass >::clearEntries()
  // this is just setSize( 0 ).
  {
    return setSize( 0 );
  }

  template< typename StoredClass >
  inline StoredClass*
  VectorlikeArray< StoredClass >::getPointer( int const soughtIndex )
  // if an element out of range is asked for, this returns NULL.
  {
    return safeGetPointer( soughtIndex );
  }

  template< typename StoredClass >
  inline StoredClass const*
  VectorlikeArray< StoredClass >::getPointer( int const soughtIndex ) const
  // if an element out of range is asked for, this returns NULL.
  {
    return safeGetPointer( soughtIndex );
  }

  template< typename StoredClass >
  inline StoredClass&
  VectorlikeArray< StoredClass >::safeAt( int const soughtIndex )
  // if an element out of range is asked for, this throws an exception.
  {
    return getReference( soughtIndex );
  }

  template< typename StoredClass >
  inline StoredClass const&
  VectorlikeArray< StoredClass >::safeAt( int const soughtIndex ) const
  // if an element out of range is asked for, this throws an exception.
  {
    return getReference( soughtIndex );
  }

  template< typename StoredClass >
  inline StoredClass&
  VectorlikeArray< StoredClass >::getReference( int const soughtIndex )
  // if an element out of range is asked for, this throws an exception.
  {
    StoredClass* soughtPointer( safeGetPointer( soughtIndex ) );
    if( NULL == soughtPointer )
    {
      throw std::out_of_range(
                           "VectorlikeArray::getReference(...) out of range" );
    }
    return *soughtPointer;
  }

  template< typename StoredClass >
  inline StoredClass const&
  VectorlikeArray< StoredClass >::getReference( int const soughtIndex ) const
  // if an element out of range is asked for, this throws an exception.
  {
    StoredClass* soughtPointer( safeGetPointer( soughtIndex ) );
    if( NULL == soughtPointer )
    {
      throw std::out_of_range(
                           "VectorlikeArray::getReference(...) out of range" );
    }
    return *soughtPointer;
  }

  template< typename StoredClass >
  inline StoredClass&
  VectorlikeArray< StoredClass >::getFront()
  // this returns the 1st element deemed current, throwing an out_of_range
  // exception if the container is empty.
  {
    return *(safeGetPointer( 0 ));
  }

  template< typename StoredClass >
  inline StoredClass const&
  VectorlikeArray< StoredClass >::getFront() const
  // this returns the 1st element deemed current, throwing an out_of_range
  // exception if the container is empty.
  {
    return *(safeGetPointer( 0 ));
  }

  template< typename StoredClass >
  inline StoredClass&
  VectorlikeArray< StoredClass >::getBack()
  // this returns the last element deemed current, throwing an out_of_range
  // exception if the container is empty.
  {
    return *(safeGetPointer( currentEndIndex ));
  }

  template< typename StoredClass >
  inline StoredClass const&
  VectorlikeArray< StoredClass >::getBack() const
  // this returns the last element deemed current, throwing an out_of_range
  // exception if the container is empty.
  {
    return *(safeGetPointer( currentEndIndex ));
  }

  template< typename StoredClass >
  inline std::vector< StoredClass* >&
  VectorlikeArray< StoredClass >::getAsVector(
                              std::vector< StoredClass* >& vectorToFill ) const
  // this fills the given vector with pointers to the current elements,
  // returning a reference to the vector.
  {
    vectorToFill.assign( pointerArray,
                         ( pointerArray + currentSizeOfArray ) );
    return vectorToFill;
  }

  template< typename StoredClass >
  inline std::list< StoredClass* >&
  VectorlikeArray< StoredClass >::getAsList(
                                  std::list< StoredClass* >& listToFill ) const
  // this fills the given list with pointers to the current elements,
  // returning a reference to the list.
  {
    listToFill.assign( pointerArray,
                       ( pointerArray + currentSizeOfArray ) );
    return listToFill;
  }

  template< typename StoredClass >
  inline VectorlikeArray< StoredClass >&
  VectorlikeArray< StoredClass >::memoryFreeingResize( int const newSize )
  // this frees up as much memory as possible held in this VectorlikeArray by
  // leaving it as big as specified by newSize.
  {
    currentSizeOfArray = newSize;
    currentEndIndex = ( newSize - 1 );
    while( indexOfEndConstructed >= newSize )
    {
      delete pointerArray[ indexOfEndConstructed ];
      --indexOfEndConstructed;
    }
    if( newSize != allocatedSizeOfArray )
    {
      allocatedSizeOfArray = newSize;
      StoredClass**
      newPointerArray = new StoredClass*[ allocatedSizeOfArray ];
      // now newCurrentPointersArray is a new array where the pointers
      // can be put.

      for( int copyIndex( indexOfEndConstructed );
           0 <= copyIndex;
           --copyIndex )
      {
        newPointerArray[ copyIndex ] = pointerArray[ copyIndex ];
      }
      // now newCurrentPointersArray has copied all the contents of
      // pointerArray, so it is OK to delete[] pointerArray:
      delete[] pointerArray;
      pointerArray = newPointerArray;
    }

    // next it is ensured that enough pointers point to constructed instances:
    while( currentEndIndex > indexOfEndConstructed )
    {
      pointerArray[ ++indexOfEndConstructed ] = new StoredClass;
      // the index is increased then the pointer at the incremented index is
      // set to point at a new StoredClass instance.
    }

    return *this;
  }

  template< typename StoredClass >
  inline StoredClass*
  VectorlikeArray< StoredClass >::unsafeGetPointer( int const soughtIndex )
  // this does not check bounds!
  {
    return pointerArray[ soughtIndex ];
  }

  template< typename StoredClass >
  inline StoredClass const*
  VectorlikeArray< StoredClass >::unsafeGetPointer(
                                                  int const soughtIndex ) const
  // this does not check bounds!
  {
    return pointerArray[ soughtIndex ];
  }

  template< typename StoredClass >
  inline StoredClass*
  VectorlikeArray< StoredClass >::safeGetPointer( int const soughtIndex )
  // if an element out of range is asked for, this returns NULL.
  {
    if( ( 0 <= soughtIndex )
        &&
        ( currentEndIndex >= soughtIndex ) )
    {
      return unsafeGetPointer( soughtIndex );
    }
    else
    {
      return NULL;
    }
  }

  template< typename StoredClass >
  inline StoredClass const*
  VectorlikeArray< StoredClass >::safeGetPointer( int const soughtIndex ) const
  // if an element out of range is asked for, this returns NULL.
  {
    if( ( 0 <= soughtIndex )
        &&
        ( currentEndIndex >= soughtIndex ) )
    {
      return unsafeGetPointer( soughtIndex );
    }
    else
    {
      return NULL;
    }
  }

  template< typename StoredClass >
  inline VectorlikeArray< StoredClass >&
  VectorlikeArray< StoredClass >::becomeCopyOf(
                              VectorlikeArray< StoredClass > const& copySource,
                        void (*firstArgumentBecomesCopyOfSecond)( StoredClass&,
                                                         StoredClass const& ) )
  /* this sets this instance of VectorlikeArray to be a copy of copySource. if
   * the function pointer that is meant to point to a function that sets a
   * StoredClass instance to be a copy of another is given as NULL, the
   * StoredClass copies are attempted to be copied with the default StoredClass
   * operator= (i.e. with newCopy = copiedInstance appropriately). it returns
   * a reference to this instance.
   */
  {
    setSize( copySource.getSize() );
    if( NULL == firstArgumentBecomesCopyOfSecond )
    {
      for( int copyIndex( copySource.getLastIndex() );
           0 <= copyIndex;
           --copyIndex )
      {
        *(pointerArray[ copyIndex ]) = copySource[ copyIndex ];
      }
    }
    else
    {
      for( int copyIndex( copySource.getLastIndex() );
           0 <= copyIndex;
           --copyIndex )
      {
        (*firstArgumentBecomesCopyOfSecond)( *(pointerArray[ copyIndex ]),
                                             copySource[ copyIndex ] );
      }
    }
    return *this;
  }

}

#endif /* VECTORLIKEARRAY_HPP_ */
