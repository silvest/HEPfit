/*
 * MapAndVectorAndBools.hpp
 *
 *  Created on: Jan 19, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef MAPANDVECTORANDBOOLS_HPP_
#define MAPANDVECTORANDBOOLS_HPP_

#include <map>
#include <vector>

namespace LHPC
{
  /* this is a class to hold a std::map of ints and pointers (intended to be
   * MassEigenstate pointers) along with a std::vector of that type of pointer,
   * and also a pointer to a std::vector< bool > (intended to hold whatever
   * flags the MassEigenstate instances are meant to have) and a const
   * reference to a bool (intended for a verbosity flag).
   */
  template< class PointerClass >
  class MapAndVectorAndBools
  {
  public:
    typedef typename std::map< int,
                               PointerClass > PointerMap;
    typedef typename std::vector< PointerClass > PointerVector;

    MapAndVectorAndBools( PointerMap& codeMap,
                          PointerVector& pointerGroup,
                          bool const boolReference );
    virtual
    ~MapAndVectorAndBools();

    MapAndVectorAndBools< PointerClass >&
    withBools( std::vector< bool > const* const boolVector );
    // this sets this->boolVector to the given pointer & returns *this.
    PointerMap&
    getMap();
    PointerVector&
    getVector();
    bool
    getBool();
    std::vector< bool > const*
    getFlags();


  protected:
    PointerMap& codeMap;
    PointerVector& pointerGroup;
    bool const boolReference;
    std::vector< bool > const* boolVector;
  };



  template< class PointerClass >
  inline
  MapAndVectorAndBools< PointerClass >::MapAndVectorAndBools(
                                                           PointerMap& codeMap,
                                                   PointerVector& pointerGroup,
                                                   bool const boolReference ) :
      codeMap( codeMap ),
      pointerGroup( pointerGroup ),
      boolReference( boolReference ),
      boolVector( NULL )
  {
    // just an initialization list.
  }

  template< class PointerClass >
  inline
  MapAndVectorAndBools< PointerClass >::~MapAndVectorAndBools()
  {
    // does nothing.
  }


  template< class PointerClass >
  inline MapAndVectorAndBools< PointerClass >&
  MapAndVectorAndBools< PointerClass >::withBools(
                                  std::vector< bool > const* const boolVector )
  // this sets this->boolVector to the given pointer & returns *this.
  {
    this->boolVector = boolVector;
    return *this;
  }

  template< class PointerClass >
  inline std::map< int,
                   PointerClass >&
  MapAndVectorAndBools< PointerClass >::getMap()
  {
    return codeMap;
  }

  template< class PointerClass >
  inline std::vector< PointerClass >&
  MapAndVectorAndBools< PointerClass >::getVector()
  {
    return pointerGroup;
  }

  template< class PointerClass >
  inline bool
  MapAndVectorAndBools< PointerClass >::getBool()
  {
    return boolReference;
  }

  template< class PointerClass >
  inline std::vector< bool > const*
  MapAndVectorAndBools< PointerClass >::getFlags()
  {
    return boolVector;
  }

}

#endif /* MAPANDVECTORANDBOOLS_HPP_ */
