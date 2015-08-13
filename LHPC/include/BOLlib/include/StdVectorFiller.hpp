/*
 * StdVectorFiller.hpp
 *
 *  Created on: Jan 20, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *
 *      This file is part of BOLlib, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.BOLlib.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef STDVECTORFILLER_HPP_
#define STDVECTORFILLER_HPP_

#include <vector>
#include <string>

namespace BOL
{
  /* this is a template class just to make it easier to initialize a
   * std::vector with a heterogeneous initial set of elements. the constructor
   * takes the first element of the vector as its argument, & subsequent
   * elements are added with operator(), which returns a reference to this
   * instance, so that the operator() calls can be chained, & then a reference
   * to its internal temporary std::vector< StoredClass > is returned with
   * v( StoredClass const& ), which takes the final element as its argument.
   * hence for example one could have:
   * std::vector< int > const
   * constantOneToFive( BOL::StdVectorFiller( 1 )( 2 )( 3 )( 4 ).end( 5 ) );
   */
  template< typename StoredClass >
  class StdVectorFiller
  {
  public:
    StdVectorFiller( StoredClass const& firstElement );
    ~StdVectorFiller();

    StdVectorFiller&
    operator()( StoredClass const& nextElement );
    std::vector< StoredClass > const&
    end( StoredClass const& lastElement );
    std::vector< StoredClass > const&
    e( StoredClass const& lastElement ) { return end( lastElement ); }


  protected:
    std::vector< StoredClass > stdVector;
  };
  typedef StdVectorFiller< int > Vi;
  typedef StdVectorFiller< double > Vd;
  typedef StdVectorFiller< bool > Vb;
  typedef StdVectorFiller< std::string > Vs;
  // these probably shouldn't be used, but they save space...
  // e.g. std::vector< int > const oneToFour( BOL::Vi( 1 )( 2 )( 3 ).e( 4 ) );



  template< typename StoredClass >
  inline
  StdVectorFiller< StoredClass >::StdVectorFiller(
                                            StoredClass const& firstElement ) :
      stdVector( 1,
                 firstElement )
  {
    // just an initialization list.
  }

  template< typename StoredClass >
  inline
  StdVectorFiller< StoredClass >::~StdVectorFiller()
  {
    // does nothing.
  }


  template< typename StoredClass >
  inline StdVectorFiller< StoredClass >&
  StdVectorFiller< StoredClass >::operator()( StoredClass const& nextElement )
  {
    stdVector.push_back( nextElement );
    return *this;
  }

  template< typename StoredClass >
  inline std::vector< StoredClass > const&
  StdVectorFiller< StoredClass >::end( StoredClass const& lastElement )
  {
    stdVector.push_back( lastElement );
    return stdVector;
  }

}

#endif /* STDVECTORFILLER_HPP_ */
