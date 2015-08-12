/*
 * SinglyIndexedMultipleEntriesBlock.hpp
 *
 *  Created on: Apr 1, 2012 (really!)
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef SINGLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_
#define SINGLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_

#include "../../MEC/ExtendedMass.hpp"
#include "IndexedBlockTemplate.hpp"
#include "InterpreterClasses/MultipleSinglyIndexed.hpp"

namespace LHPC
{
  namespace SLHA
  {
    /* this template class interprets all the blocks with the same name, though
     * differing scale values, which are interpreted as having a single int
     * index which does not have to have entries for each value (nor even
     * necessarily positive index values), allowing for multiple entries with
     * the same index.
     */
    template< class ValueClass >
    class SinglyIndexedMultipleEntriesBlock : public IndexedBlockTemplate<
                                                                    ValueClass,
                        InterpreterClass::MultipleSinglyIndexed< ValueClass > >
    {
    public:
      SinglyIndexedMultipleEntriesBlock( std::string const& blockName,
                                         ValueClass const& defaultUnsetValue,
                                         bool const isVerbose = false,
                                         int const indexDigits = 5 );
      virtual
      ~SinglyIndexedMultipleEntriesBlock();

      std::list< ValueClass* >
      operator()( int const soughtIndex );
      // this returns operator() of the lowest-scale interpreter.
      std::list< ValueClass const* >
      operator()( int const soughtIndex ) const;
      // const version of above.
      bool
      hasEntry( int const soughtIndex ) const;
      // this returns hasEntry( soughtIndex ) of the lowest-scale interpreter.
      virtual bool
      isFmassBlock() const;
      // this returns false. only a specific derived class should over-ride it
      // to return true if it is actually an interpreter for a FMASS block.
      virtual std::multimap< int, RunningConstant > const*
      getFmassMap() const;
      // this returns NULL. only a specific derived class should over-ride it
      // to return a non-NULL pointer.
      virtual bool
      isFmasserrBlock() const;
      // this returns false. only a specific derived class should over-ride it
      // to return true if it is actually an interpreter for a FMASS block.
      virtual std::multimap< int, RunningConstantError > const*
      getFmasserrMap() const;
      // this returns NULL. only a specific derived class should over-ride it
      // to return a non-NULL pointer.


    protected:
      bool isFmassBlockFlag;
      bool isFmasserrBlockFlag;
    };





    template< class ValueClass >
    inline
    SinglyIndexedMultipleEntriesBlock< ValueClass
                                          >::SinglyIndexedMultipleEntriesBlock(
                                                  std::string const& blockName,
                                           ValueClass const& defaultUnsetValue,
                                                          bool const isVerbose,
                                                      int const indexDigits ) :
        IndexedBlockTemplate< ValueClass,
                       InterpreterClass::MultipleSinglyIndexed< ValueClass > >(
                                                                     blockName,
                                                             defaultUnsetValue,
                                                                     isVerbose,
                                                         std::vector< int >( 1,
                                                               indexDigits ) ),
        isFmassBlockFlag( false ),
        isFmasserrBlockFlag( false )
    {
      if( this->nameMatches( "FMASS" ) )
      {
        isFmassBlockFlag = true;
      }
      else if( this->nameMatches( "FMASSERR" ) )
      {
        isFmasserrBlockFlag = true;
      }
    }

    template< class ValueClass >
    inline
    SinglyIndexedMultipleEntriesBlock< ValueClass
                                        >::~SinglyIndexedMultipleEntriesBlock()
    {
      // does nothing.
    }


    template< class ValueClass >
    inline std::list< ValueClass* >
    SinglyIndexedMultipleEntriesBlock< ValueClass >::operator()(
                                                        int const soughtIndex )
    // this returns operator() of the lowest-scale interpreter.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( soughtIndex );
    }

    template< class ValueClass >
    inline std::list< ValueClass const* >
    SinglyIndexedMultipleEntriesBlock< ValueClass >::operator()(
                                                  int const soughtIndex ) const
    // const version of above.
    {
      return this->dataBlocks[ this->lowestScaleIndex ]( soughtIndex );
    }

    template< class ValueClass >
    inline bool
    SinglyIndexedMultipleEntriesBlock< ValueClass >::hasEntry(
                                                  int const soughtIndex ) const
    // this returns hasEntry( soughtIndex ) of the lowest-scale interpreter.
    {
      return
      this->dataBlocks[ this->lowestScaleIndex ].hasEntry( soughtIndex );
    }

    template< class ValueClass >
    inline bool
    SinglyIndexedMultipleEntriesBlock< ValueClass >::isFmassBlock() const
    // this returns false. only a specific derived class should over-ride it
    // to return true if it is actually an interpreter for an FMASS block.
    {
      return isFmassBlockFlag;
    }

    template< class ValueClass >
    inline std::multimap< int, RunningConstant > const*
    SinglyIndexedMultipleEntriesBlock< ValueClass >::getFmassMap() const
    // this returns NULL. only a specific derived class should over-ride it
    // to return a non-NULL pointer.
    {
      return NULL;
    }

    template<>
    inline std::multimap< int, RunningConstant > const*
    SinglyIndexedMultipleEntriesBlock< RunningConstant >::getFmassMap() const
    // this over-rides the default to return a non-NULL pointer if appropriate.
    {
      if( isFmassBlockFlag )
      {
        return &(this->dataBlocks[ this->lowestScaleIndex ].getValueMap());
      }
      else
      {
        return NULL;
      }
    }

    template< class ValueClass >
    inline bool
    SinglyIndexedMultipleEntriesBlock< ValueClass >::isFmasserrBlock() const
    // this returns false. only a specific derived class should over-ride it
    // to return true if it is actually an interpreter for an FMASS block.
    {
      return isFmasserrBlockFlag;
    }

    template< class ValueClass >
    inline std::multimap< int, RunningConstantError > const*
    SinglyIndexedMultipleEntriesBlock< ValueClass >::getFmasserrMap() const
    // this returns NULL. only a specific derived class should over-ride it
    // to return a non-NULL pointer.
    {
      return NULL;
    }

    template<>
    inline std::multimap< int, RunningConstantError > const*
    SinglyIndexedMultipleEntriesBlock< RunningConstantError >::getFmasserrMap(
                                                                        ) const
    // this over-rides the default to return a non-NULL pointer if appropriate.
    {
      if( isFmasserrBlockFlag )
      {
        return &(this->dataBlocks[ this->lowestScaleIndex ].getValueMap());
      }
      else
      {
        return NULL;
      }
    }

  }  // end of SLHA namespace

}  // end of LHPC namespace

#endif /* SINGLYINDEXEDMULTIPLEENTRIESBLOCK_HPP_ */
