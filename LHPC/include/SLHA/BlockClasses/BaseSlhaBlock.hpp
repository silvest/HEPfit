/*
 * BaseSlhaBlock.hpp
 *
 *  Created on: Mar 15, 2012
 *      Author: Ben O'Leary (benjamin.oleary@gmail.com)
 *      Copyright 2012 Ben O'Leary
 *
 *      This file is part of LesHouchesParserClasses, released under the
 *      GNU General Public License. Please see the accompanying
 *      README.LHPC_CPP.txt file for a full list of files, brief documentation
 *      on how to use these classes, and further details on the license.
 */

#ifndef BASESLHABLOCK_HPP_
#define BASESLHABLOCK_HPP_

#include <string>
#include <map>
#include "BOLlib/include/BOLlib.hpp"
#include "BaseStringBlock.hpp"
#include "../../MEC/RunningConstant.hpp"
#include "../../MEC/RunningConstantError.hpp"

namespace LHPC
{
  namespace SLHA
  {
    typedef
    BOL::PushedToObserver< BlockClass::BaseStringBlock > StringBlockObserver;

    /* this abstract base class allows other classes to know basic functions
     * of the SLHA block interpreter class without needing to know the specific
     * interpretation.
     */
    class BaseSlhaBlock : public StringBlockObserver
    {
    public:
      BaseSlhaBlock( std::string const& blockName );
      virtual
      ~BaseSlhaBlock();

      std::string const&
      getName() const;
      // this returns the name in uppercase.
      bool
      nameMatches( std::string const& nameToCompare ) const;
      // this returns true if nameToCompare matches the block name ignoring
      // case.
      virtual bool
      isFmassBlock() const;
      // this returns false. only a specific derived class should over-ride it
      // to return true if it is actually an interpreter for an FMASS block.
      virtual std::multimap< int, RunningConstant > const*
      getFmassMap() const;
      // this returns NULL. only a specific derived class should over-ride it
      // to return a non-NULL pointer.
      virtual bool
      isFmasserrBlock() const;
      // this returns false. only a specific derived class should over-ride it
      // to return true if it is actually an interpreter for an FMASSERR block.
      virtual std::multimap< int, RunningConstantError > const*
      getFmasserrMap() const;
      // this returns NULL. only a specific derived class should over-ride it
      // to return a non-NULL pointer.
      virtual bool
      isMassBlock() const;
      // this returns false. only a specific derived class should over-ride it
      // to return true if it is actually an interpreter for a MASS block.
      virtual std::map< int, double > const*
      getMassMap() const;
      // this returns NULL. only a specific derived class should over-ride it
      // to return a non-NULL pointer.


    protected:
      std::string blockName;
    };




    inline std::string const&
    BaseSlhaBlock::getName() const
    // this returns the name in uppercase.
    {
      return blockName;
    }

    inline bool
    BaseSlhaBlock::nameMatches( std::string const& nameToCompare ) const
    // this returns true if nameToCompare matches the block name ignoring
    // case.
    {
      return BOL::StringParser::stringsMatchIgnoringCase( nameToCompare,
                                                          blockName );
    }

    inline bool
    BaseSlhaBlock::isFmassBlock() const
    // this returns false. only a specific derived class should over-ride it
    // to return true if it is actually an interpreter for an FMASS block.
    {
      return false;
    }

    inline std::multimap< int, RunningConstant > const*
    BaseSlhaBlock::getFmassMap() const
    {
      return NULL;
    }

    inline bool
    BaseSlhaBlock::isFmasserrBlock() const
    // this returns false. only a specific derived class should over-ride it
    // to return true if it is actually an interpreter for an FMASS block.
    {
      return false;
    }

    inline std::multimap< int, RunningConstantError > const*
    BaseSlhaBlock::getFmasserrMap() const
    {
      return NULL;
    }

    inline bool
    BaseSlhaBlock::isMassBlock() const
    // this returns false. only a specific derived class should over-ride it
    // to return true if it is actually an interpreter for a MASS block.
    {
      return false;
    }

    inline std::map< int, double > const*
    BaseSlhaBlock::getMassMap() const
    {
      return NULL;
    }

  }

}

#endif /* BASESLHABLOCK_HPP_ */
