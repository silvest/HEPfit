/*
 * SpectrumUpdater.hpp
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

#ifndef SPECTRUMUPDATER_HPP_
#define SPECTRUMUPDATER_HPP_

#include <map>
#include <vector>
#include "BOLlib/include/BOLlib.hpp"
#include "MassEigenstate.hpp"
#include "RunningConstant.hpp"
#include "RunningConstantError.hpp"

namespace LHPC
{
  // this is a class that holds information to allow a MassSpectrum instance to
  // update its masses and decays.
  class SpectrumUpdater
  {
  public:
    SpectrumUpdater();
    ~SpectrumUpdater();

    void
    setFmassMap( std::multimap< int, RunningConstant > const* fmassMap );
    void
    setFmasserrMap(
               std::multimap< int, RunningConstantError > const* fmasserrMap );
    void
    setMassMap( std::map< int, double > const* massMap );
    bool
    isHoldingDecay() const;
    void
    clearDecay();
    void
    recordDecayHeader( std::string const& decayerCode,
                       std::string const& decayWidth );
    void
    recordDecayLine( BOL::VectorlikeArray< std::string > const& decayLine );
    void
    updateMassEigenstates( MassEigenstateCodeToPointerMap& codeMap ) const;


  protected:
    std::multimap< int, RunningConstant > const* fmassMap;
    std::multimap< int, RunningConstantError > const* fmasserrMap;
    std::map< int, double > const* massMap;
    bool isHoldingDecayFlag;
    int decayerCode;
    double decayWidth;
    std::vector< std::pair< double, std::vector< int > > >
    decayChannels;
    std::pair< double, std::vector< int > > branchingRatioAndProducts;
  };





  inline void
  SpectrumUpdater::setFmassMap(
                        std::multimap< int, RunningConstant > const* fmassMap )
  {
    this->fmassMap = fmassMap;
  }

  inline void
  SpectrumUpdater::setFmasserrMap(
                std::multimap< int, RunningConstantError > const* fmasserrMap )
  {
    this->fmasserrMap = fmasserrMap;
  }

  inline void
  SpectrumUpdater::setMassMap( std::map< int, double > const* massMap )
  {
    this->massMap = massMap;
  }

  inline bool
  SpectrumUpdater::isHoldingDecay() const
  {
    return isHoldingDecayFlag;
  }

  inline void
  SpectrumUpdater::clearDecay()
  {
    decayChannels.clear();
    isHoldingDecayFlag = false;
  }

  inline void
  SpectrumUpdater::recordDecayHeader( std::string const& decayerCode,
                                      std::string const& decayWidth )
  {
    this->decayerCode = BOL::StringParser::stringToInt( decayerCode );
    this->decayWidth = BOL::StringParser::stringToDouble( decayWidth );
    isHoldingDecayFlag = true;
  }

  inline void
  SpectrumUpdater::recordDecayLine(
                         BOL::VectorlikeArray< std::string > const& decayLine )
  {
    branchingRatioAndProducts.first
    = BOL::StringParser::stringToDouble( decayLine[ 0 ] );
    branchingRatioAndProducts.second.clear();
    for( int whichWord( decayLine.getLastIndex() );
         1 < whichWord;
         --whichWord )
    {
      branchingRatioAndProducts.second.push_back(
                    BOL::StringParser::stringToInt( decayLine[ whichWord ] ) );
    }
    decayChannels.push_back( branchingRatioAndProducts );
  }

}

#endif /* SPECTRUMUPDATER_HPP_ */
