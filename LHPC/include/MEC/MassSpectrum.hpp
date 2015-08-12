/*
 * MassSpectrum.hpp
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

#ifndef MASSSPECTRUM_HPP_
#define MASSSPECTRUM_HPP_

#include <iostream>
#include <stdexcept>
#include <map>
#include "MassEigenstate.hpp"
#include "SpectrumUpdater.hpp"
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  // this class holds a set of MassEigenstate instances & can update them with
  // a pair of maps of codes to masses.
  class MassSpectrum : public BOL::PushedToObserver< SpectrumUpdater >
  {
  public:
    enum defaultFlags
    {
      escapesDetector = 0,
      isJet = 1,
      isLightLepton = 2,
      sizeOfEnum = 3
    };
    MassSpectrum( bool const isVerbose = false,
                  std::vector< bool > const* defaultFlags = NULL );
    virtual
    ~MassSpectrum();

    MassEigenstate&
    operator[]( int const pdgCode );
    MassEigenstate const&
    operator[]( int const pdgCode ) const;
    MassEigenstate*
    getMassEigenstate( int const pdgCode );
    MassEigenstate const*
    getMassEigenstate( int const pdgCode ) const;
    std::vector< MassEigenstate* >&
    getMassEigenstateSet();
    std::vector< MassEigenstate* > const&
    getMassEigenstateSet() const;
    MassEigenstate&
    ensureMassEigenstateExists( int const pdgCode );
    MassSpectrum&
    clearMassesAndDecays();
    bool
    getVerbosity() const;
    virtual void
    respondToObservedSignal();
    // the default is over-ridden to call clearMassesAndDecays().
    virtual void
    respondToPush( SpectrumUpdater const& pushedValue );


  protected:
    static std::vector< bool > const defaultBoolVector;
    static std::vector< bool > const defaultEscapesDetectorBoolVector;
    static std::vector< bool > const defaultIsJetBoolVector;
    static std::vector< bool > const defaultIsLightLeptonBoolVector;

    static MassEigenstate&
    findMassEigenstateReference( int const pdgCode,
                               MassEigenstateCodeToPointerMap const& codeMap );
    static void
    warnThatMassEigenstateWasNotFound( int const pdgCode );

    std::vector< MassEigenstate* > allMassEigenstates;
    std::vector< MassEigenstate* > unknownMassEigenstates;
    MassEigenstateCodeToPointerMap pdgCodeMap;
    bool const isVerbose;
    MassEigenstateMapAndVectorAndBools mapAndVectorAndBools;
  };



  inline MassEigenstate&
  MassSpectrum::operator[]( int const pdgCode )
  {
    return findMassEigenstateReference( pdgCode,
                                        pdgCodeMap );
  }

  inline MassEigenstate const&
  MassSpectrum::operator[]( int const pdgCode ) const
  {
    return findMassEigenstateReference( pdgCode,
                                        pdgCodeMap );
  }

  inline MassEigenstate*
  MassSpectrum::getMassEigenstate( int const pdgCode )
  {
    return MassEigenstate::findPointerWithCode( pdgCode,
                                                pdgCodeMap );
  }

  inline MassEigenstate const*
  MassSpectrum::getMassEigenstate( int const pdgCode ) const
  {
    return MassEigenstate::findPointerWithCode( pdgCode,
                                                pdgCodeMap );
  }

  inline std::vector< MassEigenstate* >&
  MassSpectrum::getMassEigenstateSet()
  {
    return allMassEigenstates;
  }

  inline std::vector< MassEigenstate* > const&
  MassSpectrum::getMassEigenstateSet() const
  {
    return allMassEigenstates;
  }

  inline MassSpectrum&
  MassSpectrum::clearMassesAndDecays()
  {
    for( int clearingIndex( allMassEigenstates.size() - 1 );
         0 <= clearingIndex;
         --clearingIndex )
    {
      allMassEigenstates[ clearingIndex ]->clearMassesAndDecays();
    }
    return *this;
  }

  inline bool
  MassSpectrum::getVerbosity() const
  {
    return isVerbose;
  }

  inline void
  MassSpectrum::respondToObservedSignal()
  // the default is over-ridden to call clearMassesAndDecays().
  {
    clearMassesAndDecays();
  }

  inline void
  MassSpectrum::respondToPush( SpectrumUpdater const& pushedValue )
  {
    pushedValue.updateMassEigenstates( pdgCodeMap );
  }

  inline MassEigenstate&
  MassSpectrum::findMassEigenstateReference( int const pdgCode,
                                MassEigenstateCodeToPointerMap const& codeMap )
  {
    MassEigenstate*
    massEigenstatePointer( MassEigenstate::findPointerWithCode( pdgCode,
                                                                codeMap ) );
    if( NULL != massEigenstatePointer )
    {
      return *massEigenstatePointer;
    }
    else
    {
      std::string
      errorMessage( "MassSpectrum::findMassEigenstateReference( " );
      errorMessage.append( BOL::StringParser::intToString( pdgCode,
                                                           1 ) );
      errorMessage.append( ") out of range." );
      warnThatMassEigenstateWasNotFound( pdgCode );
      throw std::out_of_range( errorMessage );
    }
  }

  inline void
  MassSpectrum::warnThatMassEigenstateWasNotFound( int const pdgCode )
  {
    std::cout
    << std::endl
    << "LhaParsing::error! MassEigenstate::findPointerWithCode could not"
    << " find the MassEigenstate with particle code " << pdgCode
    << ", so it is returning a NULL pointer. MassSpectrum is throwing an"
    << " out-of-range exception because of this.";
    std::cout << std::endl;
  }

}

#endif /* MASSSPECTRUM_HPP_ */
