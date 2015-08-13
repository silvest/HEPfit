/*
 * MassEigenstate.hpp
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

#ifndef MASSEIGENSTATE_HPP_
#define MASSEIGENSTATE_HPP_

#include <map>
#include <iostream>
#include "BOLlib/include/BOLlib.hpp"
#include "ExtendedMass.hpp"
#include "PointersWithValue.hpp"
#include "MapAndVectorAndBools.hpp"

namespace LHPC
{
  /* this is a class to hold the relevant properties of a particle mass
   * eigenstate as far as the SLHA format is concerned: masses, decay width,
   * branching ratios.
   * charge conjugate pairs are kept separate, but linked.
   */
  class MassEigenstate
  {
  public:
    typedef
    MapAndVectorAndBools< MassEigenstate* > MassEigenstateMapVectorBools;
    typedef std::map< int, MassEigenstate* > MassEigenstateCodeMap;
    typedef
    PointersWithValue< MassEigenstate const*, double >
    MassEigenstatesPairedWithBr;
    typedef std::vector< MassEigenstate* > MassEigenstateVector;

    static MassEigenstate*
    findPointerWithCode( int pdgCode,
                         MassEigenstateCodeMap const& codeMap );
    static bool
    pointerHasCode( MassEigenstate const* const massEigenstatePointer,
                    int const pdgCode );

    MassEigenstate( int const pdgCode,
                    MassEigenstateMapVectorBools& mapAndVectorAndBools,
                    bool const isSelfConjugate,
                    std::string const& asciiName,
                    std::string const& latexName,
                  double const defaultResetMass = BOL::UsefulStuff::notANumber,
               double const defaultDecayWidth = BOL::UsefulStuff::notANumber );
    MassEigenstate( int const firstPdgCode,
                    int const secondPdgCode,
                    MassEigenstateMapVectorBools& mapAndVectorAndBools,
                    bool const isSelfConjugate,
                    std::string const& asciiName,
                    std::string const& latexName,
                  double const defaultResetMass = BOL::UsefulStuff::notANumber,
               double const defaultDecayWidth = BOL::UsefulStuff::notANumber );
    MassEigenstate( std::vector< int > const& pdgCodes,
                    MassEigenstateMapVectorBools& mapAndVectorAndBools,
                    bool const isSelfConjugate,
                    std::string const& asciiName,
                    std::string const& latexName,
                  double const defaultResetMass = BOL::UsefulStuff::notANumber,
               double const defaultDecayWidth = BOL::UsefulStuff::notANumber );
    MassEigenstate( MassEigenstate const& copySource,
                    MassEigenstateMapVectorBools& mapAndVectorAndBools );
    // this last version copies as the charge conjugate of copySource:
    MassEigenstate( MassEigenstate& copySource,
                    std::string const& asciiName,
                    std::string const& latexName );
    ~MassEigenstate();

    bool
    isSelfConjugate() const;
    MassEigenstate&
    getChargeConjugate();
    MassEigenstate const&
    getChargeConjugate() const;
    void
    setChargeConjugate( MassEigenstate* const chargeConjugate );
    void
    setToBeSelfConjugate();
    void
    setToBeChargeConjugate( MassEigenstate* const chargeConjugate );
    bool
    hasCode( int const pdgCode ) const;
    int
    getCode() const;
    std::vector< int > const&
    getAllCodes() const;
    MassEigenstate&
    addCode( int extraCode );
    bool
    hasMassBeenRecorded() const;
    double
    getSignedMass() const;
    // this defaults to the pole mass, if available.
    double
    getAbsoluteMass() const;
    // this defaults to the pole mass, if available.
    ExtendedMass const&
    getExtendedMass() const;
    // this only returns the 1st of the extended masses.
    std::vector< ExtendedMass* > const&
    getAllRecordedMasses() const;
    void
    recordMass( double const massValue,
                double const minusUncertainty = 0.0,
                double const plusUncertainty = 0.0,
                int const schemeType = 0,
                double const evaluationScale = 0.0 );
    void
    recordMassError( double const minusUncertainty,
                     double const plusUncertainty,
                     int const schemeType,
                     double const evaluationScale );
    bool
    haveDecaysBeenRecorded() const;
    double
    getDecayWidth() const;
    std::vector< MassEigenstatesPairedWithBr* >&
    getDecaySet();
    std::vector< MassEigenstatesPairedWithBr* > const&
    getDecaySet() const;
    MassEigenstatesPairedWithBr&
    getDecay( int const whichDecay );
    MassEigenstatesPairedWithBr const&
    getDecay( int const whichDecay ) const;
    void
    setDecayWidth( double const decayWidth );
    MassEigenstate*
    clearMassesAndDecays();
    void
    recordDecay( MassEigenstatesPairedWithBr const& decayToRecord );
    void
    recordChargeConjugateOfDecay(
                            MassEigenstatesPairedWithBr const& decayToRecord );
    void
    recordDecayAsDefault( MassEigenstatesPairedWithBr const& decayToRecord );
    void
    recordChargeConjugateOfDecayAsDefault(
                            MassEigenstatesPairedWithBr const& decayToRecord );
    double
    getExactMatchBranchingRatio( MassEigenstate const& firstProduct,
                                 MassEigenstate const& secondProduct ) const;
    double
    getExactMatchBranchingRatio( MassEigenstate const& firstProduct,
                                 MassEigenstate const& secondProduct,
                                 MassEigenstate const& thirdProduct ) const;
    double
    getExactMatchBranchingRatioWithSortedList(
           std::list< MassEigenstate const* > const& sortedProductList ) const;
    double
    getExactMatchBranchingRatio(
                       std::list< MassEigenstate const* >& productList ) const;
    double
    getExactMatchBranchingRatio(
               std::vector< MassEigenstate const* > const& productList ) const;
    double
    getSubsetMatchBranchingRatioWithSortedLists(
                          std::list< MassEigenstate const* > const& soughtList,
                  std::list< MassEigenstate const* > const& vetoedList ) const;
    double
    getSubsetMatchBranchingRatio(
                                std::list< MassEigenstate const* >& soughtList,
                        std::list< MassEigenstate const* >& vetoedList ) const;
    double
    getSubsetMatchBranchingRatio(
                      std::vector< MassEigenstate const* > const& soughtVector,
              std::vector< MassEigenstate const* > const& vetoedVector ) const;
    double
    getExactMatchBranchingRatio( int const firstProductCode,
                                 int const secondProductCode ) const;
    double
    getExactMatchBranchingRatio( int const firstProductCode,
                                 int const secondProductCode,
                                 int const thirdProductCode ) const;
    double
    getExactMatchBranchingRatio(
                               std::list< int > const& productCodeList ) const;
    double
    getExactMatchBranchingRatio(
                             std::vector< int > const& productCodeList ) const;
    double
    getSubsetMatchBranchingRatio( std::list< int > const& soughtCodeList,
                                std::list< int > const& vetoedCodeList ) const;
    double
    getSubsetMatchBranchingRatio( std::vector< int > const& soughtCodeVector,
                            std::vector< int > const& vetoedCodeVector ) const;
    std::string const&
    getAsciiName() const;
    MassEigenstate&
    setAsciiName( std::string const& asciiName );
    std::string const&
    getLatexName() const;
    MassEigenstate&
    setLatexName( std::string const& latexName );
    std::vector< bool > const&
    getFlags() const;
    MassEigenstate&
    setFlags( std::vector< bool > const* const flagBools );
    MassEigenstate*
    findPointerWithCode( int const pdgCode );
    MassEigenstate const*
    findPointerWithCode( int const pdgCode ) const;
    bool
    getVerbosity() const;


  protected:
    MassEigenstate* chargeConjugate;
    bool isSelfConjugateFlag;
    std::vector< int > identifyingPdgCodes;
    MassEigenstateCodeMap& pdgCodeMap;
    std::pair< int, MassEigenstate* > mapFiller;
    bool massRecorded;
    double const defaultResetMass;
    double signedDefaultMass;
    // this defaults to the pole mass, if available.
    double absoluteDefaultMass;
    // this defaults to the pole mass, if available.
    BOL::VectorlikeArray< ExtendedMass > runningMasses;
    std::vector< ExtendedMass* > runningMassesAsVector;
    bool decaysRecorded;
    double const defaultDecayWidth;
    double decayWidth;
    BOL::VectorlikeArray< MassEigenstatesPairedWithBr > defaultDecaySet;
    BOL::VectorlikeArray< MassEigenstatesPairedWithBr > decaySet;
    std::vector< MassEigenstatesPairedWithBr* > decaySetAsVector;
    std::string asciiName;
    std::string latexName;
    bool const isVerbose;
    std::vector< bool > const* flagBools;
    MassEigenstateVector& setOfPointersOfMassEigenstateGroup;

    void
    constructorBodyFunction();
    void
    addToCodeMap( int const positiveExtraCode );
    void
    prepareToRecordDecay();
  };
  typedef
  MapAndVectorAndBools< MassEigenstate* > MassEigenstateMapAndVectorAndBools;
  typedef std::map< int, MassEigenstate* > MassEigenstateCodeToPointerMap;



  inline bool
  MassEigenstate::pointerHasCode(
                             MassEigenstate const* const massEigenstatePointer,
                                  int const pdgCode )
  {
    return massEigenstatePointer->hasCode( pdgCode );
  }


  inline bool
  MassEigenstate::isSelfConjugate() const
  {
    return isSelfConjugateFlag;
  }

  inline MassEigenstate&
  MassEigenstate::getChargeConjugate()
  {
    return *chargeConjugate;
  }

  inline MassEigenstate const&
  MassEigenstate::getChargeConjugate() const
  {
    return *chargeConjugate;
  }

  inline void
  MassEigenstate::setChargeConjugate( MassEigenstate* const chargeConjugate )
  {
    this->chargeConjugate = chargeConjugate;
  }

  inline void
  MassEigenstate::setToBeSelfConjugate()
  {
    chargeConjugate = this;
    /* not only does this need to set its chargeConjugate pointer back to
     * itself, but the codes need to be set to positive & pdgCodeMap needs to
     * only look for this instance rather than its old chargeConjugate.
     */
    for( int codeIndex( identifyingPdgCodes.size() - 1 );
         0 <= codeIndex;
         --codeIndex )
    {
      MassEigenstateCodeMap::iterator codeFinder;
      if( 0 > identifyingPdgCodes[ codeIndex ] )
      {
        identifyingPdgCodes[ codeIndex ] = -identifyingPdgCodes[ codeIndex ];
        codeFinder = pdgCodeMap.find( identifyingPdgCodes[ codeIndex ] );
        if( pdgCodeMap.end() != codeFinder )
        {
          codeFinder->second = this;
        }
      }
    }
  }

  inline bool
  MassEigenstate::hasCode( int const pdgCode ) const
  {
    bool returnBool( false );
    for( int codeSeeker( identifyingPdgCodes.size() - 1 );
         0 <= codeSeeker;
         --codeSeeker )
    {
      if( pdgCode == identifyingPdgCodes[ codeSeeker ] )
      {
        returnBool = true;
        codeSeeker = -1;
      }
    }
    return returnBool;
  }

  inline int
  MassEigenstate::getCode() const
  {
    return identifyingPdgCodes.front();
  }

  inline std::vector< int > const&
  MassEigenstate::getAllCodes() const
  {
    return identifyingPdgCodes;
  }

  inline bool
  MassEigenstate::hasMassBeenRecorded() const
  {
    return massRecorded;
  }

  inline double
  MassEigenstate::getSignedMass() const
  // this defaults to the pole mass, if available.
  {
    return signedDefaultMass;
  }

  inline double
  MassEigenstate::getAbsoluteMass() const
   // this defaults to the pole mass, if available.
  {
    return absoluteDefaultMass;
  }

  inline ExtendedMass const&
  MassEigenstate::getExtendedMass() const
  // this only returns the 1st of the extended masses.
  {
    return runningMasses.getFront();
  }

  inline std::vector< ExtendedMass* > const&
  MassEigenstate::getAllRecordedMasses() const
  {
    return runningMassesAsVector;
  }

  inline void
  MassEigenstate::recordMassError( double const minusUncertainty,
                                   double const plusUncertainty,
                                   int const schemeType,
                                   double const evaluationScale )
  {
    for( int whichMass( runningMasses.getLastIndex() );
         0 <= whichMass;
         --whichMass )
    {
      if( ( schemeType == runningMasses[ whichMass ].getScheme() )
          &&
          ( evaluationScale == runningMasses[ whichMass ].getScale() ) )
      {
        runningMasses[ whichMass ].setUncertainties( minusUncertainty,
                                                     plusUncertainty );
        whichMass = -1;
      }
    }
  }

  inline bool
  MassEigenstate::haveDecaysBeenRecorded() const
  {
    return decaysRecorded;
  }

  inline double
  MassEigenstate::getDecayWidth() const
  {
    return decayWidth;
  }

  inline void
  MassEigenstate::setDecayWidth( double const decayWidth )
  {
    this->decayWidth = decayWidth;
    decaysRecorded = true;
    // I don't think that I want to set up an elaborate check that a width
    // & 100% sum of BRs have been recorded to set decaysRecorded.
  }

  inline std::vector< MassEigenstate::MassEigenstatesPairedWithBr* >&
  MassEigenstate::getDecaySet()
  {
    return decaySetAsVector;
  }

  inline std::vector< MassEigenstate::MassEigenstatesPairedWithBr* > const&
  MassEigenstate::getDecaySet() const
  {
    return decaySetAsVector;
  }

  inline MassEigenstate::MassEigenstatesPairedWithBr&
  MassEigenstate::getDecay( int const whichDecay )
  {
    return *(decaySetAsVector[ whichDecay ]);
  }

  inline MassEigenstate::MassEigenstatesPairedWithBr const&
  MassEigenstate::getDecay( int const whichDecay ) const
  {
    return *(decaySetAsVector[ whichDecay ]);
  }

  inline MassEigenstate*
  MassEigenstate::clearMassesAndDecays()
  {
    runningMasses.clearEntries();
    runningMassesAsVector.clear();
    signedDefaultMass = defaultResetMass;
    absoluteDefaultMass = defaultResetMass;
    if( 0.0 > absoluteDefaultMass )
    {
      absoluteDefaultMass = -absoluteDefaultMass;
    }
    massRecorded = false;
    decaySet.becomeCopyOf( defaultDecaySet );
    decaySetAsVector.clear();
    for( int whichDecay( 0 );
         decaySet.getSize() > whichDecay;
         ++whichDecay )
    {
      decaySetAsVector.push_back( decaySet.getPointer( whichDecay ) );
    }
    decayWidth = defaultDecayWidth;
    decaysRecorded = false;
    return this;
  }

  inline void
  MassEigenstate::recordDecay(
                             MassEigenstatesPairedWithBr const& decayToRecord )
  {
    prepareToRecordDecay();
    decaySetAsVector.push_back( decaySet.newEnd().copyAndReturnPointer(
                                                             decayToRecord ) );
    // I decided to keep a separate std::vector in parallel rather than
    // constantly converting decaySet with getAsVector(...).
  }

  inline void
  MassEigenstate::recordChargeConjugateOfDecay(
                             MassEigenstatesPairedWithBr const& decayToRecord )
  {
    prepareToRecordDecay();
    MassEigenstatesPairedWithBr& decayRecorder( decaySet.newEnd() );
    decaySetAsVector.push_back( &decayRecorder );
    // I decided to keep a separate std::vector in parallel rather than
    // constantly converting decaySet with getAsVector(...).
    decayRecorder.clearPointers();
    for( std::list< MassEigenstate const* >::const_iterator
         productIterator( decayToRecord.getPointerList().begin() );
         decayToRecord.getPointerList().end() != productIterator;
         ++productIterator )
    {
      decayRecorder.addPointer( &((*productIterator)->getChargeConjugate()) );
    }
    decayRecorder.setPairedValueAndSortPointers(
                                              decayToRecord.getPairedValue() );
  }

  inline void
  MassEigenstate::recordDecayAsDefault(
                             MassEigenstatesPairedWithBr const& decayToRecord )
  {
    defaultDecaySet.newEnd().becomeCopyOf( decayToRecord );
  }

  inline void
  MassEigenstate::recordChargeConjugateOfDecayAsDefault(
                             MassEigenstatesPairedWithBr const& decayToRecord )
  {
    MassEigenstatesPairedWithBr& decayRecorder( defaultDecaySet.newEnd() );
    decayRecorder.clearPointers();
    for( std::list< MassEigenstate const* >::const_iterator
         productIterator( decayToRecord.getPointerList().begin() );
         decayToRecord.getPointerList().end() != productIterator;
         ++productIterator )
    {
      decayRecorder.addPointer( &((*productIterator)->getChargeConjugate()) );
    }
    decayRecorder.setPairedValueAndSortPointers(
                                              decayToRecord.getPairedValue() );
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio(
                                            MassEigenstate const& firstProduct,
                                    MassEigenstate const& secondProduct ) const
  {
    double returnValue( 0.0 );
    for( int whichDecay( decaySet.getLastIndex() );
         0 <= whichDecay;
         --whichDecay )
    {
      if( decaySet[ whichDecay ].matchesExactly( &firstProduct,
                                                 &secondProduct ) )
      {
        returnValue = decaySet[ whichDecay ].getPairedValue();
        whichDecay = -1;
      }
    }
    return returnValue;
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio(
                                            MassEigenstate const& firstProduct,
                                           MassEigenstate const& secondProduct,
                                     MassEigenstate const& thirdProduct ) const
  {
    std::list< MassEigenstate const* > comparisonList( 1,
                                                       &firstProduct );
    comparisonList.push_back( &secondProduct );
    comparisonList.push_back( &thirdProduct );
    return getExactMatchBranchingRatio( comparisonList );
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatioWithSortedList(
            std::list< MassEigenstate const* > const& sortedProductList ) const
  {
    double returnValue( 0.0 );
    for( int whichDecay( decaySet.getLastIndex() );
         0 <= whichDecay;
         --whichDecay )
    {
      if( decaySet[ whichDecay ].matchesExactly( sortedProductList ) )
      {
        returnValue = decaySet[ whichDecay ].getPairedValue();
        whichDecay = -1;
      }
    }
    return returnValue;
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio(
                        std::list< MassEigenstate const* >& productList ) const
  {
    productList.sort();
    return getExactMatchBranchingRatioWithSortedList( productList );
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio(
                std::vector< MassEigenstate const* > const& productList ) const
  {
    std::list< MassEigenstate const* > comparisonList( productList.begin(),
                                                       productList.end() );
    comparisonList.sort();
    return getExactMatchBranchingRatioWithSortedList( comparisonList );
  }

  inline double
  MassEigenstate::getSubsetMatchBranchingRatioWithSortedLists(
                    std::list< MassEigenstate const* > const& sortedSoughtList,
             std::list< MassEigenstate const* > const& sortedVetoedList ) const
  {
    double returnValue( 0.0 );
    for( int whichDecay( decaySet.getLastIndex() );
         0 <= whichDecay;
         --whichDecay )
    {
      if( decaySet[ whichDecay ].containsAllAsSubset( sortedSoughtList )
          &&
          !(decaySet[ whichDecay ].containsAnyAsSubset( sortedVetoedList )) )
      {
        returnValue += decaySet[ whichDecay ].getPairedValue();
      }
    }
    return returnValue;
  }

  inline double
  MassEigenstate::getSubsetMatchBranchingRatio(
                                std::list< MassEigenstate const* >& soughtList,
                         std::list< MassEigenstate const* >& vetoedList ) const
  {
    soughtList.sort();
    vetoedList.sort();
    return getSubsetMatchBranchingRatioWithSortedLists( soughtList,
                                                        vetoedList );
  }

  inline double
  MassEigenstate::getSubsetMatchBranchingRatio(
                      std::vector< MassEigenstate const* > const& soughtVector,
               std::vector< MassEigenstate const* > const& vetoedVector ) const
  {
    std::list< MassEigenstate const* > soughtList( soughtVector.begin(),
                                                   soughtVector.end() );
    std::list< MassEigenstate const* > vetoedList( vetoedVector.begin(),
                                                   vetoedVector.end() );
    return getSubsetMatchBranchingRatio( soughtList,
                                         vetoedList );
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio( int const firstProductCode,
                                            int const secondProductCode ) const
  {
    double returnValue( 0.0 );
    for( int whichDecay( decaySet.getLastIndex() );
         0 <= whichDecay;
         --whichDecay )
    {
      if( decaySet[ whichDecay ].matchesExactly( firstProductCode,
                                                 secondProductCode,
                                                 &pointerHasCode ) )
      {
        returnValue = decaySet[ whichDecay ].getPairedValue();
        whichDecay = -1;
      }
    }
    return returnValue;
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio( int const firstProductCode,
                                               int const secondProductCode,
                                             int const thirdProductCode ) const
  {
    std::list< MassEigenstate const* > comparisonList( 1,
                                     findPointerWithCode( firstProductCode ) );
    comparisonList.push_back( findPointerWithCode( secondProductCode )  );
    comparisonList.push_back( findPointerWithCode( thirdProductCode )  );
    return getExactMatchBranchingRatio( comparisonList );
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio(
                                std::list< int > const& productCodeList ) const
  {
    std::list< MassEigenstate const* > comparisonList;
    for( std::list< int >::const_iterator
         codeIterator( productCodeList.begin() );
         productCodeList.end() != codeIterator;
         ++codeIterator )
    {
      comparisonList.push_back( findPointerWithCode( *codeIterator )  );
    }
    return getExactMatchBranchingRatio( comparisonList );
  }

  inline double
  MassEigenstate::getExactMatchBranchingRatio(
                              std::vector< int > const& productCodeList ) const
  {
    std::list< MassEigenstate const* > comparisonList;
    for( int codeIndex( productCodeList.size() - 1 );
         0 <= codeIndex;
         --codeIndex )
    {
      comparisonList.push_back( findPointerWithCode(
                                              productCodeList[ codeIndex ] ) );
    }
    return getExactMatchBranchingRatio( comparisonList );
  }

  inline double
  MassEigenstate::getSubsetMatchBranchingRatio(
                                        std::list< int > const& soughtCodeList,
                                 std::list< int > const& vetoedCodeList ) const
  {
    std::list< MassEigenstate const* > soughtPointerList;
    for( std::list< int >::const_iterator
         codeIterator( soughtCodeList.begin() );
         soughtCodeList.end() != codeIterator;
         ++codeIterator )
    {
      soughtPointerList.push_back( findPointerWithCode( *codeIterator )  );
    }
    std::list< MassEigenstate const* > vetoedPointerList;
    for( std::list< int >::const_iterator
         codeIterator( vetoedCodeList.begin() );
         vetoedCodeList.end() != codeIterator;
         ++codeIterator )
    {
      vetoedPointerList.push_back( findPointerWithCode( *codeIterator )  );
    }
    return getSubsetMatchBranchingRatio( soughtPointerList,
                                         vetoedPointerList );
  }

  inline double
  MassEigenstate::getSubsetMatchBranchingRatio(
                                      std::vector< int > const& soughtCodeList,
                               std::vector< int > const& vetoedCodeList ) const
  {
    std::list< MassEigenstate const* > soughtPointerList;
    for( int codeIndex( soughtCodeList.size() - 1 );
         0 <= codeIndex;
         --codeIndex )
    {
      soughtPointerList.push_back( findPointerWithCode(
                                              soughtCodeList[ codeIndex ] )  );
    }
    std::list< MassEigenstate const* > vetoedPointerList;
    for( int codeIndex( vetoedCodeList.size() - 1 );
         0 <= codeIndex;
         --codeIndex )
    {
      vetoedPointerList.push_back( findPointerWithCode(
                                              vetoedCodeList[ codeIndex ] )  );
    }
    return getSubsetMatchBranchingRatio( soughtPointerList,
                                         vetoedPointerList );
  }

  inline std::string const&
  MassEigenstate::getAsciiName() const
  {
    return asciiName;
  }

  inline MassEigenstate&
  MassEigenstate::setAsciiName( std::string const& asciiName )
  {
    this->asciiName.assign( asciiName );
    return *this;
  }

  inline std::string const&
  MassEigenstate::getLatexName() const
  {
    return latexName;
  }

  inline MassEigenstate&
  MassEigenstate::setLatexName( std::string const& latexName )
  {
    this->latexName.assign( latexName );
    return *this;
  }

  inline std::vector< bool > const&
  MassEigenstate::getFlags() const
  {
    return *flagBools;
  }

  inline MassEigenstate&
  MassEigenstate::setFlags( std::vector< bool > const* const flagBools )
  {
    this->flagBools = flagBools;
    return *this;
  }

  inline MassEigenstate*
  MassEigenstate::findPointerWithCode( int const pdgCode )
  {
    return findPointerWithCode( pdgCode,
                                pdgCodeMap );
  }

  inline MassEigenstate const*
  MassEigenstate::findPointerWithCode( int const pdgCode ) const
  {
    return findPointerWithCode( pdgCode,
                                pdgCodeMap );
  }

  inline bool
  MassEigenstate::getVerbosity() const
  {
    return isVerbose;
  }

  inline void
  MassEigenstate::addToCodeMap( int const positiveExtraCode )
  {
    mapFiller.first = positiveExtraCode;
    mapFiller.second = this;
    pdgCodeMap.insert( mapFiller );
  }

  inline void
  MassEigenstate::prepareToRecordDecay()
  {
    if( !decaysRecorded )
    {
      // this ensures that any present default decays are wiped if any decays
      // are recorded.
      decaySet.clearEntries();
      decaySetAsVector.clear();
      decaysRecorded = true;
    }
  }

}

#endif /* MASSEIGENSTATE_HPP_ */
