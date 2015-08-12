/*
 * MassEigenstate.cpp
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

#include "MEC.hpp"

namespace LHPC
{
  MassEigenstate*
  MassEigenstate::findPointerWithCode( int pdgCode,
                                       MassEigenstateCodeMap const& codeMap )
  // this finds the MassEigenstate pointer which corresponds to the requested
  // code in the given map, returning NULL if there is none.
  {
    MassEigenstate* returnPointer( NULL );
    bool flipToChargeConjugate( false );
    // codeMap only knows about positive codes, but negatives codes will
    // return charge-conjugates of the MassEigenstates with the positive codes:
    if( 0 > pdgCode )
    {
      flipToChargeConjugate = true;
      pdgCode = -pdgCode;
    }
    MassEigenstateCodeMap::const_iterator
    returnFromMap( codeMap.find( pdgCode ) );
    if( codeMap.end() != returnFromMap )
      // if pdgCode was found in codeMap...
    {
      if( flipToChargeConjugate )
      {
        returnPointer = (*returnFromMap).second->chargeConjugate;
      }
      else
      {
        returnPointer = (*returnFromMap).second;
      }
    }
    return returnPointer;
  }

  MassEigenstate::MassEigenstate( int const pdgCode,
                            MassEigenstateMapVectorBools& mapAndVectorAndBools,
                                  bool const isSelfConjugate,
                                  std::string const& asciiName,
                                  std::string const& latexName,
                                  double const defaultResetMass,
                                  double const defaultDecayWidth ) :
      chargeConjugate( NULL ),
      isSelfConjugateFlag( isSelfConjugate ),
      identifyingPdgCodes( 1,
                           pdgCode ),
      pdgCodeMap( mapAndVectorAndBools.getMap() ),
      mapFiller( 0,
                 NULL ),
      massRecorded( false ),
      defaultResetMass( defaultResetMass ),
      signedDefaultMass( defaultResetMass ),
      absoluteDefaultMass( defaultResetMass ),
      runningMasses(),
      runningMassesAsVector(),
      decaysRecorded( false ),
      defaultDecayWidth( defaultDecayWidth ),
      decayWidth( defaultDecayWidth ),
      defaultDecaySet(),
      decaySet(),
      decaySetAsVector(),
      asciiName( asciiName ),
      latexName( latexName ),
      isVerbose( mapAndVectorAndBools.getBool() ),
      flagBools( mapAndVectorAndBools.getFlags() ),
      setOfPointersOfMassEigenstateGroup( mapAndVectorAndBools.getVector() )
  {
    constructorBodyFunction();
    /* constructorBodyFunction puts pointers to this instance into pdgCodeMap
     * for each of this instance's positive codes. if this instance has
     * negative codes, they will only get into pdgCodeMap through the
     * charge-conjugate's constructor (or by manually adding them, but that's
     * not recommended).
     */
  }

  MassEigenstate::MassEigenstate( int const firstPdgCode,
                                  int const secondPdgCode,
                            MassEigenstateMapVectorBools& mapAndVectorAndBools,
                                  bool const isSelfConjugate,
                                  std::string const& asciiName,
                                  std::string const& latexName,
                                  double const defaultResetMass,
                                  double const defaultDecayWidth ) :
      chargeConjugate( NULL ),
      isSelfConjugateFlag( isSelfConjugate ),
      identifyingPdgCodes( 2,
                           secondPdgCode ),
      pdgCodeMap( mapAndVectorAndBools.getMap() ),
      mapFiller( 0,
                 NULL ),
      massRecorded( false ),
      defaultResetMass( defaultResetMass ),
      signedDefaultMass( defaultResetMass ),
      absoluteDefaultMass( defaultResetMass ),
      runningMasses(),
      runningMassesAsVector(),
      decaysRecorded( false ),
      defaultDecayWidth( defaultDecayWidth ),
      decayWidth( defaultDecayWidth ),
      defaultDecaySet(),
      decaySet(),
      decaySetAsVector(),
      asciiName( asciiName ),
      latexName( latexName ),
      isVerbose( mapAndVectorAndBools.getBool() ),
      flagBools( mapAndVectorAndBools.getFlags() ),
      setOfPointersOfMassEigenstateGroup( mapAndVectorAndBools.getVector() )
  {
    identifyingPdgCodes.front() = firstPdgCode;
    // this, along with the initialization of identifyingPdgCodes, sorts out
    // the 2 particle codes.
    constructorBodyFunction();
    /* constructorBodyFunction puts pointers to this instance into pdgCodeMap
     * for each of this instance's positive codes. if this instance has
     * negative codes, they will only get into pdgCodeMap through the
     * charge-conjugate's constructor (or by manually adding them, but that's
     * not recommended).
     */
  }

  MassEigenstate::MassEigenstate( std::vector< int > const& pdgCodes,
                            MassEigenstateMapVectorBools& mapAndVectorAndBools,
                                  bool const isSelfConjugate,
                                  std::string const& asciiName,
                                  std::string const& latexName,
                                  double const defaultResetMass,
                                  double const defaultDecayWidth ) :
      chargeConjugate( NULL ),
      isSelfConjugateFlag( isSelfConjugate ),
      identifyingPdgCodes( pdgCodes.begin(),
                           pdgCodes.end() ),
      pdgCodeMap( mapAndVectorAndBools.getMap() ),
      mapFiller( 0,
                 NULL ),
      massRecorded( false ),
      defaultResetMass( defaultResetMass ),
      signedDefaultMass( defaultResetMass ),
      absoluteDefaultMass( defaultResetMass ),
      runningMasses(),
      runningMassesAsVector(),
      decaysRecorded( false ),
      defaultDecayWidth( defaultDecayWidth ),
      decayWidth( defaultDecayWidth ),
      defaultDecaySet(),
      decaySet(),
      decaySetAsVector(),
      asciiName( asciiName ),
      latexName( latexName ),
      isVerbose( mapAndVectorAndBools.getBool() ),
      flagBools( mapAndVectorAndBools.getFlags() ),
      setOfPointersOfMassEigenstateGroup( mapAndVectorAndBools.getVector() )
  {
    constructorBodyFunction();
    /* constructorBodyFunction puts pointers to this instance into pdgCodeMap
     * for each of this instance's positive codes. if this instance has
     * negative codes, they will only get into pdgCodeMap through the
     * charge-conjugate's constructor (or by manually adding them, but that's
     * not recommended).
     */
  }

  MassEigenstate::MassEigenstate( MassEigenstate const& copySource,
                         MassEigenstateMapVectorBools& mapAndVectorAndBools ) :
      chargeConjugate( NULL ),
      isSelfConjugateFlag( copySource.isSelfConjugateFlag ),
      identifyingPdgCodes( copySource.identifyingPdgCodes.begin(),
                           copySource.identifyingPdgCodes.end() ),
      pdgCodeMap( mapAndVectorAndBools.getMap() ),
      mapFiller( 0,
                 NULL ),
      massRecorded( false ),
      defaultResetMass( copySource.defaultResetMass ),
      signedDefaultMass( copySource.signedDefaultMass ),
      absoluteDefaultMass( copySource.absoluteDefaultMass ),
      runningMasses(),
      runningMassesAsVector(),
      decaysRecorded( copySource.decaysRecorded ),
      defaultDecayWidth( copySource.defaultDecayWidth ),
      decayWidth( copySource.decayWidth ),
      defaultDecaySet( copySource.defaultDecaySet ),
      decaySet( copySource.decaySet ),
      decaySetAsVector( copySource.decaySetAsVector ),
      asciiName( copySource.asciiName ),
      latexName( copySource.latexName ),
      isVerbose( mapAndVectorAndBools.getBool() ),
      flagBools( mapAndVectorAndBools.getFlags() ),
      setOfPointersOfMassEigenstateGroup( mapAndVectorAndBools.getVector() )
  {
    constructorBodyFunction();
    /* constructorBodyFunction puts pointers to this instance into pdgCodeMap
     * for each of this instance's positive codes. if this instance has
     * negative codes, they will only get into pdgCodeMap through the
     * charge-conjugate's constructor (or by manually adding them, but that's
     * not recommended).
     */
  }

  // this last version copies as the charge conjugate of copySource:
  MassEigenstate::MassEigenstate( MassEigenstate& copySource,
                                  std::string const& asciiName,
                                  std::string const& latexName ) :
      chargeConjugate( &copySource ),
      isSelfConjugateFlag( copySource.isSelfConjugateFlag ),
      identifyingPdgCodes( copySource.identifyingPdgCodes.size(),
                           -(copySource.getCode()) ),
      pdgCodeMap( copySource.pdgCodeMap ),
      mapFiller( 0,
                 NULL ),
      massRecorded( false ),
      defaultResetMass( copySource.defaultResetMass ),
      signedDefaultMass( copySource.signedDefaultMass ),
      absoluteDefaultMass( copySource.absoluteDefaultMass ),
      runningMasses(),
      runningMassesAsVector(),
      decaysRecorded( false ),
      defaultDecayWidth( copySource.defaultDecayWidth ),
      decayWidth( copySource.decayWidth ),
      defaultDecaySet(),
      decaySet(),
      decaySetAsVector(),
      asciiName( asciiName ),
      latexName( latexName ),
      isVerbose( copySource.isVerbose ),
      flagBools( copySource.flagBools ),
      setOfPointersOfMassEigenstateGroup(
                                copySource.setOfPointersOfMassEigenstateGroup )
  {
    // this for loop doesn't copy the front element because it was already
    // copied into all the entries when identifyingPdgCodes was initialized.
    for( int codeIndex( copySource.identifyingPdgCodes.size() - 1 );
         0 < codeIndex;
         --codeIndex )
    {
      identifyingPdgCodes[ codeIndex ]
      = -(copySource.identifyingPdgCodes[ codeIndex ]);
      // this gets the negatives of copySource's codes.
    }
    copySource.setChargeConjugate( this );
    constructorBodyFunction();
    /* constructorBodyFunction sorts out getting copySource's negative codes
     * into pdgCodeMap as the positive codes for this instance, with this
     * instance's pointer.
     */
  }

  MassEigenstate::~MassEigenstate()
  {
    // does nothing, since there was no dynamic allocation.
  }


  MassEigenstate&
  MassEigenstate::addCode( int extraCode )
  {
    if( isSelfConjugate() )
    {
      if( 0 > extraCode )
      {
        extraCode = -extraCode;
      }
      // self-conjugate states was forced into having only positive codes.
    }
    else if( NULL != chargeConjugate )
    {
      // the charge-conjugate should get the opposite-signed codes, &
      // pdgCodeMap should know which is mapped to by the positive code:
      chargeConjugate->identifyingPdgCodes.push_back( -extraCode );
      if( 0 > extraCode )
      {
        chargeConjugate->addToCodeMap( -extraCode );
      }
    }
    // extraCode can now be added to identifyingPdgCodes:
    identifyingPdgCodes.push_back( extraCode );
    // if the code is positive, pdgCodeMap should know that it maps to this
    // instance:
    if( 0 <= extraCode )
    {
      addToCodeMap( extraCode );
    }
    return *this;
  }

  void
  MassEigenstate::setToBeChargeConjugate(
                                        MassEigenstate* const chargeConjugate )
  {
    this->chargeConjugate = chargeConjugate;
    /* not only does the pointer to the charge conjugate need to be set, but
     * the codes need to be signed. chargeConjugate's codes are assumed to be
     * the right signs.
     */
    for( int codeIndex( identifyingPdgCodes.size() - 1 );
         0 <= codeIndex;
         --codeIndex )
    {
      MassEigenstateCodeMap::iterator codeFinder;
      if( chargeConjugate->hasCode( identifyingPdgCodes[ codeIndex ] ) )
        // if chargeConjugate should have the positive code...
      {
        // ensure that pdgCodeMap associates chargeConjugate with the positive
        // code:
        codeFinder = pdgCodeMap.find( identifyingPdgCodes[ codeIndex ] );
        if( pdgCodeMap.end() != codeFinder )
        {
          codeFinder->second = chargeConjugate;
        }
        identifyingPdgCodes[ codeIndex ] = -(identifyingPdgCodes[ codeIndex ]);
        // now this instance can set its code to be negative.
      }
    }
  }

  void
  MassEigenstate::recordMass( double const massValue,
                              double const minusUncertainty,
                              double const plusUncertainty,
                              int const schemeType,
                              double const evaluationScale )
  {
    runningMasses.newEnd().setValues( massValue,
                                      minusUncertainty,
                                      plusUncertainty,
                                      schemeType,
                                      evaluationScale );
    runningMassesAsVector.push_back( runningMasses.getPointer(
                                              runningMasses.getLastIndex() ) );
    // I decided to keep a separate std::vector in parallel rather than
    // constantly converting runningMasses with getAsVector(...).
    if( !massRecorded
        ||
        ( 0 == schemeType ) )
      // the pole mass becomes the default mass, unless there is no mass
      // already recorded.
    {
      signedDefaultMass = massValue;
      if( 0 > massValue )
      {
        absoluteDefaultMass = -massValue;
      }
      else
      {
        absoluteDefaultMass = massValue;
      }
      massRecorded = true;
    }
  }

  void
  MassEigenstate::constructorBodyFunction()
  {
    if( isSelfConjugate() )
    {
      chargeConjugate = this;
      for( int codeIndex( identifyingPdgCodes.size() - 1 );
           0 <= codeIndex;
           --codeIndex )
      {
        if( 0 > identifyingPdgCodes[ codeIndex ] )
        {
          identifyingPdgCodes[ codeIndex ] = -identifyingPdgCodes[ codeIndex ];
        }
        addToCodeMap( identifyingPdgCodes[ codeIndex ] );
        // self-conjugate particles are set to have only positive codes, &
        // the codes are mapped to this MassEigenstate.
      }
    }
    else
    {
      for( int codeIndex( identifyingPdgCodes.size() - 1 );
           0 <= codeIndex;
           --codeIndex )
      {
        if( 0 < identifyingPdgCodes[ codeIndex ] )
        {
          addToCodeMap( identifyingPdgCodes[ codeIndex ] );
          /* only this instance's positive codes are mapped to its pointer in
           * pdgCodeMap. if this instance has negative codes, they will only
           * get into pdgCodeMap through the charge-conjugate's constructor (or
           * by manually adding them, but that's not recommended).
           */
        }
      }
    }
    setOfPointersOfMassEigenstateGroup.push_back( this );
  }

}
