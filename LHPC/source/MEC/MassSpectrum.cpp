/*
 * MassSpectrum.cpp
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
#include "BOLlib/include/BOLlib.hpp"

namespace LHPC
{
  std::vector< bool > const
  MassSpectrum::defaultBoolVector( (unsigned int)(MassSpectrum::sizeOfEnum),
                                   false );
  std::vector< bool > const
  MassSpectrum::defaultEscapesDetectorBoolVector(
                  BOL::StdVectorFiller< bool >( true )( false ).end( false ) );
  std::vector< bool > const
  MassSpectrum::defaultIsJetBoolVector(
                  BOL::StdVectorFiller< bool >( false )( true ).end( false ) );
  std::vector< bool > const
  MassSpectrum::defaultIsLightLeptonBoolVector(
                  BOL::StdVectorFiller< bool >( false )( false ).end( true ) );

  MassSpectrum::MassSpectrum( bool const isVerbose,
                              std::vector< bool > const* defaultFlags ) :
      allMassEigenstates(),
      unknownMassEigenstates(),
      pdgCodeMap(),
      isVerbose( isVerbose ),
      mapAndVectorAndBools( pdgCodeMap,
                            allMassEigenstates,
                            isVerbose )
  {
    if( NULL == defaultFlags )
    {
      defaultFlags = &defaultBoolVector;
    }
    mapAndVectorAndBools.withBools( defaultFlags );
  }

  MassSpectrum::~MassSpectrum()
  {
    for( int deletionIndex( unknownMassEigenstates.size() - 1 );
         0 <= deletionIndex;
         --deletionIndex )
    {
      delete unknownMassEigenstates[ deletionIndex ];
    }
  }


  MassEigenstate&
  MassSpectrum::ensureMassEigenstateExists( int const pdgCode )
  /* this looks for a MassEigenstate with code pdgCode, & if it doesn't find
   * one, it creates a new MassEigenstate instance with a name that is just
   * pdgCode as a string, & also its charge conjugate, & notes their pointers
   * in unknownMassEigenstates so that they can be deleted by this MassSpectrum
   * instance's destructor.
   */
  {
    MassEigenstate*
    massEigenstatePointer( MassEigenstate::findPointerWithCode( pdgCode,
                                                                pdgCodeMap ) );
    if( NULL == massEigenstatePointer )
      // if the requested MassEigenstate could not be found...
    {
      std::string codeAsName( BOL::StringParser::intToString( pdgCode,
                                                              1 ) );
      // the code is turned into a name.
      massEigenstatePointer = new MassEigenstate( pdgCode,
                          mapAndVectorAndBools.withBools( &defaultBoolVector ),
                                                  false,
                                                  codeAsName,
                                                  codeAsName );
      // the MassEigenstate is constructed with default flags.
      unknownMassEigenstates.push_back( massEigenstatePointer );
      // the pointer is noted so that it can be deleted when appropriate.
      codeAsName.assign( BOL::StringParser::intToString( -pdgCode,
                                                         1 ) );
      // the charge conjugate's name is made, & then used for the creation of
      // the charge conjugate itself, while storing its pointer for deletion:
      unknownMassEigenstates.push_back( new MassEigenstate(
                                                        *massEigenstatePointer,
                                                            codeAsName,
                                                            codeAsName ) );
    }
    return *massEigenstatePointer;
  }

}
