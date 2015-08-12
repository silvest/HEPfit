/*
 * StandardModel.cpp
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
  namespace MassSpectrumClass
  {
    StandardModel::StandardModel( bool const isVerbose,
                                  bool const neutrinosAreMajorana,
                                  std::vector< bool >* const defaultFlags ) :
        MassSpectrum( isVerbose,
                      defaultFlags ),
        neutralColorlessScalarOne( PDGIX::neutralColorlessScalarOne,
                                   PDGVII::neutralColorlessScalarOne,
                                   mapAndVectorAndBools,
                                   true,
                                   "h01",
                                   "$h_{1}^{0}$" ),
        positronOne( PDGIX::positronOne,
                     -PDGVII::electronOne,
                     mapAndVectorAndBools,
                     false,
                     "e1bar",
                     "${\\bar{e}}_{1}$",
                     PdgData::electronMass,
                     0.0 ),
        antipositronOne( positronOne,
                         "e1",
                         "$e_{1}$" ),
        positronTwo( PDGIX::positronTwo,
                     -PDGVII::electronTwo,
                     mapAndVectorAndBools,
                     false,
                     "e2bar",
                     "${\\bar{e}}_{2}$",
                     PdgData::muonMass,
                     0.0 ),
        antipositronTwo( positronTwo,
                         "e2",
                         "$e_{2}$" ),
        positronThree( PDGIX::positronThree,
                       -PDGVII::electronThree,
                       mapAndVectorAndBools,
                       false,
                       "e3bar",
                       "${\\bar{e}}_{3}$",
                       PdgData::tauLeptonMass,
                       0.0 ),
        antipositronThree( positronThree,
                           "e3",
                           "$e_{3}$" ),
        antineutrinoOne( PDGIX::antineutrinoOne,
                         -PDGVII::neutrinoOne,
                         mapAndVectorAndBools,
                         false,
                         "v1bar",
                         "${\\bar{{\\nu}}}_{1}$",
                         PdgData::electronNeutrinoMass,
                         0.0 ),
        neutrinoOne( antineutrinoOne,
                     "v1",
                     "${\\nu}_{1}$" ),
        antineutrinoTwo( PDGIX::antineutrinoTwo,
                         -PDGVII::neutrinoTwo,
                         mapAndVectorAndBools,
                         false,
                         "v2bar",
                         "${\\bar{{\\nu}}}_{2}$",
                         PdgData::muonNeutrinoMass,
                         0.0 ),
        neutrinoTwo( antineutrinoTwo,
                     "v2",
                     "${\\nu}_{2}$" ),
        antineutrinoThree( PDGIX::antineutrinoThree,
                           -PDGVII::neutrinoThree,
                           mapAndVectorAndBools,
                           false,
                           "v3bar",
                           "${\\bar{{\\nu}}}_{3}$",
                           PdgData::tauNeutrinoMass,
                           0.0 ),
        neutrinoThree( antineutrinoThree,
                       "v3",
                       "${\\nu}_{3}$" ),
        antidownOne( PDGIX::antidownOne,
                     -PDGVII::downOne,
                     mapAndVectorAndBools,
                     false,
                     "d1bar",
                     "${\\bar{d}}_{1}$",
                     PdgData::downMass,
                     0.0 ),
        downOne( antidownOne,
                 "d1",
                 "$d_{1}$" ),
        antidownTwo( PDGIX::antidownTwo,
                     -PDGVII::downTwo,
                     mapAndVectorAndBools,
                     false,
                     "d2bar",
                     "${\\bar{d}}_{2}$",
                     PdgData::strangeMass,
                     0.0 ),
        downTwo( antidownTwo,
                 "d2",
                 "$d_{2}$" ),
        antidownThree( PDGIX::antidownThree,
                       -PDGVII::downThree,
                       mapAndVectorAndBools,
                       false,
                       "d3bar",
                       "${\\bar{d}}_{3}$",
                       PdgData::bottomMass,
                       0.0 ),
        downThree( antidownThree,
                   "d3",
                   "$d_{3}$" ),
        upOne( PDGIX::upOne,
               PDGVII::upOne,
               mapAndVectorAndBools,
               false,
               "u1",
               "$u_{1}$",
               PdgData::upMass,
               0.0 ),
        antiupOne( upOne,
                   "u1bar",
                   "${\\bar{u}}_{1}$" ),
        upTwo( PDGIX::upTwo,
               PDGVII::upTwo,
               mapAndVectorAndBools,
               false,
               "u2",
               "$u_{2}$",
               PdgData::charmMass,
               0.0 ),
        antiupTwo( upTwo,
                   "u2bar",
                   "${\\bar{u}}_{2}$" ),
        upThree( PDGIX::upThree,
                 PDGVII::upThree,
                 mapAndVectorAndBools,
                 false,
                 "u3",
                 "$u_{3}$",
                 PdgData::topMass,
                 PdgData::topDecayWidth ),
        antiupThree( upThree,
                     "u3bar",
                     "${\\bar{u}}_{3}$" ),
        photonBoson( PDGIX::photonBoson,
                     PDGVII::photonBoson,
                     mapAndVectorAndBools,
                     true,
                     "A",
                     "${\\gamma}$",
                     PdgData::photonMass,
                     0.0 ),
        zBosonOne( PDGIX::zBosonOne,
                   PDGVII::zBosonOne,
                   mapAndVectorAndBools,
                   true,
                   "Z",
                   "$Z$",
                   PdgData::zMass,
                   PdgData::zDecayWidth ),
        wPlusBosonOne( PDGIX::wPlusBosonOne,
                       PDGVII::wPlusBosonOne,
                       mapAndVectorAndBools,
                       false,
                       "Wp",
                       "$W^{+}$",
                       PdgData::wPlusMass,
                       PdgData::wPlusDecayWidth ),
        wMinusBosonOne( wPlusBosonOne,
                        "Wm",
                        "$W^{-}$" ),
        gluonBoson( PDGIX::gluonBoson,
                    PDGVII::gluonBoson,
                    mapAndVectorAndBools,
                    true,
                    "g",
                    "$g$",
                    PdgData::gluonMass,
                    0.0 ),
        positiveLeptonPointers( 3,
                                &positronOne ),
        negativeLeptonPointers( 3,
                                &antipositronOne ),
        antineutrinoPointers( 3,
                              &antineutrinoOne ),
        neutrinoPointers( 3,
                          &neutrinoOne ),
        downAntiquarkPointers( 3,
                               &antidownOne ),
        downQuarkPointers( 3,
                           &downOne ),
        upQuarkPointers( 3,
                         &upOne ),
        upAntiquarkPointers( 3,
                             &antiupOne ),
        defaultDecayFiller()
    {
      positiveLeptonPointers[ 1 ] = &positronTwo;
      positiveLeptonPointers[ 2 ] = &positronThree;
      negativeLeptonPointers[ 1 ] = &antipositronTwo;
      negativeLeptonPointers[ 2 ] = &antipositronThree;
      antineutrinoPointers[ 1 ] = &antineutrinoTwo;
      antineutrinoPointers[ 2 ] = &antineutrinoThree;
      neutrinoPointers[ 1 ] = &neutrinoTwo;
      neutrinoPointers[ 2 ] = &neutrinoThree;
      downAntiquarkPointers[ 1 ] = &antidownTwo;
      downAntiquarkPointers[ 2 ] = &antidownThree;
      downQuarkPointers[ 1 ] = &downTwo;
      downQuarkPointers[ 2 ] = &downThree;
      upQuarkPointers[ 1 ] = &upTwo;
      upQuarkPointers[ 2 ] = &upThree;
      upAntiquarkPointers[ 1 ] = &antiupTwo;
      upAntiquarkPointers[ 2 ] = &antiupThree;
      if( neutrinosAreMajorana )
      {
        setNeutrinosToMajorana();
      }
      if( &defaultBoolVector == mapAndVectorAndBools.getFlags() )
      {
        positronOne.setFlags( &defaultIsLightLeptonBoolVector );
        antipositronOne.setFlags( &defaultIsLightLeptonBoolVector );
        positronTwo.setFlags( &defaultIsLightLeptonBoolVector );
        antipositronTwo.setFlags( &defaultIsLightLeptonBoolVector );
        for( int whichParticle( 2 );
             0 <= whichParticle;
             --whichParticle )
          // the number of generations is hard-coded in for the Standard Model.
        {
          antineutrinoPointers[ whichParticle ]->setFlags(
                                           &defaultEscapesDetectorBoolVector );
          neutrinoPointers[ whichParticle ]->setFlags(
                                           &defaultEscapesDetectorBoolVector );
          downAntiquarkPointers[ whichParticle ]->setFlags(
                                                     &defaultIsJetBoolVector );
          downQuarkPointers[ whichParticle ]->setFlags(
                                                     &defaultIsJetBoolVector );
          upQuarkPointers[ whichParticle ]->setFlags(
                                                     &defaultIsJetBoolVector );
          upAntiquarkPointers[ whichParticle ]->setFlags(
                                                     &defaultIsJetBoolVector );
        }
        gluonBoson.setFlags( &defaultIsJetBoolVector );
      }

      // now the PDG decays of the top quark & W & Z bosons are recorded.

      // top decays:
      defaultDecayFiller.addPointer( &wPlusBosonOne );
      defaultDecayFiller.addPointer( &downThree );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                                 PdgData::topToWPlusBottomBr );
      upThree.recordDecayAsDefault( defaultDecayFiller );
      antiupThree.recordChargeConjugateOfDecayAsDefault( defaultDecayFiller );

      // Z decays:
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &positronOne );
      defaultDecayFiller.addPointer( &antipositronOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                          PdgData::zToElectronAntielectronBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &positronTwo );
      defaultDecayFiller.addPointer( &antipositronTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                                  PdgData::zToMuonAntimuonBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &positronThree );
      defaultDecayFiller.addPointer( &antipositronThree );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                           PdgData::zToTauLeptonAntileptonBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &neutrinoOne );
      defaultDecayFiller.addPointer( &antineutrinoOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                  PdgData::zToElectronNeutrinoAntineutrinoBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &neutrinoTwo );
      defaultDecayFiller.addPointer( &antineutrinoTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                      PdgData::zToMuonNeutrinoAntineutrinoBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &neutrinoThree );
      defaultDecayFiller.addPointer( &antineutrinoThree );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                       PdgData::zToTauNeutrinoAntineutrinoBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &downOne );
      defaultDecayFiller.addPointer( &antidownOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                                  PdgData::zToDownAntidownBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &downTwo );
      defaultDecayFiller.addPointer( &antidownTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                            PdgData::zToStrangeAntistrangeBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &downThree );
      defaultDecayFiller.addPointer( &antidownThree );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                              PdgData::zToBottomAntibottomBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &upOne );
      defaultDecayFiller.addPointer( &antiupOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                                      PdgData::zToUpAntiupBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &upTwo );
      defaultDecayFiller.addPointer( &antiupTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                                PdgData::zToCharmAnticharmBr );
      zBosonOne.recordDecayAsDefault( defaultDecayFiller );

      // W decays:
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &positronOne );
      defaultDecayFiller.addPointer( &neutrinoOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                      PdgData::wPlusToNeutrinoAntielectronBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &positronTwo );
      defaultDecayFiller.addPointer( &neutrinoTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                          PdgData::wPlusToNeutrinoAntimuonBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &positronThree );
      defaultDecayFiller.addPointer( &neutrinoThree );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                     PdgData::wPlusToNeutrinoTauAntileptonBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &upOne );
      defaultDecayFiller.addPointer( &antidownOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                                PdgData::wPlusToUpAntidownBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &upOne );
      defaultDecayFiller.addPointer( &antidownTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                             PdgData::wPlusToUpAntistrangeBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &upTwo );
      defaultDecayFiller.addPointer( &antidownOne );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                             PdgData::wPlusToCharmAntidownBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
      defaultDecayFiller.clearPointers();
      defaultDecayFiller.addPointer( &upTwo );
      defaultDecayFiller.addPointer( &antidownTwo );
      defaultDecayFiller.setPairedValueAndSortPointers(
                                          PdgData::wPlusToCharmAntistrangeBr );
      wPlusBosonOne.recordDecayAsDefault( defaultDecayFiller );
      wMinusBosonOne.recordChargeConjugateOfDecayAsDefault(
                                                          defaultDecayFiller );
    }

    StandardModel::~StandardModel()
    {
      // does nothing.
    }

  }

}
