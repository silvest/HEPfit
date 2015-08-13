/*
 * SevenSlhaCodes.cpp
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

#include "PDG.hpp"

namespace LHPC
{
  // Standard Model (SM) particles:
  int const SevenDigitSlhaCodes::neutralColorlessScalarOne( 25 );
  int const SevenDigitSlhaCodes::higgsBoson( neutralColorlessScalarOne );
  int const SevenDigitSlhaCodes::positronOne( -11 );
  int const SevenDigitSlhaCodes::electronOne( -positronOne );
  int const SevenDigitSlhaCodes::positiveElectron( positronOne );
  int const SevenDigitSlhaCodes::negativeElectron( -positronOne );
  int const SevenDigitSlhaCodes::positronTwo( -13 );
  int const SevenDigitSlhaCodes::electronTwo( -positronTwo );
  int const SevenDigitSlhaCodes::positiveMuon( positronTwo );
  int const SevenDigitSlhaCodes::negativeMuon( -positronTwo );
  int const SevenDigitSlhaCodes::positronThree( -15 );
  int const SevenDigitSlhaCodes::electronThree( -positronThree );
  int const SevenDigitSlhaCodes::positiveTau( positronThree );
  int const SevenDigitSlhaCodes::negativeTau( -positronThree );
  int const SevenDigitSlhaCodes::antineutrinoOne( -12 );
  int const SevenDigitSlhaCodes::neutrinoOne( -antineutrinoOne );
  int const SevenDigitSlhaCodes::electronAntineutrino( antineutrinoOne );
  int const SevenDigitSlhaCodes::electronNeutrino( -antineutrinoOne );
  int const SevenDigitSlhaCodes::antineutrinoTwo( -14 );
  int const SevenDigitSlhaCodes::neutrinoTwo( -antineutrinoTwo );
  int const SevenDigitSlhaCodes::muonAntineutrino( antineutrinoTwo );
  int const SevenDigitSlhaCodes::muonNeutrino( -antineutrinoTwo );
  int const SevenDigitSlhaCodes::antineutrinoThree( -16 );
  int const SevenDigitSlhaCodes::neutrinoThree( -antineutrinoThree );
  int const SevenDigitSlhaCodes::tauAntineutrino( antineutrinoThree );
  int const SevenDigitSlhaCodes::tauNeutrino( -antineutrinoThree );
  int const SevenDigitSlhaCodes::neutrinoOneMajorana( -antineutrinoOne );
  int const SevenDigitSlhaCodes::neutrinoTwoMajorana( -antineutrinoTwo );
  int const SevenDigitSlhaCodes::neutrinoThreeMajorana( -antineutrinoThree );
  int const SevenDigitSlhaCodes::antidownOne( -1 );
  int const SevenDigitSlhaCodes::downOne( -antidownOne );
  int const SevenDigitSlhaCodes::downAntiquark( antidownOne );
  int const SevenDigitSlhaCodes::downQuark( -antidownOne );
  int const SevenDigitSlhaCodes::antidownTwo( -3 );
  int const SevenDigitSlhaCodes::downTwo( -antidownTwo );
  int const SevenDigitSlhaCodes::strangeAntiquark( antidownTwo );
  int const SevenDigitSlhaCodes::strangeQuark( -antidownTwo );
  int const SevenDigitSlhaCodes::antidownThree( -5 );
  int const SevenDigitSlhaCodes::downThree( -antidownThree );
  int const SevenDigitSlhaCodes::bottomAntiquark( antidownThree );
  int const SevenDigitSlhaCodes::bottomQuark( -antidownThree );
  int const SevenDigitSlhaCodes::upOne( 2 );
  int const SevenDigitSlhaCodes::antiupOne( -upOne );
  int const SevenDigitSlhaCodes::upQuark( upOne );
  int const SevenDigitSlhaCodes::upAntiquark( -upOne );
  int const SevenDigitSlhaCodes::upTwo( 4 );
  int const SevenDigitSlhaCodes::antiupTwo( -upTwo );
  int const SevenDigitSlhaCodes::charmQuark( upTwo );
  int const SevenDigitSlhaCodes::charmAntiquark( -upTwo );
  int const SevenDigitSlhaCodes::upThree( 6 );
  int const SevenDigitSlhaCodes::antiupThree( -upThree );
  int const SevenDigitSlhaCodes::topQuark( upThree );
  int const SevenDigitSlhaCodes::topAntiquark( -upThree );
  int const SevenDigitSlhaCodes::photonBoson( 22 );
  int const SevenDigitSlhaCodes::zBosonOne( 23 );
  int const SevenDigitSlhaCodes::zBoson( zBosonOne );
  int const SevenDigitSlhaCodes::wPlusBosonOne( 24 );
  int const SevenDigitSlhaCodes::wPlus( wPlusBosonOne );
  int const SevenDigitSlhaCodes::wMinus( -wPlusBosonOne );
  int const SevenDigitSlhaCodes::gluonBoson( 21 );

  // almost-SM particles:
  int const SevenDigitSlhaCodes::gravitonBoson( 39 );

  // Minimal Supersymmetric Standard Model (MSSM) particles except those in the
  // SM above, assuming R-parity conservation:
  int const SevenDigitSlhaCodes::higgsScalarOne( neutralColorlessScalarOne );
  int const SevenDigitSlhaCodes::lightHiggs( neutralColorlessScalarOne );
  int const SevenDigitSlhaCodes::neutralColorlessScalarTwo( 35 );
  int const SevenDigitSlhaCodes::higgsScalarTwo( neutralColorlessScalarTwo );
  int const SevenDigitSlhaCodes::heavyHiggs( neutralColorlessScalarTwo );
  int const SevenDigitSlhaCodes::neutralColorlessPseudoscalarOne( 36 );
  int const
  SevenDigitSlhaCodes::higgsPseudoscalar( neutralColorlessPseudoscalarOne );
  int const SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonOne( 37 );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonOne(
                                          -positiveColorlessSpinZeroBosonOne );
  int const SevenDigitSlhaCodes::positiveHiggs(
                                           positiveColorlessSpinZeroBosonOne );
  int const SevenDigitSlhaCodes::negativeHiggs(
                                          -positiveColorlessSpinZeroBosonOne );
  int const SevenDigitSlhaCodes::spositronOne( -1000011 );
  int const SevenDigitSlhaCodes::selectronOne( -spositronOne );
  int const SevenDigitSlhaCodes::antiselectronL( spositronOne );
  int const SevenDigitSlhaCodes::selectronL( -spositronOne );
  int const SevenDigitSlhaCodes::positiveSelectronL( spositronOne );
  int const SevenDigitSlhaCodes::negativeSelectronL( -spositronOne );
  int const SevenDigitSlhaCodes::spositronTwo( -1000013 );
  int const SevenDigitSlhaCodes::selectronTwo( -spositronTwo );
  int const SevenDigitSlhaCodes::antismuonL( spositronTwo );
  int const SevenDigitSlhaCodes::smuonL( -spositronTwo );
  int const SevenDigitSlhaCodes::positiveSmuonL( spositronTwo );
  int const SevenDigitSlhaCodes::negativeSmuonL( -spositronTwo );
  int const SevenDigitSlhaCodes::spositronThree( -1000015 );
  int const SevenDigitSlhaCodes::selectronThree( -spositronThree );
  int const SevenDigitSlhaCodes::antistauOne( spositronThree );
  int const SevenDigitSlhaCodes::stauOne( -spositronThree );
  int const SevenDigitSlhaCodes::positiveStauOne( spositronThree );
  int const SevenDigitSlhaCodes::negativeStauOne( -spositronThree );
  int const SevenDigitSlhaCodes::spositronFour( -2000011 );
  int const SevenDigitSlhaCodes::selectronFour( -spositronFour );
  int const SevenDigitSlhaCodes::antiselectronR( spositronFour );
  int const SevenDigitSlhaCodes::selectronR( -spositronFour );
  int const SevenDigitSlhaCodes::positiveSelectronR( spositronFour );
  int const SevenDigitSlhaCodes::negativeSelectronR( -spositronFour );
  int const SevenDigitSlhaCodes::spositronFive( -2000013 );
  int const SevenDigitSlhaCodes::selectronFive( -spositronFive );
  int const SevenDigitSlhaCodes::antismuonR( spositronFive );
  int const SevenDigitSlhaCodes::smuonR( -spositronFive );
  int const SevenDigitSlhaCodes::positiveSmuonR( spositronFive );
  int const SevenDigitSlhaCodes::negativeSmuonR( -spositronFive );
  int const SevenDigitSlhaCodes::spositronSix( -2000015 );
  int const SevenDigitSlhaCodes::selectronSix( -spositronSix );
  int const SevenDigitSlhaCodes::antistauTwo( spositronSix );
  int const SevenDigitSlhaCodes::stauTwo( -spositronSix );
  int const SevenDigitSlhaCodes::positiveStauTwo( spositronSix );
  int const SevenDigitSlhaCodes::negativeStauTwo( -spositronSix );
  int const SevenDigitSlhaCodes::antisneutrinoOne( -1000012 );
  int const SevenDigitSlhaCodes::sneutrinoOne( -antisneutrinoOne );
  int const SevenDigitSlhaCodes::electronAntisneutrinoL( antisneutrinoOne );
  int const SevenDigitSlhaCodes::electronSneutrinoL( -antisneutrinoOne );
  int const SevenDigitSlhaCodes::antisneutrinoTwo( -1000014 );
  int const SevenDigitSlhaCodes::sneutrinoTwo( -antisneutrinoTwo );
  int const SevenDigitSlhaCodes::muonAntisneutrinoL( antisneutrinoTwo );
  int const SevenDigitSlhaCodes::muonSneutrinoL( -antisneutrinoTwo );
  int const SevenDigitSlhaCodes::antisneutrinoThree( -1000016 );
  int const SevenDigitSlhaCodes::sneutrinoThree( -antisneutrinoThree );
  int const SevenDigitSlhaCodes::tauAntisneutrinoL( antisneutrinoThree );
  int const SevenDigitSlhaCodes::tauSneutrinoL( -antisneutrinoThree );
  int const SevenDigitSlhaCodes::sneutrinoScalarOne( 1000012 );
  int const SevenDigitSlhaCodes::sneutrinoScalarTwo( 1000014 );
  int const SevenDigitSlhaCodes::sneutrinoScalarThree( 1000016 );
  int const SevenDigitSlhaCodes::sneutrinoPseudoscalarOne( 1000017 );
  int const SevenDigitSlhaCodes::sneutrinoPseudoscalarTwo( 1000018 );
  int const SevenDigitSlhaCodes::sneutrinoPseudoscalarThree( 1000019 );
  int const SevenDigitSlhaCodes::antisdownOne( -1000001 );
  int const SevenDigitSlhaCodes::sdownOne( -antisdownOne );
  int const SevenDigitSlhaCodes::antisdownL( antisdownOne );
  int const SevenDigitSlhaCodes::sdownL( -antisdownOne );
  int const SevenDigitSlhaCodes::antisdownTwo( -1000003 );
  int const SevenDigitSlhaCodes::sdownTwo( -antisdownTwo );
  int const SevenDigitSlhaCodes::antisstrangeL( antisdownTwo );
  int const SevenDigitSlhaCodes::sstrangeL( -antisdownTwo );
  int const SevenDigitSlhaCodes::antisdownThree( -1000005 );
  int const SevenDigitSlhaCodes::sdownThree( -antisdownThree );
  int const SevenDigitSlhaCodes::antisbottomOne( antisdownThree );
  int const SevenDigitSlhaCodes::sbottomOne( -antisdownThree );
  int const SevenDigitSlhaCodes::antisdownFour( -2000001 );
  int const SevenDigitSlhaCodes::sdownFour( -antisdownFour );
  int const SevenDigitSlhaCodes::antisdownR( antisdownFour );
  int const SevenDigitSlhaCodes::sdownR( -antisdownFour );
  int const SevenDigitSlhaCodes::antisdownFive( -2000003 );
  int const SevenDigitSlhaCodes::sdownFive( -antisdownFive );
  int const SevenDigitSlhaCodes::antisstrangeR( antisdownFive );
  int const SevenDigitSlhaCodes::sstrangeR( -antisdownFive );
  int const SevenDigitSlhaCodes::antisdownSix( -2000005 );
  int const SevenDigitSlhaCodes::sdownSix( -antisdownSix );
  int const SevenDigitSlhaCodes::antisbottomTwo( antisdownSix );
  int const SevenDigitSlhaCodes::sbottomTwo( -antisdownSix );
  int const SevenDigitSlhaCodes::supOne( 1000002 );
  int const SevenDigitSlhaCodes::antisupOne( -supOne );
  int const SevenDigitSlhaCodes::supL( supOne );
  int const SevenDigitSlhaCodes::antisupL( -supOne );
  int const SevenDigitSlhaCodes::supTwo( 1000004 );
  int const SevenDigitSlhaCodes::antisupTwo( -supTwo );
  int const SevenDigitSlhaCodes::scharmL( supTwo );
  int const SevenDigitSlhaCodes::antischarmL( -supTwo );
  int const SevenDigitSlhaCodes::supThree( 1000006 );
  int const SevenDigitSlhaCodes::antisupThree( -supThree );
  int const SevenDigitSlhaCodes::stopOne( supThree );
  int const SevenDigitSlhaCodes::antistopOne( -supThree );
  int const SevenDigitSlhaCodes::supFour( 2000002 );
  int const SevenDigitSlhaCodes::antisupFour( -supFour );
  int const SevenDigitSlhaCodes::supR( supFour );
  int const SevenDigitSlhaCodes::antisupR( -supFour );
  int const SevenDigitSlhaCodes::supFive( 2000004 );
  int const SevenDigitSlhaCodes::antisupFive( -supFive );
  int const SevenDigitSlhaCodes::scharmR( supFive );
  int const SevenDigitSlhaCodes::antischarmR( -supFive );
  int const SevenDigitSlhaCodes::supSix( 2000006 );
  int const SevenDigitSlhaCodes::antisupSix( -supSix );
  int const SevenDigitSlhaCodes::stopTwo( supSix );
  int const SevenDigitSlhaCodes::antistopTwo( -supSix );
  int const SevenDigitSlhaCodes::neutralinoOne( 1000022 );
  int const SevenDigitSlhaCodes::neutralinoTwo( 1000023 );
  int const SevenDigitSlhaCodes::neutralinoThree( 1000025 );
  int const SevenDigitSlhaCodes::neutralinoFour( 1000035 );
  int const SevenDigitSlhaCodes::positiveCharginoOne( 1000024 );
  int const SevenDigitSlhaCodes::negativeCharginoOne( -positiveCharginoOne );
  int const SevenDigitSlhaCodes::positiveCharginoTwo( 1000037 );
  int const SevenDigitSlhaCodes::negativeCharginoTwo( -positiveCharginoTwo );
  int const SevenDigitSlhaCodes::gluinoFermion( 1000021 );

  // extra MSSM particles without parities (R, CP):
  int const SevenDigitSlhaCodes::cpvHiggsOne( neutralColorlessScalarOne );
  int const SevenDigitSlhaCodes::cpvHiggsTwo( neutralColorlessScalarTwo );
  int const
  SevenDigitSlhaCodes::cpvHiggsThree( neutralColorlessPseudoscalarOne );
  int const
  SevenDigitSlhaCodes::rpvHiggsScalarOne( neutralColorlessScalarOne );
  int const
  SevenDigitSlhaCodes::rpvHiggsScalarTwo( neutralColorlessScalarTwo );
  int const SevenDigitSlhaCodes::rpvHiggsScalarThree( sneutrinoOne );
  int const SevenDigitSlhaCodes::rpvHiggsScalarFour( sneutrinoTwo );
  int const SevenDigitSlhaCodes::rpvHiggsScalarFive( sneutrinoThree );
  int const SevenDigitSlhaCodes::rpvHiggsPseudoscalarOne(
                                             neutralColorlessPseudoscalarOne );
  int const
  SevenDigitSlhaCodes::rpvHiggsPseudoscalarTwo( sneutrinoPseudoscalarOne );
  int const
  SevenDigitSlhaCodes::rpvHiggsPseudoscalarThree( sneutrinoPseudoscalarTwo );
  int const
  SevenDigitSlhaCodes::rpvHiggsPseudoscalarFour( sneutrinoPseudoscalarThree );
  int const
  SevenDigitSlhaCodes::positiveHiggsOne( positiveColorlessSpinZeroBosonOne );
  int const
  SevenDigitSlhaCodes::negativeHiggsOne( -positiveColorlessSpinZeroBosonOne );
  int const
  SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonTwo( spositronOne );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonTwo(
                                          -positiveColorlessSpinZeroBosonTwo );
  int const
  SevenDigitSlhaCodes::positiveHiggsTwo( positiveColorlessSpinZeroBosonTwo );
  int const
  SevenDigitSlhaCodes::negativeHiggsTwo( -positiveColorlessSpinZeroBosonTwo );
  int const
  SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonThree( spositronTwo );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonThree(
                                        -positiveColorlessSpinZeroBosonThree );
  int const SevenDigitSlhaCodes::positiveHiggsThree(
                                         positiveColorlessSpinZeroBosonThree );
  int const SevenDigitSlhaCodes::negativeHiggsThree(
                                        -positiveColorlessSpinZeroBosonThree );
  int const
  SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonFour( spositronThree );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonFour(
                                         -positiveColorlessSpinZeroBosonFour );
  int const
  SevenDigitSlhaCodes::positiveHiggsFour( positiveColorlessSpinZeroBosonFour );
  int const
  SevenDigitSlhaCodes::negativeHiggsFour( -positiveColorlessSpinZeroBosonFour );
  int const
  SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonFive( spositronFour );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonFive(
                                         -positiveColorlessSpinZeroBosonFive );
  int const
  SevenDigitSlhaCodes::positiveHiggsFive( positiveColorlessSpinZeroBosonFive );
  int const
  SevenDigitSlhaCodes::negativeHiggsFive( -positiveColorlessSpinZeroBosonFive );
  int const
  SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonSix( spositronFive );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonSix(
                                          -positiveColorlessSpinZeroBosonSix );
  int const
  SevenDigitSlhaCodes::positiveHiggsSix( positiveColorlessSpinZeroBosonSix );
  int const
  SevenDigitSlhaCodes::negativeHiggsSix( -positiveColorlessSpinZeroBosonSix );
  int const
  SevenDigitSlhaCodes::positiveColorlessSpinZeroBosonSeven( spositronSix );
  int const SevenDigitSlhaCodes::negativeColorlessSpinZeroBosonSeven(
                                        -negativeColorlessSpinZeroBosonSeven );
  int const SevenDigitSlhaCodes::positiveHiggsSeven(
                                         negativeColorlessSpinZeroBosonSeven );
  int const SevenDigitSlhaCodes::negativeHiggsSeven(
                                        -negativeColorlessSpinZeroBosonSeven );
  int const SevenDigitSlhaCodes::neutrinoFour( neutralinoOne );
  int const SevenDigitSlhaCodes::neutrinoFive( neutralinoTwo );
  int const SevenDigitSlhaCodes::neutrinoSix( neutralinoThree );
  int const SevenDigitSlhaCodes::neutrinoSeven( neutralinoFour );
  int const SevenDigitSlhaCodes::positronFour( positiveCharginoOne );
  int const SevenDigitSlhaCodes::electronFour( -positronFour );
  int const SevenDigitSlhaCodes::positiveElectronFour( positronFour );
  int const SevenDigitSlhaCodes::negativeElectronFour( -positronFour );
  int const SevenDigitSlhaCodes::positronFive( positiveCharginoTwo );
  int const SevenDigitSlhaCodes::electronFive( -positronFive );
  int const SevenDigitSlhaCodes::positiveElectronFive( positronFive );
  int const SevenDigitSlhaCodes::negativeElectronFive( -positronFive );

  // Next-to-Minimal Supersymmetric Standard Model (NMSSM) particles except
  // those in the MSSM above, with & without parities:
  int const
  SevenDigitSlhaCodes::nmssmHiggsScalarOne( neutralColorlessScalarOne );
  int const
  SevenDigitSlhaCodes::nmssmHiggsScalarTwo( neutralColorlessScalarTwo );
  int const SevenDigitSlhaCodes::nmssmHiggsScalarThree( 45 );
  int const SevenDigitSlhaCodes::nmssmHiggsPseudoscalarOne(
                                             neutralColorlessPseudoscalarOne );
  int const SevenDigitSlhaCodes::nmssmHiggsPseudoscalarTwo( 46 );
  int const SevenDigitSlhaCodes::neutralinoFive( 1000045 );
  int const SevenDigitSlhaCodes::neutrinoEight( neutralinoFive );

}
