/*
 * NineDigitSlhaCodes.cpp
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
  int const NineDigitSlhaCodes::neutralColorlessScalarOne( 101000001 );
  int const NineDigitSlhaCodes::higgsBoson( neutralColorlessScalarOne );
  int const NineDigitSlhaCodes::positronOne( 110000601 );
  int const NineDigitSlhaCodes::electronOne( -positronOne );
  int const NineDigitSlhaCodes::positiveElectron( positronOne );
  int const NineDigitSlhaCodes::negativeElectron( -positronOne );
  int const NineDigitSlhaCodes::positronTwo( 110000602 );
  int const NineDigitSlhaCodes::electronTwo( -positronTwo );
  int const NineDigitSlhaCodes::positiveMuon( positronTwo );
  int const NineDigitSlhaCodes::negativeMuon( -positronTwo );
  int const NineDigitSlhaCodes::positronThree( 110000603 );
  int const NineDigitSlhaCodes::electronThree( -positronThree );
  int const NineDigitSlhaCodes::positiveTau( positronThree );
  int const NineDigitSlhaCodes::negativeTau( -positronThree );
  int const NineDigitSlhaCodes::antineutrinoOne( 110000001 );
  int const NineDigitSlhaCodes::neutrinoOne( -antineutrinoOne );
  int const NineDigitSlhaCodes::electronAntineutrino( antineutrinoOne );
  int const NineDigitSlhaCodes::electronNeutrino( -antineutrinoOne );
  int const NineDigitSlhaCodes::antineutrinoTwo( 110000002 );
  int const NineDigitSlhaCodes::neutrinoTwo( -antineutrinoTwo );
  int const NineDigitSlhaCodes::muonAntineutrino( antineutrinoTwo );
  int const NineDigitSlhaCodes::muonNeutrino( -antineutrinoTwo );
  int const NineDigitSlhaCodes::antineutrinoThree( 110000003 );
  int const NineDigitSlhaCodes::neutrinoThree( -antineutrinoThree );
  int const NineDigitSlhaCodes::tauAntineutrino( antineutrinoThree );
  int const NineDigitSlhaCodes::tauNeutrino( -antineutrinoThree );
  int const NineDigitSlhaCodes::neutrinoOneMajorana( 111000001 );
  int const NineDigitSlhaCodes::neutrinoTwoMajorana( 111000002 );
  int const NineDigitSlhaCodes::neutrinoThreeMajorana( 111000003 );
  int const NineDigitSlhaCodes::antidownOne( 110890201 );
  int const NineDigitSlhaCodes::downOne( -antidownOne );
  int const NineDigitSlhaCodes::downAntiquark( antidownOne );
  int const NineDigitSlhaCodes::downQuark( -antidownOne );
  int const NineDigitSlhaCodes::antidownTwo( 110890202 );
  int const NineDigitSlhaCodes::downTwo( -antidownTwo );
  int const NineDigitSlhaCodes::strangeAntiquark( antidownTwo );
  int const NineDigitSlhaCodes::strangeQuark( -antidownTwo );
  int const NineDigitSlhaCodes::antidownThree( 110890203 );
  int const NineDigitSlhaCodes::downThree( -antidownThree );
  int const NineDigitSlhaCodes::bottomAntiquark( antidownThree );
  int const NineDigitSlhaCodes::bottomQuark( -antidownThree );
  int const NineDigitSlhaCodes::upOne( 110100401 );
  int const NineDigitSlhaCodes::antiupOne( -upOne );
  int const NineDigitSlhaCodes::upQuark( upOne );
  int const NineDigitSlhaCodes::upAntiquark( -upOne );
  int const NineDigitSlhaCodes::upTwo( 110100402 );
  int const NineDigitSlhaCodes::antiupTwo( -upTwo );
  int const NineDigitSlhaCodes::charmQuark( upTwo );
  int const NineDigitSlhaCodes::charmAntiquark( -upTwo );
  int const NineDigitSlhaCodes::upThree( 110100403 );
  int const NineDigitSlhaCodes::antiupThree( -upThree );
  int const NineDigitSlhaCodes::topQuark( upThree );
  int const NineDigitSlhaCodes::topAntiquark( -upThree );
  int const NineDigitSlhaCodes::photonBoson( 121000001 );
  int const NineDigitSlhaCodes::zBosonOne( 122000001 );
  int const NineDigitSlhaCodes::zBoson( zBosonOne );
  int const NineDigitSlhaCodes::wPlusBosonOne( 120000601 );
  int const NineDigitSlhaCodes::wPlus( wPlusBosonOne );
  int const NineDigitSlhaCodes::wMinus( -wPlusBosonOne );
  int const NineDigitSlhaCodes::gluonBoson( 121110001 );

  // almost-SM particles:
  int const NineDigitSlhaCodes::gravitonBoson( 141000001 );

  // Minimal Supersymmetric Standard Model (MSSM) particles except those in the
  // SM above, assuming R-parity conservation:
  int const NineDigitSlhaCodes::higgsScalarOne( neutralColorlessScalarOne );
  int const NineDigitSlhaCodes::lightHiggs( neutralColorlessScalarOne );
  int const NineDigitSlhaCodes::neutralColorlessScalarTwo( 101000002 );
  int const NineDigitSlhaCodes::higgsScalarTwo( neutralColorlessScalarTwo );
  int const NineDigitSlhaCodes::heavyHiggs( neutralColorlessScalarTwo );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarOne( 102000001 );
  int const
  NineDigitSlhaCodes::higgsPseudoscalar( neutralColorlessPseudoscalarOne );
  int const NineDigitSlhaCodes::positiveColorlessSpinZeroBosonOne( 100000601 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonOne(
                                          -positiveColorlessSpinZeroBosonOne );
  int const NineDigitSlhaCodes::positiveHiggs(
                                           positiveColorlessSpinZeroBosonOne );
  int const NineDigitSlhaCodes::negativeHiggs(
                                          -positiveColorlessSpinZeroBosonOne );
  int const NineDigitSlhaCodes::spositronOne( 200000601 );
  int const NineDigitSlhaCodes::selectronOne( -spositronOne );
  int const NineDigitSlhaCodes::antiselectronL( spositronOne );
  int const NineDigitSlhaCodes::selectronL( -spositronOne );
  int const NineDigitSlhaCodes::positiveSelectronL( spositronOne );
  int const NineDigitSlhaCodes::negativeSelectronL( -spositronOne );
  int const NineDigitSlhaCodes::spositronTwo( 200000602 );
  int const NineDigitSlhaCodes::selectronTwo( -spositronTwo );
  int const NineDigitSlhaCodes::antismuonL( spositronTwo );
  int const NineDigitSlhaCodes::smuonL( -spositronTwo );
  int const NineDigitSlhaCodes::positiveSmuonL( spositronTwo );
  int const NineDigitSlhaCodes::negativeSmuonL( -spositronTwo );
  int const NineDigitSlhaCodes::spositronThree( 200000603 );
  int const NineDigitSlhaCodes::selectronThree( -spositronThree );
  int const NineDigitSlhaCodes::antistauOne( spositronThree );
  int const NineDigitSlhaCodes::stauOne( -spositronThree );
  int const NineDigitSlhaCodes::positiveStauOne( spositronThree );
  int const NineDigitSlhaCodes::negativeStauOne( -spositronThree );
  int const NineDigitSlhaCodes::spositronFour( 200000604 );
  int const NineDigitSlhaCodes::selectronFour( -spositronFour );
  int const NineDigitSlhaCodes::antiselectronR( spositronFour );
  int const NineDigitSlhaCodes::selectronR( -spositronFour );
  int const NineDigitSlhaCodes::positiveSelectronR( spositronFour );
  int const NineDigitSlhaCodes::negativeSelectronR( -spositronFour );
  int const NineDigitSlhaCodes::spositronFive( 200000605 );
  int const NineDigitSlhaCodes::selectronFive( -spositronFive );
  int const NineDigitSlhaCodes::antismuonR( spositronFive );
  int const NineDigitSlhaCodes::smuonR( -spositronFive );
  int const NineDigitSlhaCodes::positiveSmuonR( spositronFive );
  int const NineDigitSlhaCodes::negativeSmuonR( -spositronFive );
  int const NineDigitSlhaCodes::spositronSix( 200000606 );
  int const NineDigitSlhaCodes::selectronSix( -spositronSix );
  int const NineDigitSlhaCodes::antistauTwo( spositronSix );
  int const NineDigitSlhaCodes::stauTwo( -spositronSix );
  int const NineDigitSlhaCodes::positiveStauTwo( spositronSix );
  int const NineDigitSlhaCodes::negativeStauTwo( -spositronSix );
  int const NineDigitSlhaCodes::antisneutrinoOne( 200000001 );
  int const NineDigitSlhaCodes::sneutrinoOne( -antisneutrinoOne );
  int const NineDigitSlhaCodes::electronAntisneutrinoL( antisneutrinoOne );
  int const NineDigitSlhaCodes::electronSneutrinoL( -antisneutrinoOne );
  int const NineDigitSlhaCodes::antisneutrinoTwo( 200000002 );
  int const NineDigitSlhaCodes::sneutrinoTwo( -antisneutrinoTwo );
  int const NineDigitSlhaCodes::muonAntisneutrinoL( antisneutrinoTwo );
  int const NineDigitSlhaCodes::muonSneutrinoL( -antisneutrinoTwo );
  int const NineDigitSlhaCodes::antisneutrinoThree( 200000003 );
  int const NineDigitSlhaCodes::sneutrinoThree( -antisneutrinoThree );
  int const NineDigitSlhaCodes::tauAntisneutrinoL( antisneutrinoThree );
  int const NineDigitSlhaCodes::tauSneutrinoL( -antisneutrinoThree );
  int const NineDigitSlhaCodes::sneutrinoScalarOne( 201000001 );
  int const NineDigitSlhaCodes::sneutrinoScalarTwo( 201000002 );
  int const NineDigitSlhaCodes::sneutrinoScalarThree( 201000003 );
  int const NineDigitSlhaCodes::sneutrinoPseudoscalarOne( 202000001 );
  int const NineDigitSlhaCodes::sneutrinoPseudoscalarTwo( 202000002 );
  int const NineDigitSlhaCodes::sneutrinoPseudoscalarThree( 202000003 );
  int const NineDigitSlhaCodes::antisdownOne( 200890201 );
  int const NineDigitSlhaCodes::sdownOne( -antisdownOne );
  int const NineDigitSlhaCodes::antisdownL( antisdownOne );
  int const NineDigitSlhaCodes::sdownL( -antisdownOne );
  int const NineDigitSlhaCodes::antisdownTwo( 200890202 );
  int const NineDigitSlhaCodes::sdownTwo( -antisdownTwo );
  int const NineDigitSlhaCodes::antisstrangeL( antisdownTwo );
  int const NineDigitSlhaCodes::sstrangeL( -antisdownTwo );
  int const NineDigitSlhaCodes::antisdownThree( 200890203 );
  int const NineDigitSlhaCodes::sdownThree( -antisdownThree );
  int const NineDigitSlhaCodes::antisbottomOne( antisdownThree );
  int const NineDigitSlhaCodes::sbottomOne( -antisdownThree );
  int const NineDigitSlhaCodes::antisdownFour( 200890204 );
  int const NineDigitSlhaCodes::sdownFour( -antisdownFour );
  int const NineDigitSlhaCodes::antisdownR( antisdownFour );
  int const NineDigitSlhaCodes::sdownR( -antisdownFour );
  int const NineDigitSlhaCodes::antisdownFive( 200890205 );
  int const NineDigitSlhaCodes::sdownFive( -antisdownFive );
  int const NineDigitSlhaCodes::antisstrangeR( antisdownFive );
  int const NineDigitSlhaCodes::sstrangeR( -antisdownFive );
  int const NineDigitSlhaCodes::antisdownSix( 200890206 );
  int const NineDigitSlhaCodes::sdownSix( -antisdownSix );
  int const NineDigitSlhaCodes::antisbottomTwo( antisdownSix );
  int const NineDigitSlhaCodes::sbottomTwo( -antisdownSix );
  int const NineDigitSlhaCodes::supOne( 200100401 );
  int const NineDigitSlhaCodes::antisupOne( -supOne );
  int const NineDigitSlhaCodes::supL( supOne );
  int const NineDigitSlhaCodes::antisupL( -supOne );
  int const NineDigitSlhaCodes::supTwo( 200100402 );
  int const NineDigitSlhaCodes::antisupTwo( -supTwo );
  int const NineDigitSlhaCodes::scharmL( supTwo );
  int const NineDigitSlhaCodes::antischarmL( -supTwo );
  int const NineDigitSlhaCodes::supThree( 200100403 );
  int const NineDigitSlhaCodes::antisupThree( -supThree );
  int const NineDigitSlhaCodes::stopOne( supThree );
  int const NineDigitSlhaCodes::antistopOne( -supThree );
  int const NineDigitSlhaCodes::supFour( 200100404 );
  int const NineDigitSlhaCodes::antisupFour( -supFour );
  int const NineDigitSlhaCodes::supR( supFour );
  int const NineDigitSlhaCodes::antisupR( -supFour );
  int const NineDigitSlhaCodes::supFive( 200100405 );
  int const NineDigitSlhaCodes::antisupFive( -supFive );
  int const NineDigitSlhaCodes::scharmR( supFive );
  int const NineDigitSlhaCodes::antischarmR( -supFive );
  int const NineDigitSlhaCodes::supSix( 200100406 );
  int const NineDigitSlhaCodes::antisupSix( -supSix );
  int const NineDigitSlhaCodes::stopTwo( supSix );
  int const NineDigitSlhaCodes::antistopTwo( -supSix );
  int const NineDigitSlhaCodes::neutralinoOne( 211000001 );
  int const NineDigitSlhaCodes::neutralinoTwo( 211000002 );
  int const NineDigitSlhaCodes::neutralinoThree( 211000003 );
  int const NineDigitSlhaCodes::neutralinoFour( 211000004 );
  int const NineDigitSlhaCodes::positiveCharginoOne( 210000601 );
  int const NineDigitSlhaCodes::negativeCharginoOne( -positiveCharginoOne );
  int const NineDigitSlhaCodes::positiveCharginoTwo( 210000602 );
  int const NineDigitSlhaCodes::negativeCharginoTwo( -positiveCharginoTwo );
  int const NineDigitSlhaCodes::gluinoFermion( 211110001 );

  // extra MSSM particles without parities (R, CP):
  int const NineDigitSlhaCodes::neutralColorlessScalarThree( 101000003 );
  int const
  NineDigitSlhaCodes::higgsScalarThree( neutralColorlessScalarThree );
  int const NineDigitSlhaCodes::neutralColorlessScalarFour( 101000004 );
  int const NineDigitSlhaCodes::higgsScalarFour( neutralColorlessScalarFour );
  int const NineDigitSlhaCodes::neutralColorlessScalarFive( 101000005 );
  int const NineDigitSlhaCodes::higgsScalarFive( neutralColorlessScalarFive );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarTwo( 102000002 );
  int const
  NineDigitSlhaCodes::higgsPseudoscalarTwo( neutralColorlessPseudoscalarTwo );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarThree( 102000003 );
  int const NineDigitSlhaCodes::higgsPseudoscalarThree(
                                           neutralColorlessPseudoscalarThree );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarFour( 102000004 );
  int const NineDigitSlhaCodes::higgsPseudoscalarFour(
                                            neutralColorlessPseudoscalarFour );
  int const NineDigitSlhaCodes::neutralColorlessSpinZeroBosonOne( 100000001 );
  int const
  NineDigitSlhaCodes::cpvHiggsOne( neutralColorlessSpinZeroBosonOne );
  int const NineDigitSlhaCodes::neutralColorlessSpinZeroBosonTwo( 100000002 );
  int const
  NineDigitSlhaCodes::cpvHiggsTwo( neutralColorlessSpinZeroBosonTwo );
  int const
  NineDigitSlhaCodes::neutralColorlessSpinZeroBosonThree( 100000003 );
  int const
  NineDigitSlhaCodes::cpvHiggsThree( neutralColorlessSpinZeroBosonThree );
  int const NineDigitSlhaCodes::neutralColorlessSpinZeroBosonFour( 100000004 );
  int const
  NineDigitSlhaCodes::cpvHiggsFour( neutralColorlessSpinZeroBosonFour );
  int const NineDigitSlhaCodes::neutralColorlessSpinZeroBosonFive( 100000005 );
  int const
  NineDigitSlhaCodes::cpvHiggsFive( neutralColorlessSpinZeroBosonFive );
  int const NineDigitSlhaCodes::neutralColorlessSpinZeroBosonSix( 100000006 );
  int const
  NineDigitSlhaCodes::cpvHiggsSix( neutralColorlessSpinZeroBosonSix );
  int const
  NineDigitSlhaCodes::positiveHiggsOne( positiveColorlessSpinZeroBosonOne );
  int const
  NineDigitSlhaCodes::negativeHiggsOne( -positiveColorlessSpinZeroBosonOne );
  int const NineDigitSlhaCodes::positiveColorlessSpinZeroBosonTwo( 100000602 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonTwo(
                                          -positiveColorlessSpinZeroBosonTwo );
  int const
  NineDigitSlhaCodes::positiveHiggsTwo( positiveColorlessSpinZeroBosonTwo );
  int const
  NineDigitSlhaCodes::negativeHiggsTwo( -positiveColorlessSpinZeroBosonTwo );
  int const
  NineDigitSlhaCodes::positiveColorlessSpinZeroBosonThree( 100000603 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonThree(
                                        -positiveColorlessSpinZeroBosonThree );
  int const NineDigitSlhaCodes::positiveHiggsThree(
                                         positiveColorlessSpinZeroBosonThree );
  int const NineDigitSlhaCodes::negativeHiggsThree(
                                        -positiveColorlessSpinZeroBosonThree );
  int const
  NineDigitSlhaCodes::positiveColorlessSpinZeroBosonFour( 100000604 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonFour(
                                         -positiveColorlessSpinZeroBosonFour );
  int const
  NineDigitSlhaCodes::positiveHiggsFour( positiveColorlessSpinZeroBosonFour );
  int const
  NineDigitSlhaCodes::negativeHiggsFour( -positiveColorlessSpinZeroBosonFour );
  int const
  NineDigitSlhaCodes::positiveColorlessSpinZeroBosonFive( 100000605 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonFive(
                                         -positiveColorlessSpinZeroBosonFive );
  int const
  NineDigitSlhaCodes::positiveHiggsFive( positiveColorlessSpinZeroBosonFive );
  int const
  NineDigitSlhaCodes::negativeHiggsFive( -positiveColorlessSpinZeroBosonFive );
  int const NineDigitSlhaCodes::positiveColorlessSpinZeroBosonSix( 100000606 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonSix(
                                          -positiveColorlessSpinZeroBosonSix );
  int const
  NineDigitSlhaCodes::positiveHiggsSix( positiveColorlessSpinZeroBosonSix );
  int const
  NineDigitSlhaCodes::negativeHiggsSix( -positiveColorlessSpinZeroBosonSix );
  int const
  NineDigitSlhaCodes::positiveColorlessSpinZeroBosonSeven( 100000607 );
  int const NineDigitSlhaCodes::negativeColorlessSpinZeroBosonSeven(
                                        -positiveColorlessSpinZeroBosonSeven );
  int const NineDigitSlhaCodes::positiveHiggsSeven(
                                         positiveColorlessSpinZeroBosonSeven );
  int const NineDigitSlhaCodes::negativeHiggsSeven(
                                        -positiveColorlessSpinZeroBosonSeven );
  int const NineDigitSlhaCodes::neutrinoFourMajorana( 111000004 );
  int const NineDigitSlhaCodes::neutrinoFiveMajorana( 111000005 );
  int const NineDigitSlhaCodes::neutrinoSixMajorana( 111000006 );
  int const NineDigitSlhaCodes::neutrinoSevenMajorana( 111000007 );
  int const NineDigitSlhaCodes::positronFour( 110000604 );
  int const NineDigitSlhaCodes::electronFour( -positronFour );
  int const NineDigitSlhaCodes::positiveElectronFour( positronFour );
  int const NineDigitSlhaCodes::negativeElectronFour( -positronFour );
  int const NineDigitSlhaCodes::positronFive( 110000605 );
  int const NineDigitSlhaCodes::electronFive( -positronFive );
  int const NineDigitSlhaCodes::positiveElectronFive( positronFive );
  int const NineDigitSlhaCodes::negativeElectronFive( -positronFive );

  // Next-to-Minimal Supersymmetric Standard Model (NMSSM) particles except
  // those in the MSSM above, with & without parities:
  int const NineDigitSlhaCodes::neutralColorlessScalarSix( 101000006 );
  int const NineDigitSlhaCodes::higgsScalarSix( neutralColorlessScalarSix );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarFive( 102000005 );
  int const NineDigitSlhaCodes::higgsPseudoscalarFive(
                                            neutralColorlessPseudoscalarFive );
  int const
  NineDigitSlhaCodes::neutralColorlessSpinZeroBosonSeven( 100000007 );
  int const
  NineDigitSlhaCodes::cpvHiggsSeven( neutralColorlessSpinZeroBosonSeven );
  int const
  NineDigitSlhaCodes::neutralColorlessSpinZeroBosonEight( 100000008 );
  int const
  NineDigitSlhaCodes::cpvHiggsEight( neutralColorlessSpinZeroBosonEight );
  int const NineDigitSlhaCodes::neutralinoFive( 211000004 );
  int const NineDigitSlhaCodes::neutrinoEightMajorana( 111000008 );

  // particles from the addition 3 generations of right-handed neutrino
  // superfields:
  int const NineDigitSlhaCodes::antisneutrinoFour( 200000004 );
  int const NineDigitSlhaCodes::sneutrinoFour( -antisneutrinoFour );
  int const NineDigitSlhaCodes::antisneutrinoFive( 200000005 );
  int const NineDigitSlhaCodes::sneutrinoFive( -antisneutrinoFive );
  int const NineDigitSlhaCodes::antisneutrinoSix( 200000006 );
  int const NineDigitSlhaCodes::sneutrinoSix( -antisneutrinoSix );
  int const NineDigitSlhaCodes::sneutrinoScalarFour( 201000004 );
  int const NineDigitSlhaCodes::sneutrinoScalarFive( 201000005 );
  int const NineDigitSlhaCodes::sneutrinoScalarSix( 201000006 );
  int const NineDigitSlhaCodes::sneutrinoPseudoscalarFour( 202000004 );
  int const NineDigitSlhaCodes::sneutrinoPseudoscalarFive( 202000005 );
  int const NineDigitSlhaCodes::sneutrinoPseudoscalarSix( 202000006 );
  int const NineDigitSlhaCodes::neutralColorlessScalarSeven( 101000007 );
  int const
  NineDigitSlhaCodes::higgsScalarSeven( neutralColorlessScalarSeven );
  int const NineDigitSlhaCodes::neutralColorlessScalarEight( 101000008 );
  int const
  NineDigitSlhaCodes::higgsScalarEight( neutralColorlessScalarEight );
  int const NineDigitSlhaCodes::neutralColorlessScalarNine( 101000009 );
  int const NineDigitSlhaCodes::higgsScalarNine( neutralColorlessScalarNine );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarSix( 102000006 );
  int const
  NineDigitSlhaCodes::higgsPseudoscalarSix( neutralColorlessPseudoscalarSix );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarSeven( 102000007 );
  int const NineDigitSlhaCodes::higgsPseudoscalarSeven(
                                           neutralColorlessPseudoscalarSeven );
  int const NineDigitSlhaCodes::neutralColorlessPseudoscalarEight( 102000008 );
  int const NineDigitSlhaCodes::higgsPseudoscalarEight(
                                           neutralColorlessPseudoscalarEight );
  int const NineDigitSlhaCodes::antineutrinoFour( 110000004 );
  int const NineDigitSlhaCodes::neutrinoFour( -antineutrinoFour );
  int const NineDigitSlhaCodes::antineutrinoFive( 110000005 );
  int const NineDigitSlhaCodes::neutrinoFive( -antineutrinoFive );
  int const NineDigitSlhaCodes::antineutrinoSix( 110000006 );
  int const NineDigitSlhaCodes::neutrinoSix( -antineutrinoSix );
  int const NineDigitSlhaCodes::neutrinoNineMajorana( 111000009 );
  int const NineDigitSlhaCodes::neutrinoTenMajorana( 111000010 );
  int const NineDigitSlhaCodes::neutrinoElevenMajorana( 111000011 );

}
