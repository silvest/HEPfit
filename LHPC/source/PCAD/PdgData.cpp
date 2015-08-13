/*
 * PdgData.cpp
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
  // Standard Model particle masses:
  double const PdgData::downMass( 0.0050 );
  double const PdgData::upMass( 0.0025 );
  double const PdgData::strangeMass( 0.100 );
  double const PdgData::charmMass( 1.29 );
  double const PdgData::bottomMass( 4.19 );
  double const PdgData::topMass( 172.9 );
  double const PdgData::electronMass( 0.000510998910 );
  double const PdgData::electronNeutrinoMass( 0.0 );
  double const PdgData::muonMass( 0.1056583668 );
  double const PdgData::muonNeutrinoMass( 0.0 );
  double const PdgData::tauLeptonMass( 1.77682 );
  double const PdgData::tauNeutrinoMass( 0.0 );
  double const PdgData::gluonMass( 0.0 );
  double const PdgData::photonMass( 0.0 );
  double const PdgData::zMass( 91.1876 );
  double const PdgData::wPlusMass( 80.399 );

  double const PdgData::CkmUpDown( 0.97418 );
  double const PdgData::CkmUpStrange( 0.2255 );
  double const
  PdgData::CkmUpDownSquared( PdgData::CkmUpDown * PdgData::CkmUpDown );
  double const PdgData::CkmUpStrangeSquared( PdgData::CkmUpStrange
                                             * PdgData::CkmUpStrange );
  double const
  PdgData::CkmUpDownSquaredFraction( PdgData::CkmUpDownSquared
                                     / ( PdgData::CkmUpDownSquared
                                         + PdgData::CkmUpStrangeSquared ) );
  double const
  PdgData::CkmUpStrangeSquaredFraction( PdgData::CkmUpStrangeSquared
                                        / ( PdgData::CkmUpDownSquared
                                            + PdgData::CkmUpStrangeSquared ) );
  double const PdgData::CkmCharmDown( 0.23 );
  double const PdgData::CkmCharmStrange( 1.04 );
  double const PdgData::CkmCharmDownSquared( PdgData::CkmCharmDown
                                             * PdgData::CkmCharmDown );
  double const PdgData::CkmCharmStrangeSquared( PdgData::CkmCharmStrange
                                                * PdgData::CkmCharmStrange );
  double const
  PdgData::CkmCharmDownSquaredFraction( PdgData::CkmCharmDownSquared
                                        / ( PdgData::CkmCharmDownSquared
                                         + PdgData::CkmCharmStrangeSquared ) );
  double const
  PdgData::CkmCharmStrangeSquaredFraction( PdgData::CkmCharmStrangeSquared
                                           / ( PdgData::CkmCharmDownSquared
                                         + PdgData::CkmCharmStrangeSquared ) );

  double const PdgData::zDecayWidth( 2.4952 );
  double const PdgData::zToElectronAntielectronBr( 0.03363 );
  double const PdgData::zToMuonAntimuonBr( 0.03366 );
  double const PdgData::zToTauLeptonAntileptonBr( 0.03367 );
  double const PdgData::zToInvisibleBr( 0.20 );

  // this code assumes that the invisible decays of Z bosons are equally
  // divided between the 3 flavors of neutrino:
  double const
  PdgData::zToElectronNeutrinoAntineutrinoBr( PdgData::zToInvisibleBr / 3.0 );
  double const PdgData::zToMuonNeutrinoAntineutrinoBr(
                                  PdgData::zToElectronNeutrinoAntineutrinoBr );
  double const PdgData::zToTauNeutrinoAntineutrinoBr(
                                  PdgData::zToElectronNeutrinoAntineutrinoBr );

  // this code assumes that the rest of the hadronic decay width of the Z boson
  // is equally divided between down, up and strange:
  double const PdgData::zToCharmAnticharmBr( 0.1203 );
  double const PdgData::zToBottomAntibottomBr( 0.1512 );
  double const
  PdgData::zToDownAntidownBr ( ( 0.6991 - PdgData::zToCharmAnticharmBr
                                 - PdgData::zToBottomAntibottomBr ) / 3.0 );
  double const
  PdgData::zToUpAntiupBr( PdgData::zToDownAntidownBr );
  double const
  PdgData::zToStrangeAntistrangeBr( PdgData::zToDownAntidownBr );

  double const PdgData::wPlusDecayWidth( 2.085 );
  double const PdgData::wPlusToNeutrinoAntielectronBr( 0.1075 );
  double const PdgData::wPlusToNeutrinoAntimuonBr( 0.1057 );
  double const PdgData::wPlusToNeutrinoTauAntileptonBr( 0.1125 );
  double const PdgData::wPlusToHadronsBr( 0.6760 );

  /* this code assumes that the BRs of the W^+ into charm + antidown &
   * charm + antistrange account for all the BR of the charm + X,
   * in the ratio ( |CkmCharmDown|^2 ) to ( |CkmCharmStrange|^2 ):
   */
  double const PdgData::wPlusToCharmPlusXBr( 0.334 );
  double const
  PdgData::wPlusToCharmAntidownBr( PdgData::CkmCharmDownSquaredFraction
                                   * PdgData::wPlusToCharmPlusXBr );
  double const
  PdgData::wPlusToCharmAntistrangeBr( PdgData::CkmCharmStrangeSquaredFraction
                                      * PdgData::wPlusToCharmPlusXBr );

  /* this code assumes that the rest of the hadronic decay width of the W is
   * divided between up + antidown and up + antistrange
   * in the ratio ( |CkmUpDown|^2 ) to ( |CkmUpStrange|^2 ):
   */
  double const PdgData::wPlusToCharmlessPlusXBr( PdgData::wPlusToHadronsBr
                                              - PdgData::wPlusToCharmPlusXBr );
  double const PdgData::wPlusToUpAntidownBr( PdgData::CkmUpDownSquaredFraction
                                          * PdgData::wPlusToCharmlessPlusXBr );
  double const
  PdgData::wPlusToUpAntistrangeBr( PdgData::CkmUpStrangeSquaredFraction
                                   * PdgData::wPlusToCharmlessPlusXBr );

  double const PdgData::topDecayWidth( 1.99 );
  double const PdgData::topToWPlusBottomBr( 1.0 );

  /* currently tau leptons are treated as stable by the constructor, but if
   * they are to be implemented, these values should be used.
   * all these values were taken from the PDG on 2009-11-10.
   */
  double const
  PdgData::ReducedPlankConstantTimesSpeedOfLightOverTenToTheFifteenSeconds(
                                                        0.000000000658211899 );
  double const PdgData::tauLeptonDecayWidth(
     ReducedPlankConstantTimesSpeedOfLightOverTenToTheFifteenSeconds / 290.6 );
  double const PdgData::tauLeptonToNeutrinosElectronBr( 0.1782 );
  double const PdgData::tauLeptonToNeutrinosMuonBr( 0.1739 );
  /* this code assumes that the rest of the decay width of the tau lepton is
   * divided between down + antiup and strange + antiup in the ratio
   * ( |CkmUpDown|^2 ) to ( |CkmUpStrange|^2 ) as in the case of the decays of
   * the W^+:
   */
  double const PdgData::tauLeptonToNeutrinoHadronBr( 1.0
                                      - PdgData::tauLeptonToNeutrinosElectronBr
                                       - PdgData::tauLeptonToNeutrinosMuonBr );
  double const
  PdgData::tauLeptonToNeutrinoDownAntiupBr( PdgData::CkmUpDownSquaredFraction
                                      * PdgData::tauLeptonToNeutrinoHadronBr );
  double const PdgData::tauLeptonToNeutrinoStrangeAntiupBr(
                                           PdgData::CkmUpStrangeSquaredFraction
                                      * PdgData::tauLeptonToNeutrinoHadronBr );
}
