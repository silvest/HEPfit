/*
 * PdgData.hpp
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

#ifndef PDGDATA_HPP_
#define PDGDATA_HPP_

namespace LHPC
{
  /* this class holds data on known particles using the central values of the
   * Particle Data Group (PDG: http://pdg.lbl.gov/, all data taken from this
   * website on 2012-03-16).
   */
  class PdgData
  {
  public:
    // Standard Model particle masses:
    static double const downMass;
    static double const upMass;
    static double const strangeMass;
    static double const charmMass;
    static double const bottomMass;
    static double const topMass;
    static double const electronMass;
    static double const electronNeutrinoMass;
    static double const muonMass;
    static double const muonNeutrinoMass;
    static double const tauLeptonMass;
    static double const tauNeutrinoMass;
    static double const gluonMass;
    static double const photonMass;
    static double const zMass;
    static double const wPlusMass;

    // CKM matrix elements, required for dividing up the hadronic decays:
    static double const CkmUpDown;
    static double const CkmUpStrange;
    static double const CkmUpDownSquared;
    static double const CkmUpStrangeSquared;
    static double const CkmUpDownSquaredFraction;
    static double const CkmUpStrangeSquaredFraction;
    static double const CkmCharmDown;
    static double const CkmCharmStrange;
    static double const CkmCharmDownSquared;
    static double const CkmCharmStrangeSquared;
    static double const CkmCharmDownSquaredFraction;
    static double const CkmCharmStrangeSquaredFraction;

    // some SM particle decay widths & branching ratios:
    static double const zDecayWidth;
    static double const zToElectronAntielectronBr;
    static double const zToMuonAntimuonBr;
    static double const zToTauLeptonAntileptonBr;
    static double const zToInvisibleBr;
    static double const zToElectronNeutrinoAntineutrinoBr;
    static double const zToMuonNeutrinoAntineutrinoBr;
    static double const zToTauNeutrinoAntineutrinoBr;
    static double const zToCharmAnticharmBr;
    static double const zToBottomAntibottomBr;
    static double const zToDownAntidownBr;
    static double const zToUpAntiupBr;
    static double const zToStrangeAntistrangeBr;

    static double const wPlusDecayWidth;
    static double const wPlusToNeutrinoAntielectronBr;
    static double const wPlusToNeutrinoAntimuonBr;
    static double const wPlusToNeutrinoTauAntileptonBr;
    static double const wPlusToHadronsBr;
    static double const wPlusToCharmPlusXBr;
    static double const wPlusToCharmAntidownBr;
    static double const wPlusToCharmAntistrangeBr;
    static double const wPlusToCharmlessPlusXBr;
    static double const wPlusToUpAntidownBr;
    static double const wPlusToUpAntistrangeBr;

    static double const topDecayWidth;
    static double const topToWPlusBottomBr;

    /* currently tau leptons are treated as stable in EWScaleSpectrum, but if
     * they are to be implemented, these values should be used.
     * all these values were taken from the PDG on 2009-11-10.
     */
    static double const
    ReducedPlankConstantTimesSpeedOfLightOverTenToTheFifteenSeconds;
    static double const tauLeptonDecayWidth;
    static double const tauLeptonToNeutrinosElectronBr;
    static double const tauLeptonToNeutrinosMuonBr;
    /* this code assumes that the rest of the decay width of the tau lepton is
     * divided between down + antiup and strange + antiup in the ratio
     * ( |CkmUpDown|^2 ) to ( |CkmUpStrange|^2 ) as in the case of the decays
     * of the W^+:
     */
    static double const tauLeptonToNeutrinoHadronBr;
    static double const tauLeptonToNeutrinoDownAntiupBr;
    static double const tauLeptonToNeutrinoStrangeAntiupBr;
  };

}

#endif /* PDGDATA_HPP_ */
