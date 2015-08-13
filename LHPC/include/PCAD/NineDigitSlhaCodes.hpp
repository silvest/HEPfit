/*
 * NineDigitSlhaCodes.hpp
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

#ifndef NINEDIGITSLHACODES_HPP_
#define NINEDIGITSLHACODES_HPP_

namespace LHPC
{
  // this is a class to basically just hold a lot of synonyms for particle
  // codes in the proposed 9-digits-in-base-10 format.
  class NineDigitSlhaCodes
  {
  public:
    // Standard Model (SM) particles:
    static int const neutralColorlessScalarOne;
    static int const higgsBoson;
    static int const positronOne;
    static int const electronOne;
    static int const positiveElectron;
    static int const negativeElectron;
    static int const positronTwo;
    static int const electronTwo;
    static int const positiveMuon;
    static int const negativeMuon;
    static int const positronThree;
    static int const electronThree;
    static int const positiveTau;
    static int const negativeTau;
    static int const antineutrinoOne;
    static int const neutrinoOne;
    static int const electronAntineutrino;
    static int const electronNeutrino;
    static int const antineutrinoTwo;
    static int const neutrinoTwo;
    static int const muonAntineutrino;
    static int const muonNeutrino;
    static int const antineutrinoThree;
    static int const neutrinoThree;
    static int const tauAntineutrino;
    static int const tauNeutrino;
    static int const neutrinoOneMajorana;
    static int const neutrinoTwoMajorana;
    static int const neutrinoThreeMajorana;
    static int const antidownOne;
    static int const downOne;
    static int const downAntiquark;
    static int const downQuark;
    static int const antidownTwo;
    static int const downTwo;
    static int const strangeAntiquark;
    static int const strangeQuark;
    static int const antidownThree;
    static int const downThree;
    static int const bottomAntiquark;
    static int const bottomQuark;
    static int const upOne;
    static int const antiupOne;
    static int const upQuark;
    static int const upAntiquark;
    static int const upTwo;
    static int const antiupTwo;
    static int const charmQuark;
    static int const charmAntiquark;
    static int const upThree;
    static int const antiupThree;
    static int const topQuark;
    static int const topAntiquark;
    static int const photonBoson;
    static int const zBosonOne;
    static int const zBoson;
    static int const wPlusBosonOne;
    static int const wPlus;
    static int const wMinus;
    static int const gluonBoson;

    // almost-SM particles:
    static int const gravitonBoson;

    // Minimal Supersymmetric Standard Model (MSSM) particles except those in
    // the SM above, assuming R-parity conservation:
    static int const higgsScalarOne;
    static int const lightHiggs;
    static int const neutralColorlessScalarTwo;
    static int const higgsScalarTwo;
    static int const heavyHiggs;
    static int const neutralColorlessPseudoscalarOne;
    static int const higgsPseudoscalar;
    static int const positiveColorlessSpinZeroBosonOne;
    static int const negativeColorlessSpinZeroBosonOne;
    static int const positiveHiggs;
    static int const negativeHiggs;
    static int const spositronOne;
    static int const selectronOne;
    static int const antiselectronL;
    static int const selectronL;
    static int const positiveSelectronL;
    static int const negativeSelectronL;
    static int const spositronTwo;
    static int const selectronTwo;
    static int const antismuonL;
    static int const smuonL;
    static int const positiveSmuonL;
    static int const negativeSmuonL;
    static int const spositronThree;
    static int const selectronThree;
    static int const antistauOne;
    static int const stauOne;
    static int const positiveStauOne;
    static int const negativeStauOne;
    static int const spositronFour;
    static int const selectronFour;
    static int const antiselectronR;
    static int const selectronR;
    static int const positiveSelectronR;
    static int const negativeSelectronR;
    static int const spositronFive;
    static int const selectronFive;
    static int const antismuonR;
    static int const smuonR;
    static int const positiveSmuonR;
    static int const negativeSmuonR;
    static int const spositronSix;
    static int const selectronSix;
    static int const antistauTwo;
    static int const stauTwo;
    static int const positiveStauTwo;
    static int const negativeStauTwo;
    static int const antisneutrinoOne;
    static int const sneutrinoOne;
    static int const electronAntisneutrinoL;
    static int const electronSneutrinoL;
    static int const antisneutrinoTwo;
    static int const sneutrinoTwo;
    static int const muonAntisneutrinoL;
    static int const muonSneutrinoL;
    static int const antisneutrinoThree;
    static int const sneutrinoThree;
    static int const tauAntisneutrinoL;
    static int const tauSneutrinoL;
    static int const sneutrinoScalarOne;
    static int const sneutrinoScalarTwo;
    static int const sneutrinoScalarThree;
    static int const sneutrinoPseudoscalarOne;
    static int const sneutrinoPseudoscalarTwo;
    static int const sneutrinoPseudoscalarThree;
    static int const antisdownOne;
    static int const sdownOne;
    static int const antisdownL;
    static int const sdownL;
    static int const antisdownTwo;
    static int const sdownTwo;
    static int const antisstrangeL;
    static int const sstrangeL;
    static int const antisdownThree;
    static int const sdownThree;
    static int const antisbottomOne;
    static int const sbottomOne;
    static int const antisdownFour;
    static int const sdownFour;
    static int const antisdownR;
    static int const sdownR;
    static int const antisdownFive;
    static int const sdownFive;
    static int const antisstrangeR;
    static int const sstrangeR;
    static int const antisdownSix;
    static int const sdownSix;
    static int const antisbottomTwo;
    static int const sbottomTwo;
    static int const supOne;
    static int const antisupOne;
    static int const supL;
    static int const antisupL;
    static int const supTwo;
    static int const antisupTwo;
    static int const scharmL;
    static int const antischarmL;
    static int const supThree;
    static int const antisupThree;
    static int const stopOne;
    static int const antistopOne;
    static int const supFour;
    static int const antisupFour;
    static int const supR;
    static int const antisupR;
    static int const supFive;
    static int const antisupFive;
    static int const scharmR;
    static int const antischarmR;
    static int const supSix;
    static int const antisupSix;
    static int const stopTwo;
    static int const antistopTwo;
    static int const neutralinoOne;
    static int const neutralinoTwo;
    static int const neutralinoThree;
    static int const neutralinoFour;
    static int const positiveCharginoOne;
    static int const negativeCharginoOne;
    static int const positiveCharginoTwo;
    static int const negativeCharginoTwo;
    static int const gluinoFermion;

    // extra MSSM particles without parities (R, CP):
    static int const neutralColorlessScalarThree;
    static int const higgsScalarThree;
    static int const neutralColorlessScalarFour;
    static int const higgsScalarFour;
    static int const neutralColorlessScalarFive;
    static int const higgsScalarFive;
    static int const neutralColorlessPseudoscalarTwo;
    static int const higgsPseudoscalarTwo;
    static int const neutralColorlessPseudoscalarThree;
    static int const higgsPseudoscalarThree;
    static int const neutralColorlessPseudoscalarFour;
    static int const higgsPseudoscalarFour;
    static int const neutralColorlessSpinZeroBosonOne;
    static int const cpvHiggsOne;
    static int const neutralColorlessSpinZeroBosonTwo;
    static int const cpvHiggsTwo;
    static int const neutralColorlessSpinZeroBosonThree;
    static int const cpvHiggsThree;
    static int const neutralColorlessSpinZeroBosonFour;
    static int const cpvHiggsFour;
    static int const neutralColorlessSpinZeroBosonFive;
    static int const cpvHiggsFive;
    static int const neutralColorlessSpinZeroBosonSix;
    static int const cpvHiggsSix;
    static int const positiveHiggsOne;
    static int const negativeHiggsOne;
    static int const positiveColorlessSpinZeroBosonTwo;
    static int const negativeColorlessSpinZeroBosonTwo;
    static int const positiveHiggsTwo;
    static int const negativeHiggsTwo;
    static int const positiveColorlessSpinZeroBosonThree;
    static int const negativeColorlessSpinZeroBosonThree;
    static int const positiveHiggsThree;
    static int const negativeHiggsThree;
    static int const positiveColorlessSpinZeroBosonFour;
    static int const negativeColorlessSpinZeroBosonFour;
    static int const positiveHiggsFour;
    static int const negativeHiggsFour;
    static int const positiveColorlessSpinZeroBosonFive;
    static int const negativeColorlessSpinZeroBosonFive;
    static int const positiveHiggsFive;
    static int const negativeHiggsFive;
    static int const positiveColorlessSpinZeroBosonSix;
    static int const negativeColorlessSpinZeroBosonSix;
    static int const positiveHiggsSix;
    static int const negativeHiggsSix;
    static int const positiveColorlessSpinZeroBosonSeven;
    static int const negativeColorlessSpinZeroBosonSeven;
    static int const positiveHiggsSeven;
    static int const negativeHiggsSeven;
    static int const neutrinoFourMajorana;
    static int const neutrinoFiveMajorana;
    static int const neutrinoSixMajorana;
    static int const neutrinoSevenMajorana;
    static int const positronFour;
    static int const electronFour;
    static int const positiveElectronFour;
    static int const negativeElectronFour;
    static int const positronFive;
    static int const electronFive;
    static int const positiveElectronFive;
    static int const negativeElectronFive;

    // Next-to-Minimal Supersymmetric Standard Model (NMSSM) particles except
    // those in the MSSM above, with & without parities:
    static int const neutralColorlessScalarSix;
    static int const higgsScalarSix;
    static int const neutralColorlessPseudoscalarFive;
    static int const higgsPseudoscalarFive;
    static int const neutralColorlessSpinZeroBosonSeven;
    static int const cpvHiggsSeven;
    static int const neutralColorlessSpinZeroBosonEight;
    static int const cpvHiggsEight;
    static int const neutralinoFive;
    static int const neutrinoEightMajorana;

    // particles from the addition 3 generations of right-handed neutrino
    // superfields:
    static int const antisneutrinoFour;
    static int const sneutrinoFour;
    static int const antisneutrinoFive;
    static int const sneutrinoFive;
    static int const antisneutrinoSix;
    static int const sneutrinoSix;
    static int const sneutrinoScalarFour;
    static int const sneutrinoScalarFive;
    static int const sneutrinoScalarSix;
    static int const sneutrinoPseudoscalarFour;
    static int const sneutrinoPseudoscalarFive;
    static int const sneutrinoPseudoscalarSix;
    static int const neutralColorlessScalarSeven;
    static int const higgsScalarSeven;
    static int const neutralColorlessScalarEight;
    static int const higgsScalarEight;
    static int const neutralColorlessScalarNine;
    static int const higgsScalarNine;
    static int const neutralColorlessPseudoscalarSix;
    static int const higgsPseudoscalarSix;
    static int const neutralColorlessPseudoscalarSeven;
    static int const higgsPseudoscalarSeven;
    static int const neutralColorlessPseudoscalarEight;
    static int const higgsPseudoscalarEight;
    static int const antineutrinoFour;
    static int const neutrinoFour;
    static int const antineutrinoFive;
    static int const neutrinoFive;
    static int const antineutrinoSix;
    static int const neutrinoSix;
    static int const neutrinoNineMajorana;
    static int const neutrinoTenMajorana;
    static int const neutrinoElevenMajorana;
  };
  typedef NineDigitSlhaCodes PDGIX;

}

#endif /* NINEDIGITSLHACODES_HPP_ */
