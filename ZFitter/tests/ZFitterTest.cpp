/* 
 * File:   ZFitterTest.cpp
 * Author: mishima
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <StandardModel.h>
#include "ZFitter.h"
#include "ZFitterObservables.h"


/* Test function */
void test_ZFitterClass(ZFitterObservables& ZFO) {

    /* Outputs */
    double XS[12];
    double DXS[12];
    double AFB[12];
    double TAUPOL;
    double TAUAFB;
    double XSPL[12];
    double XSMI[12];
    double C1U, C1D, C2U, C2D; 

    int indexFermion;

    /* sqrt(s) = Mz */
    double sqrt_s = ZFO.getZMASS();
    double s = sqrt_s*sqrt_s;

    /* print input parameters */
    std::cout << std::endl << "##### call ZFitter::printInputs() #####"
              << std::endl << std::endl;
    ZFO.printInputs();

    /* set flags */
    std::cout << "##### call ZFitter::setAllFlags() #####" << std::endl << std::endl;
    // default: AMT4=4, ISPP=2
    int AFBC = 1, SCAL = 0, SCRE = 0, AMT4 = 6, BORN = 0,
        BOXD = 1, CONV = 1, FINR = 1, FOT2 = 3, GAMS = 1,
        DIAG = 1, INTF = 1, BARB = 2, PART = 0, POWR = 1,
        PRNT = 0, ALEM = 3, QCDC = 3, VPOL = 1, WEAK = 1,
        FTJR = 1, EXPR = 0, EXPF = 0, HIGS = 0, AFMT = 3,
        CZAK = 1, PREC = 10,HIG2 = 0, ALE2 = 3, GFER = 2,
        ISPP = 1, FSRS = 1, MISC = 0, MISD = 1, IPFC = 5,
        IPSC = 0, IPTO = 3, FBHO = 0, FSPP = 0, FUNA = 0,
        ASCR = 1, SFSR = 1, ENUE = 1, TUPV = 1, DMWW = 0,
        DSWW = 0;
    int flags[46] = {AFBC, SCAL, SCRE, AMT4, BORN, BOXD, CONV, FINR, FOT2, GAMS,
                     DIAG, INTF, BARB, PART, POWR, PRNT, ALEM, QCDC, VPOL, WEAK,
                     FTJR, EXPR, EXPF, HIGS, AFMT, CZAK, PREC, HIG2, ALE2, GFER,
                     ISPP, FSRS, MISC, MISD, IPFC, IPSC, IPTO, FBHO, FSPP, FUNA,
                     ASCR, SFSR, ENUE, TUPV, DMWW, DSWW};
    ZFO.setAllFlags(flags, 1);

    /* set cuts */
    std::cout << "##### call ZFitter::setAllCuts() #####" << std::endl << std::endl;
    int ICUT[12];
    double ACOL[12], EMIN[12], S_PR[12], ANG0[12], ANG1[12], SPP[12];
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        ICUT[indexFermion] = -1;
        ACOL[indexFermion] = 0.0;
        EMIN[indexFermion] = 0.0;
        S_PR[indexFermion] = 0.01*s;
        ANG0[indexFermion] = 0.0;
        ANG1[indexFermion] = 180.0;
        SPP[indexFermion]  = 0.25*s;
    }
    ICUT[11] = 0;
    ACOL[11] = 15.0;
    EMIN[11] = 10.0;
    S_PR[11] = 0.0;
    ANG0[11] = 35.0;
    ANG1[11] = 145.0;
    SPP[11]  = 0.0;
    ZFO.setAllCuts(ICUT, ACOL, EMIN, S_PR, ANG0, ANG1, SPP, 1);
    std::cout << std::endl; 

    /* calculate EW common blocks */
    std::cout << "##### call ZFitter::calcCommonBlocks() #####" << std::endl;
    ZFO.calcCommonBlocks();
    std::cout << std::endl;

    /* print input parameters */
    std::cout << std::endl << "##### call ZFitter::printInputs() #####"
              << std::endl << std::endl;
    ZFO.printInputs();
    std::cout << "  Note: DAL5H is calculated in calcCommonBlocks(), "
              << "if ALEM=0 or 3. " << std::endl;
    
    /* ptint constants defiend in ZFITTER */
    std::cout << std::endl << "##### call ZFitter::printConstants() #####"
              << std::endl << std::endl;
    ZFO.printConstants();

    /* print intermediate retults */
    std::cout << "##### call ZFitter::printIntermediateResults() #####"
              << std::endl << std::endl;
    ZFO.printIntermediateResults();
    
    /* Atomic Patity Violation */
    std::cout << "##### call ZFitter::calcAPV() #####" << std::endl << std::endl;
    ZFO.calcAPV(&C1U, &C1D, &C2U, &C2D);
    std::cout << "  C1U = " << C1U << " "
              << "  C1D = " << C1D << " "
              << "  C2U = " << C2U << " "
              << "  C2D = " << C2D << std::endl << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        ZFO.calcXS_AFB(indexFermion, sqrt_s,
                      &XS[indexFermion], &AFB[indexFermion]);
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion)
                  << "  " << XS[indexFermion] << "  ";
        if (indexFermion!=0&&indexFermion!=10) std::cout << AFB[indexFermion];
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total and partial decay widths */
        double Gamma_Z = ZFO.getCommonWIDTHS(11);
        double Gamma_e = ZFO.getCommonWIDTHS(1);
        double Gamma_f = Gamma_e;
        if (indexFermion!=11) Gamma_f = ZFO.getCommonWIDTHS(indexFermion);
        ZFO.calcXS(indexFermion, sqrt_s, Gamma_Z, Gamma_e, Gamma_f,
                  &XS[indexFermion]);
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion)
                  << "  " << XS[indexFermion]
                  << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB_2() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZFO.getCommonWIDTHS(11);
        /* get the effective vector and axial-vector couplings */
        double GAE = sqrt(ZFO.getCommonAROTFZ(1))*0.5;
        double GVE = ZFO.getCommonARVEFZ(1)*GAE;
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZFO.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZFO.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = GAE;
            GVF = ZFO.getCommonARVEFZ(1)*GAF;
        }
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion) << "  ";
        if (indexFermion!=0&&indexFermion!=10) {
            ZFO.calcXS_AFB_2(indexFermion, sqrt_s, Gamma_Z, 0,
                            GVE, GAE, GVF, GAF,
                            &XS[indexFermion], &AFB[indexFermion]);
            std::cout << XS[indexFermion] << "  " << AFB[indexFermion];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB_3() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZFO.getCommonWIDTHS(11);
        /* get the effective vector and axial-vector couplings */
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZFO.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZFO.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = sqrt(ZFO.getCommonAROTFZ(1))*0.5;
            GVF = ZFO.getCommonARVEFZ(1)*GAF;
        }
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion) << "  ";
        double GVF2 = GVF*GVF;
        double GAF2 = GAF*GAF;
        if (indexFermion==1||indexFermion==2||indexFermion==3||indexFermion==11) {
            ZFO.calcXS_AFB_3(indexFermion, sqrt_s, Gamma_Z, 0, GVF2, GAF2,
                            &XS[indexFermion], &AFB[indexFermion]);
            std::cout << XS[indexFermion] << "  " << AFB[indexFermion];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

   std::cout << "##### call ZFitter::calcXS_AFB_4() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZFO.getCommonWIDTHS(11);
        /* get the effective vector and axial-vector couplings */
        double GAE = sqrt(ZFO.getCommonAROTFZ(1))*0.5;
        double GVE = ZFO.getCommonARVEFZ(1)*GAE;
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZFO.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZFO.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = GAE;
            GVF = ZFO.getCommonARVEFZ(1)*GAF;
        }
        double PFOUR = GVE*GAE*GVF*GAF;
        double PVAE2 = GVE*GVE + GAE*GAE;
        double PVAF2 = GVF*GVF + GAF*GAF;
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion) << "  ";
        if (indexFermion==1||indexFermion==2||indexFermion==3) {
            ZFO.calcXS_AFB_4(indexFermion, sqrt_s, Gamma_Z, PFOUR, PVAE2, PVAF2,
                            &XS[indexFermion], &AFB[indexFermion]);
            std::cout << XS[indexFermion] << "  " << AFB[indexFermion];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcTauPol() #####" << std::endl
              << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Tau polarization, Tau polarization asymmetry"
              << std::endl;
    ZFO.calcTauPol(sqrt_s, &TAUPOL, &TAUAFB);
    std::cout << std::setw(9) << ZFO.convertINDF(3) << "  "
              << TAUPOL << "  " << TAUAFB << std::endl
              << std::endl;

    std::cout << "##### call ZFitter::calcTauPol_2() #####" << std::endl
              << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Tau polarization, Tau polarization asymmetry"
              << std::endl;
    /* get the total width */
    double Gamma_Z = ZFO.getCommonWIDTHS(11);
    /* get the effective vector and axial-vector couplings */
    double GAE = sqrt(ZFO.getCommonAROTFZ(1))*0.5;
    double GVE = ZFO.getCommonARVEFZ(1) * GAE;
    double GAF = sqrt(ZFO.getCommonAROTFZ(3))*0.5;
    double GVF = ZFO.getCommonARVEFZ(3) * GAF;
    ZFO.calcTauPol_2(sqrt_s, Gamma_Z, 0, GVE, GAE, GVF, GAF, &TAUPOL, &TAUAFB);
    std::cout << std::setw(9) << ZFO.convertINDF(3) << "  "
              << TAUPOL << "  " << TAUAFB << std::endl
              << std::endl;

    std::cout << "##### call ZFitter::calcALR() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Left-right asymmetries " << std::endl;
    double POL = 0.63;
    double ALRI;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        ZFO.calcALR(indexFermion, sqrt_s, POL,
                   &XSPL[indexFermion], &XSMI[indexFermion]);
        ALRI = (XSMI[indexFermion] - XSPL[indexFermion])
                / (XSMI[indexFermion] + XSPL[indexFermion]) / POL;
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion) << "  ";
        if (indexFermion!=8) std::cout << ALRI;
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcDXS() "
              << "with the flag ISPP=2 #####" << std::endl << std::endl;
    ZFO.flag("ISPP", 2); // default: 2
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Differential cross sections" << std::endl;
    for (indexFermion = 0; indexFermion < 11; indexFermion++) {
        std::cout << std::setw(9) << ZFO.convertINDF(indexFermion) << "  ";
        for (int i=0; i<=6; i++) {
            double CSA = -1.0 + 2.0/6.0*(double)i; // cos(theta)
            ZFO.calcDXS(indexFermion, sqrt_s, CSA, &DXS[indexFermion]);
            std::cout << std::setw(9) << DXS[indexFermion] << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::printPO() #####" << std::endl << std::endl;
    ZFO.printPO();

    
    /* Outputs from ZFITTER subroutines (for tests) */
    //
    //std::cout << "##### call ZFitter::info(0) #####"
    //          << std::endl << std::endl;
    //std::cout << " ------------- Flag info ---------------" << std::endl;
    //ZFO.FlagInfo(); // output flag info
    //
    //std::cout << std::endl << "##### call ZFitter::info(1) #####"
    //          << std::endl << std::endl;
    //std::cout << " ------------------------ Cut info --------------------------";
    //ZFO.CutInfo(); // print cut info

}


/* Test with ZFITTER subroutine ZFTEST() */
void test_ZFTEST(ZFitterObservables& ZFO) {
    ZFO.test(0);
    //ZFO.test(1);
}


void setSMparameters(StandardModel& SM_i) {
    std::map<std::string, double> Parameters;
    // 17 parameters defined in StandardModel
    Parameters["GF"] = 1.16637E-5;
    Parameters["mneutrino_1"] = 0.0;
    Parameters["mneutrino_2"] = 0.0;
    Parameters["mneutrino_3"] = 0.0;
    Parameters["melectron"] = 0.54857990943e-3;
    Parameters["mmu"] = 0.1134289256;
    Parameters["mtau"] = 1.77682;
    Parameters["lambda"] = 0.2253;
    Parameters["A"] = 0.808;
    Parameters["rhob"] = 0.132;
    Parameters["etab"] = 0.341;
    Parameters["ale"] = 1.0/137.035999679;
    Parameters["dAle5Mz"] = 0.02758; // not used when ALEM=3
    Parameters["mHl"] = 130.0;
    Parameters["muw"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
    // 26 parameters defined in QCD    
    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    Parameters["mtop"] = 174.0;
    Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
    Parameters["MBd"] = 0.0;
    Parameters["MBs"] = 0.0;
    Parameters["MBp"] = 0.0;
    Parameters["MK0"] = 0.0;
    Parameters["MKp"] = 0.0;
    Parameters["FBs"] = 0.0;
    Parameters["FBsoFBd"] = 0.0;
    Parameters["BBsoBBd"] = 0.0;
    Parameters["BBs1"] = 0.0;
    Parameters["BBs2"] = 0.0;
    Parameters["BBs3"] = 0.0;
    Parameters["BBs4"] = 0.0;
    Parameters["BBs5"] = 0.0;
    Parameters["BBsscale"] = 0.0;
    Parameters["BBsscheme"] = 0.0;

    /* TEST for ZFITTER */
    Parameters["Mz"] = 91.1876;
    Parameters["mtop"] = 178.0;    
    Parameters["mHl"] = 100.0;
    Parameters["AlsMz"] = 0.117;
    Parameters["dAle5Mz"] = 0.027572;
    
    SM_i.Init(Parameters);
}



int main(int argc, char** argv) {

    StandardModel* myModel;
    myModel = new StandardModel();
    setSMparameters(*myModel);
    ZFitterObservables ZFO(*myModel);

    if (argc>1) {
        if (strcmp(argv[1], "-t")==0) {
            /* test with Subroutine ZFTEST */
            test_ZFTEST(ZFO);
        } else {
            std::cout << "use the option -t for ZFTEST()" << std::endl;
        }
    } else {
        /* compute EW precision observables */
        test_ZFitterClass(ZFO);
    }
    
    return (EXIT_SUCCESS);
}

