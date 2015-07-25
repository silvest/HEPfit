/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <StandardModel.h>
#include "ZFitterWrapper.h"

using namespace std;


void printPO(ZFitterWrapper& ZF_i) {
    cout << setw(15) << "m_W [GeV]" << setw(13) << ZF_i.Mw() << endl
         << setw(15) << "Gamma_W [GeV]" << setw(13)
         << ZF_i.Gamma_W() << endl
         << setw(15) << "sin^2(th_W)" << setw(13)
         << ZF_i.sw2() << endl
         << setw(15) << "sin^2(teff_e)" << setw(13)
         << ZF_i.s2teff_f(1) << endl
         << setw(15) << "sin^2(teff_mu)" << setw(13)
         << ZF_i.s2teff_f(2) << endl
         << setw(15) << "sin^2(teff_tau)" << setw(13)
         << ZF_i.s2teff_f(3) << endl
         << setw(15) << "sin^2(teff_b)" << setw(13)
         << ZF_i.s2teff_f(9) << endl
         << setw(15) << "sin^2(teff_c)" << setw(13)
         << ZF_i.s2teff_f(6) << endl
         << setw(15) << "sin^2(teff_s)" << setw(13)
         << ZF_i.s2teff_f(7) << endl
         << setw(15) << "Gamma_inv [GeV]" << setw(13)
         << ZF_i.Gamma_inv() << endl
         << setw(15) << "Gamma_had [GeV]" << setw(13)
         << ZF_i.Gamma_had() << endl
         << setw(15) << "Gamma_Z [GeV]" << setw(13)
         << ZF_i.Gamma_Z() << endl
         << endl;
}


/* Test function */
void test_ZFitterClass(ZFitterWrapper& ZF) {

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
    double sqrt_s = ZF.getModel().getMz();
    
    /* print input parameters */
    cout << endl << "##### call ZFitter::printInputs() #####" << endl << endl;
    ZF.printInputs();
    //cout << "  Note: DAL5H is calculated in calcCommonBlocks(), "
    //     << "if ALEM=0 or 3. " << endl;

    cout << "##### Pseudo-observables #####" << endl << endl;
    printPO(ZF);
    cout << endl;

    /* print constants defined in ZFITTER */
    cout << endl << "##### call ZFitter::printConstants() #####" << endl << endl;
    ZF.printConstants();
    
    /* print intermediate results */
    cout << "##### call ZFitter::printIntermediateResults() #####"
              << endl << endl;
    ZF.printIntermediateResults();

    /* Atomic Parity Violation */
    cout << "##### call ZFitter::calcAPV() #####" << endl << endl;
    ZF.calcAPV(&C1U, &C1D, &C2U, &C2D);
    cout << "  C1U = " << C1U << " "
         << "  C1D = " << C1D << " "
         << "  C2U = " << C2U << " "
         << "  C2D = " << C2D << endl << endl;

    cout << "##### call ZFitter::calcXS_AFB() #####" << endl << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Cross sections, FB asymmetries" << endl;
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        ZF.calcXS_AFB(indexFermion, sqrt_s,
                      &XS[indexFermion], &AFB[indexFermion]);
        cout << setw(9) << ZF.convertINDF(indexFermion)
             << "  " << XS[indexFermion] << "  ";
        if (indexFermion!=0&&indexFermion!=10) cout << AFB[indexFermion];
        cout << endl;
    }
    cout << endl;

    cout << "##### call ZFitter::calcXS() #####" << endl << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Cross sections" << endl;
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        /* get the total and partial decay widths */
        double Gamma_Z = ZF.Gamma_Z();
        double Gamma_e = ZF.Gamma_f(1);
        double Gamma_f = Gamma_e;
        if (indexFermion!=11) Gamma_f = ZF.Gamma_f(indexFermion);
        ZF.calcXS(indexFermion, sqrt_s, Gamma_Z, Gamma_e, Gamma_f,
                  &XS[indexFermion]);
        cout << setw(9) << ZF.convertINDF(indexFermion)
             << "  " << XS[indexFermion]
             << endl;
    }
    cout << endl;

    cout << "##### call ZFitter::calcXS_AFB_2() #####" << endl << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Cross sections, FB asymmetries" << endl;
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.Gamma_Z();
        /* get the effective vector and axial-vector couplings */
        double GAE = sqrt(ZF.rhoZ_f(1).real())*0.5;
        double GVE = ZF.gZ_f(1).real()*GAE;
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.rhoZ_f(indexFermion).real())*0.5;
            GVF = ZF.gZ_f(indexFermion).real()*GAF;
        } else {
            GAF = GAE;
            GVF = ZF.gZ_f(1).real()*GAF;
        }
        cout << setw(9) << ZF.convertINDF(indexFermion) << "  ";
        if (indexFermion!=0&&indexFermion!=10) {
            ZF.calcXS_AFB_2(indexFermion, sqrt_s, Gamma_Z, 0,
                            GVE, GAE, GVF, GAF,
                            &XS[indexFermion], &AFB[indexFermion]);
            cout << XS[indexFermion] << "  " << AFB[indexFermion];
        }
        cout << endl;
    }
    cout << endl;

    cout << "##### call ZFitter::calcXS_AFB_3() #####" << endl << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Cross sections, FB asymmetries" << endl;
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.Gamma_Z();
        /* get the effective vector and axial-vector couplings */
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.rhoZ_f(indexFermion).real())*0.5;
            GVF = ZF.gZ_f(indexFermion).real()*GAF;
        } else {
            GAF = sqrt(ZF.rhoZ_f(1).real())*0.5;
            GVF = ZF.gZ_f(1).real()*GAF;
        }
        cout << setw(9) << ZF.convertINDF(indexFermion) << "  ";
        double GVF2 = GVF*GVF;
        double GAF2 = GAF*GAF;
        if (indexFermion==1||indexFermion==2||indexFermion==3||indexFermion==11) {
            ZF.calcXS_AFB_3(indexFermion, sqrt_s, Gamma_Z, 0, GVF2, GAF2,
                            &XS[indexFermion], &AFB[indexFermion]);
            cout << XS[indexFermion] << "  " << AFB[indexFermion];
        }
        cout << endl;
    }
    cout << endl;

   cout << "##### call ZFitter::calcXS_AFB_4() #####" << endl << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Cross sections, FB asymmetries" << endl;
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.Gamma_Z();
        /* get the effective vector and axial-vector couplings */
        double GAE = sqrt(ZF.rhoZ_f(1).real())*0.5;
        double GVE = ZF.gZ_f(1).real()*GAE;
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.rhoZ_f(indexFermion).real())*0.5;
            GVF = ZF.gZ_f(indexFermion).real()*GAF;
        } else {
            GAF = GAE;
            GVF = ZF.gZ_f(1).real()*GAF;
        }
        double PFOUR = GVE*GAE*GVF*GAF;
        double PVAE2 = GVE*GVE + GAE*GAE;
        double PVAF2 = GVF*GVF + GAF*GAF;
        cout << setw(9) << ZF.convertINDF(indexFermion) << "  ";
        if (indexFermion==1||indexFermion==2||indexFermion==3) {
            ZF.calcXS_AFB_4(indexFermion, sqrt_s, Gamma_Z, PFOUR, PVAE2, PVAF2,
                            &XS[indexFermion], &AFB[indexFermion]);
            cout << XS[indexFermion] << "  " << AFB[indexFermion];
        }
        cout << endl;
    }
    cout << endl;

    cout << "##### call ZFitter::calcTauPol() #####" << endl
              << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Tau polarization, Tau polarization asymmetry"
              << endl;
    ZF.calcTauPol(sqrt_s, &TAUPOL, &TAUAFB);
    cout << setw(9) << ZF.convertINDF(3) << "  "
              << TAUPOL << "  " << TAUAFB << endl
              << endl;

    cout << "##### call ZFitter::calcTauPol_2() #####" << endl
              << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Tau polarization, Tau polarization asymmetry"
              << endl;
    /* get the total width */
    double Gamma_Z = ZF.Gamma_Z();
    /* get the effective vector and axial-vector couplings */
    double GAE = sqrt(ZF.rhoZ_f(1).real())*0.5;
    double GVE = ZF.gZ_f(1).real() * GAE;
    double GAF = sqrt(ZF.rhoZ_f(3).real())*0.5;
    double GVF = ZF.gZ_f(3).real() * GAF;
    ZF.calcTauPol_2(sqrt_s, Gamma_Z, 0, GVE, GAE, GVF, GAF, &TAUPOL, &TAUAFB);
    cout << setw(9) << ZF.convertINDF(3) << "  "
              << TAUPOL << "  " << TAUAFB << endl
              << endl;

    cout << "##### call ZFitter::calcALR() #####" << endl << endl;
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Left-right asymmetries " << endl;
    double POL = 0.63;
    double ALRI;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        ZF.calcALR(indexFermion, sqrt_s, POL,
                   &XSPL[indexFermion], &XSMI[indexFermion]);
        ALRI = (XSMI[indexFermion] - XSPL[indexFermion])
                / (XSMI[indexFermion] + XSPL[indexFermion]) / POL;
        cout << setw(9) << ZF.convertINDF(indexFermion) << "  ";
        if (indexFermion!=8) cout << ALRI;
        cout << endl;
    }
    cout << endl;

    cout << "##### call ZFitter::calcDXS() "
              << "with the flag ISPP=2 #####" << endl << endl;
    ZF.setFlag("ISPP", 2); // default: 2
    cout << "  sqrt(s) = " << sqrt_s << endl << endl;
    cout << "  Channel, Differential cross sections" << endl;
    for (indexFermion = 0; indexFermion < 11; indexFermion++) {
        cout << setw(9) << ZF.convertINDF(indexFermion) << "  ";
        for (int i=0; i<=6; i++) {
            double CSA = -1.0 + 2.0/6.0*(double)i; // cos(theta)
            ZF.calcDXS(indexFermion, sqrt_s, CSA, &DXS[indexFermion]);
            cout << setw(9) << DXS[indexFermion] << "  ";
        }
        cout << endl;
    }
    cout << endl;
    
    /* Outputs from ZFITTER subroutines (for tests) */
    //
    //cout << "##### call ZFitter::info(0) #####"
    //          << endl << endl;
    //cout << " ------------- Flag info ---------------" << endl;
    //ZF.FlagInfo(); // output flag info
    //
    //cout << endl << "##### call ZFitter::info(1) #####"
    //          << endl << endl;
    //cout << " ------------------------ Cut info --------------------------";
    //ZF.CutInfo(); // print cut info

}


/* Test with ZFITTER subroutine ZFTEST() */
void test_ZFTEST(ZFitterWrapper& ZF) {
    ZF.test(0);
    //ZF.test(1);
}


void setSMparameters(StandardModel& SM_i) {
    map<string, double> Parameters;
    // 15+5 parameters defined in StandardModel
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
    //
    Parameters["phiEpsK"] = 0.0;
    Parameters["DeltaMK"] = 0.0;
    Parameters["KbarEpsK"] = 0.0;
    Parameters["Dmk"] = 0.0;
    Parameters["SM_M12D" ] = 0.0;

    // 26+16 parameters defined in QCD    
    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    Parameters["mtop"] = 174.0;
    Parameters["mbottom"] = 4.28;
    //
    Parameters["mut"] = 175.0;
    Parameters["mub"] = 5.0;
    Parameters["muc"] = 1.5;
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
    //
    Parameters["MD"] = 0.0;
    Parameters["FD"] = 0.0;
    Parameters["BD1"] = 0.0;
    Parameters["BD2"] = 0.0;
    Parameters["BD3"] = 0.0;
    Parameters["BD4"] = 0.0;
    Parameters["BD5"] = 0.0;
    Parameters["BDscale"] = 0.0;
    Parameters["BDscheme"] = 0.0;
    Parameters["BK1"] = 0.0;
    Parameters["BK2"] = 0.0;
    Parameters["BK3"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK5"] = 0.0;
    Parameters["BKscale"] = 0.0;
    Parameters["BKscheme"] = 0.0;

//    /* TEST for ZFITTER */
//    Parameters["Mz"] = 91.1876;
//    Parameters["mtop"] = 178.0;    
//    Parameters["mHl"] = 100.0;
//    Parameters["AlsMz"] = 0.117;
//    Parameters["dAle5Mz"] = 0.027572;
    
    /** To make comparisons with ZFITTER codes **/
    Parameters["GF"] = 1.16637E-5;  // for GFER=2
    Parameters["melectron"] = 0.51099907e-3;
    Parameters["mmu"] = 0.105658389;
    Parameters["mtau"] = 1.77705;
    Parameters["mup"] = 0.062;
    Parameters["mdown"] = 0.083;
    Parameters["mcharm"] = 1.50; // In ZFitter, 1.5 is the pole mass.
    Parameters["mstrange"] = 0.215;
    Parameters["mbottom"] = 4.70; // In ZFitter, 4.7 is the pole mass.
    Parameters["ale"] = 1.0/137.0359895;
    
    /* TEST for Table 6.1 in hep-ph/0507146*/
    /* flags: AMT4=6, ALEM=2 */
    Parameters["Mz"] = 91.1875;
    Parameters["mtop"] = 175.0;    
    Parameters["mHl"] = 150.0;
    Parameters["AlsMz"] = 0.118;
    Parameters["dAle5Mz"] = 0.02758;    
    
    /* to make mb(Mz) and mc(Mz) similar to ZFitter ones */
    /* mcMz = 0.56381685, mbMz = 2.8194352 */
    /* muMz = 0.062, mdMz = 0.083, msMz = 0.215 */
    Parameters["mbottom"] = 4.122;
    Parameters["mub"] = 4.122;    
    Parameters["mcharm"] = 1.171;
    Parameters["muc"] = 1.171;
    
    SM_i.Init(Parameters);
}



int main(int argc, char** argv) {

    try {    
        StandardModel* myModel;
        myModel = new StandardModel();
        myModel->InitializeModel();
        setSMparameters(*myModel);
        ZFitterWrapper ZF(*myModel);
        
        if (argc>1) {
            if (strcmp(argv[1], "-t")==0) {
                /* test with Subroutine ZFTEST */
                test_ZFTEST(ZF);
            } else {
                cout << "use the option -t for ZFTEST()" << endl;
            }
        } else {
            /* compute EW precision observables */
            test_ZFitterClass(ZF);
        }

        cout << "Test finished" << endl;

        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}

