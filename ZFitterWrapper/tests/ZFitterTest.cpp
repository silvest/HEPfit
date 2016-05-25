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
        Parameters["tBd"] = 0.0;
    Parameters["tKp"] = 0.0;
    Parameters["tBp"] = 0.0;
    Parameters["tBs"] = 0.0;
    Parameters["tKstar"] = 0.0;
    Parameters["tphi"] = 0.0;
    Parameters["MBd"]=5.279;
    Parameters["tBd"]=1.519;
    Parameters["MBs"]=5.366; 
    Parameters["tBs"]=1.512; 
    Parameters["DGs_Gs"]=0.130; 
    Parameters["MBp"]=5.279; 
    Parameters["tBp"]=1.638; 
    Parameters["MK0"]=0.4976;
    Parameters["MKp"]=0.4937;
    Parameters["MKstar"]=0.89166; 
    Parameters["tKstar"]=1.;
    Parameters["Mphi"]=1.019461;
    Parameters["tphi"]=1.; 
    Parameters["FK"]= 0.1561;
    Parameters["FBs"]=0.2277;
    Parameters["FKstar"]=0.225; 
    Parameters["FKstarp"]= 0.185; 
    Parameters["Fphi"]=0.254; 
    Parameters["Fphip"]= 0.215; 
    Parameters["FBsoFBd"]= 1.202; 
    Parameters["csi"]= 1.06;
    Parameters["FBsSqrtBBs1"]=0.875; 
    Parameters["FBsSqrtBBs2"]=0.73;
    Parameters["FBsSqrtBBs3"]=0.89;
    Parameters["FBsSqrtBBs4"]=0.93;
    Parameters["FBsSqrtBBs5"]=1.57;
    Parameters["BBsscale"]=4.29;

    Parameters["BBsscheme"]= 0.; 
//    Parameters["lambda"]=0.2255; //HEPFit default
//    Parameters["A"]=0.82;   //HEPFit default
//    Parameters["rhob"]=0.132;   //HEPFit default
//    Parameters["etab"]=0.351;   //HEPFit default
    Parameters["lambda"] = 0.2253;  //ZFitter
    Parameters["A"] = 0.808;  //ZFitter
    Parameters["rhob"] = 0.132;  //ZFitter
    Parameters["etab"] = 0.341;  //ZFitter
    Parameters["muw"]=80.385;
    Parameters["phiEpsK"]= 43.51; 
    Parameters["KbarEpsK"]=.97;
    Parameters["DeltaMK"]= 3.483e-15; 
    Parameters["Dmk"]=1.7415e-15;
    Parameters["SM_M12D"]= 0.;
    Parameters["MD"]= 1.865;
    Parameters["FD"]= 0.212;
    Parameters["BD1"]=0.75;
    Parameters["BD2"]=0.97;
    Parameters["BD3"]=1.51;
    Parameters["BD4"]=1.22;
    Parameters["BD5"]=1.62;
    Parameters["BDscale"]= 2.;
    Parameters["BDscheme"]=2.;
    Parameters["BK1"]=0.51;
    Parameters["BK2"]=0.73;
    Parameters["BK3"]=1.29;
    Parameters["BK4"]=1.04;
    Parameters["BK5"]=0.76;
    Parameters["BKscale"]= 2.;
    Parameters["BKscheme"]=2.; 
    Parameters["EpsK"]=0.; 
    Parameters["BK(1/2)1"]=0.;
    Parameters["BK(1/2)2"]=0.;
    Parameters["BK(1/2)3"]=0.;
    Parameters["BK(1/2)4"]=0.;
    Parameters["BK(1/2)5"]=0.;
    Parameters["BK(1/2)6"]=0.;
    Parameters["BK(1/2)7"]=0.;
    Parameters["BK(1/2)8"]=0.;
    Parameters["BK(1/2)9"]=0.;
    Parameters["BK(1/2)10"]= 0.;
    Parameters["BKd_scale"]= 0.;
    Parameters["BKd_scheme"]=0.;
    Parameters["BK(3/2)1"]=0.;
    Parameters["BK(3/2)2"]=0.;
    Parameters["BK(3/2)3"]=0.;
    Parameters["BK(3/2)4"]=0.;
    Parameters["BK(3/2)5"]=0.;
    Parameters["BK(3/2)6"]=0.;
    Parameters["BK(3/2)7"]=0.;
    Parameters["BK(3/2)8"]=0.;
    Parameters["BK(3/2)9"]=0.;
    Parameters["BK(3/2)10"]= 0.; 
    Parameters["ReA0_Kd"]= 0.;
    Parameters["ReA2_Kd"]= 0.; 
    Parameters["Omega_eta_etap"]=0.; 
    Parameters["Br_Kp_P0enu"]= 0.0;
    Parameters["Br_Kp_munu"]=0.; 
    Parameters["Br_B_Xcenu"]=0.1061; 
    Parameters["DeltaP_cu"]= 0.04;
    Parameters["IB_Kl"]= 0.;
    Parameters["IB_Kp"]= 0.;
    Parameters["tKl"]=5.116e4;
    Parameters["tKp"]=1.2830e4; 

    Parameters["absh_0"]=0.000;
    Parameters["absh_p"]=0.000;
    Parameters["absh_m"]=0.000;

    Parameters["argh_0"]=0.000;
    Parameters["argh_p"]=0.000;
    Parameters["argh_m"]=0.000;

    Parameters["absh_0_1"]= 0.000;
    Parameters["absh_p_1"]= 0.000;
    Parameters["absh_m_1"]= 0.000;

    Parameters["argh_0_1"]= 0.000;
    Parameters["argh_p_1"]= 0.000;
    Parameters["argh_m_1"]= 0.000;

    Parameters["absh_0_2"]= 0.000;
    Parameters["absh_p_2"]= 0.000;
    Parameters["absh_m_2"]= 0.000;

    Parameters["argh_0_2"]= 0.000;
    Parameters["argh_p_2"]= 0.000;
    Parameters["argh_m_2"]= 0.000;

    Parameters["absh_0_MP"]=0.0005; 

    Parameters["argh_0_MP"]=0.0005; 

    Parameters["absh_0_1_MP"]=0.; 

    Parameters["argh_0_1_MP"]=0.; 

    Parameters["a_0V"]=0.365642;
    Parameters["a_1V"]=-1.08352;
    Parameters["a_2V"]=2.46546; 
    Parameters["MRV"]=5.415;
    Parameters["a_0A0"]= 0.390635;
    Parameters["a_1A0"]= -1.15441;
    Parameters["a_2A0"]= 2.08102; 
    Parameters["MRA0"]=5.366; 
    Parameters["a_0A1"]= 0.289097;
    Parameters["a_1A1"]= 0.30781;
    Parameters["a_2A1"]= 0.722586;
    Parameters["MRA1"]=5.829;
    Parameters["a_0A12"]=0.280786;
    Parameters["a_1A12"]=0.571374;
    Parameters["a_2A12"]=0.138278;
    Parameters["MRA12"]= 5.829; 
    Parameters["a_0T1"]= 0.30833; 
    Parameters["a_1T1"]= -0.961551; 
    Parameters["a_2T1"]= 2.00707; 
    Parameters["MRT1"]=5.415;
    Parameters["a_0T2"]= 0.30833; 
    Parameters["a_1T2"]= 0.423666;
    Parameters["a_2T2"]= 2.02091; 
    Parameters["MRT2"]=5.829; 
    Parameters["a_0T23"]=0.793279;
    Parameters["a_1T23"]=1.26374; 
    Parameters["a_2T23"]=1.96215; 
    Parameters["MRT23"]= 5.829; 

    Parameters["a_0Vphi"]=0.0326803; 
    Parameters["a_1Vphi"]=0.299559;
    Parameters["a_2Vphi"]=1.53848; 
    Parameters["MRVphi"]=5.415;
    Parameters["a_0A0phi"]= 0.0350432; 
    Parameters["a_1A0phi"]= 0.304956;
    Parameters["a_2A0phi"]= 1.52243; 
    Parameters["MRA0phi"]=5.366; 
    Parameters["a_0A1phi"]= 0.0271721; 
    Parameters["a_1A1phi"]= 0.218872;
    Parameters["a_2A1phi"]= 0.834282;
    Parameters["MRA1phi"]=5.829; 
    Parameters["a_0A12phi"]=0.0222117; 
    Parameters["a_1A12phi"]=0.175381;
    Parameters["a_2A12phi"]=0.999736;
    Parameters["MRA12phi"]= 5.829;
    Parameters["a_0T1phi"]= 0.0300715; 
    Parameters["a_1T1phi"]= 0.2296;
    Parameters["a_2T1phi"]= 1.3049;
    Parameters["MRT1phi"]=5.415; 
    Parameters["a_0T2phi"]= 0.0300715; 
    Parameters["a_1T2phi"]= 0.2139;
    Parameters["a_2T2phi"]= 1.13133; 
    Parameters["MRT2phi"]=5.829; 
    Parameters["a_0T23phi"]=0.0614055; 
    Parameters["a_1T23phi"]=0.473111;
    Parameters["a_2T23phi"]=2.5058;
    Parameters["MRT23phi"]= 5.829; 

    Parameters["r_1_fplus"]=0.162; 
    Parameters["r_2_fplus"]=0.173; 
    Parameters["m_fit2_fplus"]= 29.2681;
    Parameters["r_1_fT"]=0.161;
    Parameters["r_2_fT"]=0.198; 
    Parameters["m_fit2_fT"]=29.2681; 
    Parameters["r_2_f0"]=0.330;
    Parameters["m_fit2_f0"]=37.46; 

    Parameters["bsgamma_E0"]= 1.6; 
    Parameters["BLNPcorr"]= 0.; 
    Parameters["Gambino_mukin"]=1.;
    Parameters["Gambino_BRsem"]=10.67; 
    Parameters["Gambino_Mbkin"]=4.564; 
    Parameters["Gambino_Mcatmuc"]=1.087; 
    Parameters["Gambino_mupi2"]=0.470;
    Parameters["Gambino_rhoD3"]=0.171;
    Parameters["Gambino_muG2"]= 0.309;
    Parameters["Gambino_rhoLS3"]=-0.135;

    Parameters["lambdaB"]=0.350;
    Parameters["alpha1kst"]=0.2;
    Parameters["alpha2kst"]=0.05;
    Parameters["alpha2phi"]=0.;
    Parameters["alpha1kp"]=0.2;
    Parameters["alpha2kp"]=0.05;
    Parameters["delMw"]=0.;
    Parameters["delSin2th_l"]=0.;
    Parameters["delGammaZ"]=0.;
    Parameters["delR0b"]=0.;

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

