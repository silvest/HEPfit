/* 
 * File:   ZFitterTest.cpp
 * Author: mishima
 *
 * Created on Feb 15, 2011, 11:25:35 PM
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "../src/ZFitter.h"


// GeV^{-2} --> nb
const double GeVminus2_to_nb = pow(10.0, -6.0)
                               / pow(10.0, -28.0)
                               / pow(299792458.0, -2.0)
                               / pow(6.58211899*pow(10.0,-22.0), -2.0)
                               * pow(10.0, 9.0);


/* Test for functions in Zfitter class */
void test_ZFitterClass(ZFitter& ZF) {

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

    /* sqrt(s) = mZ */
    double sqrt_s = ZF.getZMASS();
    double s = sqrt_s*sqrt_s;

    /* print input parameters */
    std::cout << std::endl << "##### call ZFitter::printInputs() #####"
              << std::endl << std::endl;
    ZF.printInputs();

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
    ZF.setAllFlags(flags, 1);

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
    ZF.setAllCuts(ICUT, ACOL, EMIN, S_PR, ANG0, ANG1, SPP, 1);
    std::cout << std::endl; 

    /* calculate EW common blocks */
    std::cout << "##### call ZFitter::calcCommonBlocks() #####" << std::endl;
    ZF.calcCommonBlocks();
    std::cout << std::endl;

    /* print input parameters */
    std::cout << std::endl << "##### call ZFitter::printInputs() #####"
              << std::endl << std::endl;
    ZF.printInputs();
    std::cout << "  Note: DAL5H is calculated in calcCommonBlocks(), "
              << "if ALEM=0 or 3. " << std::endl;
    
    /* ptint constants defiend in ZFITTER */
    std::cout << std::endl << "##### call ZFitter::printConstants() #####"
              << std::endl << std::endl;
    ZF.printConstants();

    /* print intermediate retults */
    std::cout << "##### call ZFitter::printResults() #####"
              << std::endl << std::endl;
    ZF.printResults();
    
    /* Atomic Patity Violation */
    std::cout << "##### call ZFitter::calcAPV() #####" << std::endl << std::endl;
    ZF.calcAPV(&C1U, &C1D, &C2U, &C2D);
    std::cout << "  C1U = " << C1U << " "
              << "  C1D = " << C1D << " "
              << "  C2U = " << C2U << " "
              << "  C2D = " << C2D << std::endl << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        ZF.calcXS_AFB(indexFermion, sqrt_s,
                      &XS[indexFermion], &AFB[indexFermion]);
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion)
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
        double Gamma_Z = ZF.getCommonWIDTHS(11);
        double Gamma_e = ZF.getCommonWIDTHS(1);
        double Gamma_f = Gamma_e;
        if (indexFermion!=11) Gamma_f = ZF.getCommonWIDTHS(indexFermion);
        ZF.calcXS(indexFermion, sqrt_s, Gamma_Z, Gamma_e, Gamma_f,
                  &XS[indexFermion]);
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion)
                  << "  " << XS[indexFermion]
                  << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB_2() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.getCommonWIDTHS(11);
        /* get the effective vector and axial-vector couplings */
        double GAE = sqrt(ZF.getCommonAROTFZ(1))*0.5;
        double GVE = ZF.getCommonARVEFZ(1)*GAE;
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZF.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = GAE;
            GVF = ZF.getCommonARVEFZ(1)*GAF;
        }
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion) << "  ";
        if (indexFermion!=0&&indexFermion!=10) {
            ZF.calcXS_AFB_2(indexFermion, sqrt_s, Gamma_Z, 0,
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
        double Gamma_Z = ZF.getCommonWIDTHS(11);
        /* get the effective vector and axial-vector couplings */
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZF.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = sqrt(ZF.getCommonAROTFZ(1))*0.5;
            GVF = ZF.getCommonARVEFZ(1)*GAF;
        }
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion) << "  ";
        double GVF2 = GVF*GVF;
        double GAF2 = GAF*GAF;
        if (indexFermion==1||indexFermion==2||indexFermion==3||indexFermion==11) {
            ZF.calcXS_AFB_3(indexFermion, sqrt_s, Gamma_Z, 0, GVF2, GAF2,
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
        double Gamma_Z = ZF.getCommonWIDTHS(11);
        /* get the effective vector and axial-vector couplings */
        double GAE = sqrt(ZF.getCommonAROTFZ(1))*0.5;
        double GVE = ZF.getCommonARVEFZ(1)*GAE;
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZF.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = GAE;
            GVF = ZF.getCommonARVEFZ(1)*GAF;
        }
        double PFOUR = GVE*GAE*GVF*GAF;
        double PVAE2 = GVE*GVE + GAE*GAE;
        double PVAF2 = GVF*GVF + GAF*GAF;
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion) << "  ";
        if (indexFermion==1||indexFermion==2||indexFermion==3) {
            ZF.calcXS_AFB_4(indexFermion, sqrt_s, Gamma_Z, PFOUR, PVAE2, PVAF2,
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
    ZF.calcTauPol(sqrt_s, &TAUPOL, &TAUAFB);
    std::cout << std::setw(9) << ZF.convertINDF(3) << "  "
              << TAUPOL << "  " << TAUAFB << std::endl
              << std::endl;

    std::cout << "##### call ZFitter::calcTauPol_2() #####" << std::endl
              << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Tau polarization, Tau polarization asymmetry"
              << std::endl;
    /* get the total width */
    double Gamma_Z = ZF.getCommonWIDTHS(11);
    /* get the effective vector and axial-vector couplings */
    double GAE = sqrt(ZF.getCommonAROTFZ(1))*0.5;
    double GVE = ZF.getCommonARVEFZ(1) * GAE;
    double GAF = sqrt(ZF.getCommonAROTFZ(3))*0.5;
    double GVF = ZF.getCommonARVEFZ(3) * GAF;
    ZF.calcTauPol_2(sqrt_s, Gamma_Z, 0, GVE, GAE, GVF, GAF, &TAUPOL, &TAUAFB);
    std::cout << std::setw(9) << ZF.convertINDF(3) << "  "
              << TAUPOL << "  " << TAUAFB << std::endl
              << std::endl;

    std::cout << "##### call ZFitter::calcALR() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Left-right asymmetries " << std::endl;
    double POL = 0.63;
    double ALRI;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        ZF.calcALR(indexFermion, sqrt_s, POL,
                   &XSPL[indexFermion], &XSMI[indexFermion]);
        ALRI = (XSMI[indexFermion] - XSPL[indexFermion])
                / (XSMI[indexFermion] + XSPL[indexFermion]) / POL;
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion) << "  ";
        if (indexFermion!=8) std::cout << ALRI;
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcDXS() "
              << "with the flag ISPP=2 #####" << std::endl << std::endl;
    ZF.flag("ISPP", 2); // default: 2
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  Channel, Differential cross sections" << std::endl;
    for (indexFermion = 0; indexFermion < 11; indexFermion++) {
        std::cout << std::setw(9) << ZF.convertINDF(indexFermion) << "  ";
        for (int i=0; i<=6; i++) {
            double CSA = -1.0 + 2.0/6.0*(double)i; // cos(theta)
            ZF.calcDXS(indexFermion, sqrt_s, CSA, &DXS[indexFermion]);
            std::cout << std::setw(9) << DXS[indexFermion] << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcPO() #####" << std::endl << std::endl;
    double mW, Gamma_W, sw2;
    double Gamma_inv, Gamma_had, Gamma_total;
    double sigma0_e, sigma0_mu, sigma0_tau, sigma0_had;
    double R0_e, R0_mu, R0_tau, R0_b, R0_c;
    double A_e, A_mu, A_tau, A_b, A_c;
    double AFB0_e, AFB0_mu, AFB0_tau, AFB0_b, AFB0_c;
    double s2teff_e, s2teff_mu, s2teff_tau, s2teff_b, s2teff_c;
    ZF.calcPO(&mW, &Gamma_W, &sw2,
              &Gamma_inv, &Gamma_had, &Gamma_total,
              &sigma0_e, &sigma0_mu, &sigma0_tau, &sigma0_had,
              &R0_e, &R0_mu, &R0_tau, &R0_b, &R0_c,
              &A_e, &A_mu, &A_tau, &A_b, &A_c,
              &AFB0_e, &AFB0_mu, &AFB0_tau, &AFB0_b, &AFB0_c,
              &s2teff_e, &s2teff_mu, &s2teff_tau, &s2teff_b, &s2teff_c);
    std::cout << std::setw(15) << "m_W [GeV]" << std::setw(13)<< mW << std::endl
              << std::setw(15) << "Gamma_W [GeV]" << std::setw(13) << Gamma_W << std::endl
              << std::setw(15) << "sin^2(th_W)" << std::setw(13) << sw2 << std::endl
              << std::setw(15) << "Gamma_inv [GeV]" << std::setw(13) << Gamma_inv << std::endl
              << std::setw(15) << "Gamma_had [GeV]" << std::setw(13) << Gamma_had << std::endl
              << std::setw(15) << "Gamma_Z [GeV]" << std::setw(13) << Gamma_total << std::endl
              << std::setw(15) << "sigma0_e [nb]" << std::setw(13)
              << sigma0_e*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_mu [nb]" << std::setw(13)
              << sigma0_mu*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_tau [nb]" << std::setw(13)
              << sigma0_tau*GeVminus2_to_nb << std::endl
              << std::setw(15) << "sigma0_had [nb]" << std::setw(13)
              << sigma0_had*GeVminus2_to_nb << std::endl
              << std::setw(15) << "R0_e" << std::setw(13) << R0_e << std::endl
              << std::setw(15) << "R0_mu" << std::setw(13) << R0_mu << std::endl
              << std::setw(15) << "R0_tau" << std::setw(13) << R0_tau << std::endl
              << std::setw(15) << "R0_b" << std::setw(13) << R0_b << std::endl
              << std::setw(15) << "R0_c" << std::setw(13) << R0_c << std::endl
              << std::setw(15) << "A_e" << std::setw(13) << A_e << std::endl
              << std::setw(15) << "A_mu" << std::setw(13) << A_mu << std::endl
              << std::setw(15) << "A_tau" << std::setw(13) << A_tau << std::endl
              << std::setw(15) << "A_b" << std::setw(13) << A_b << std::endl
              << std::setw(15) << "A_c" << std::setw(13) << A_c << std::endl
              << std::setw(15) << "AFB0_e" << std::setw(13) << AFB0_e << std::endl
              << std::setw(15) << "AFB0_mu" << std::setw(13) << AFB0_mu << std::endl
              << std::setw(15) << "AFB0_tau" << std::setw(13) << AFB0_tau << std::endl
              << std::setw(15) << "AFB0_b" << std::setw(13) << AFB0_b << std::endl
              << std::setw(15) << "AFB0_c" << std::setw(13) << AFB0_c << std::endl
              << std::setw(15) << "sin^2(teff_e)" << std::setw(13) << s2teff_e << std::endl
              << std::setw(15) << "sin^2(teff_mu)" << std::setw(13) << s2teff_mu << std::endl
              << std::setw(15) << "sin^2(teff_tau)" << std::setw(13) << s2teff_tau << std::endl
              << std::setw(15) << "sin^2(teff_b)" << std::setw(13) << s2teff_b << std::endl
              << std::setw(15) << "sin^2(teff_c)" << std::setw(13) << s2teff_c << std::endl
              << std::endl;

    
    /* Outputs from ZFITTER subroutines (for tests) */
    //
    //std::cout << "##### call ZFitter::info(0) #####"
    //          << std::endl << std::endl;
    //std::cout << " ------------- Flag info ---------------" << std::endl;
    //ZF.info(0); // output flag info
    //
    //std::cout << std::endl << "##### call ZFitter::info(1) #####"
    //          << std::endl << std::endl;
    //std::cout << " ------------------------ Cut info --------------------------";
    //ZF.info(1); // print cut info

}


/* Test with ZFITTER subroutine ZFTEST() */
void test_ZFTEST(ZFitter& ZF) {
    ZF.test(0);
    //ZF.test(1);
}



int main(int argc, char** argv) {

    /* set parameters */
    double mZ = 91.1876;
    double mt = 178.0;
    double mh = 100.0;
    double alphas = 0.117;
    double dalpha5h = 0.027572;
    double vtb = 1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    /* call a constructor of ZFitter */
    ZFitter ZF(mZ, mt, mh, alphas, dalpha5h, vtb,
               mu_constituent, md_constituent);

    if (argc>1) {
        if (strcmp(argv[1], "-t")==0) {
            /* test with Subroutine ZFTEST in ZFITTER */
            test_ZFTEST(ZF);
        } else {
            std::cout << "use the option -t for ZFTEST" << std::endl;
        }
    } else {
        /* compute EW precision observables */
        test_ZFitterClass(ZF);
    }
    
    return (EXIT_SUCCESS);
}

