/* 
 * File:   ZFitterTest.cpp
 * Author: mishima
 *
 * Created on Feb 15, 2011, 11:25:35 PM
 */

#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <cmath>
#include "../src/ZFitter.h"


/* Test with ZFITTER subroutine ZFTEST() */
void test_ZFTEST(ZFitter& ZF) {
    ZF.test(0);
}


/* set flags (see Appendix B.2 in hep-ph/0507146) */
void setFlags(ZFitter& ZF) {
    ZF.flag("AFBC", 1);
    ZF.flag("SCAL", 0);
    ZF.flag("SCRE", 0);
    ZF.flag("AMT4", 6); // default: 4
    ZF.flag("BORN", 0);
    ZF.flag("BOXD", 1);
    ZF.flag("CONV", 1);
    ZF.flag("FINR", 1);
    ZF.flag("FOT2", 3);
    ZF.flag("GAMS", 1);
    ZF.flag("DIAG", 1);
    ZF.flag("INTF", 1);
    ZF.flag("BARB", 2);
    ZF.flag("PART", 0);
    ZF.flag("POWR", 1);
    ZF.flag("PRNT", 1); // default: 0
    ZF.flag("ALEM", 3);
    ZF.flag("QCDC", 3);
    ZF.flag("VPOL", 1);
    ZF.flag("WEAK", 1);
    ZF.flag("FTJR", 1);
    ZF.flag("EXPR", 0);
    ZF.flag("EXPF", 0);
    ZF.flag("HIGS", 0);
    ZF.flag("AFMT", 3);
    ZF.flag("CZAK", 1);
    ZF.flag("PREC", 10);
    ZF.flag("HIG2", 0);
    ZF.flag("ALE2", 3);
    ZF.flag("GFER", 2);
    ZF.flag("ISPP", 1); // default: 2
    ZF.flag("FSRS", 1);
    ZF.flag("MISC", 0);
    ZF.flag("MISD", 1);
    ZF.flag("IPFC", 5);
    ZF.flag("IPSC", 0);
    ZF.flag("IPTO", 3);
    ZF.flag("FBHO", 0);
    ZF.flag("FSPP", 0);
    ZF.flag("FUNA", 0);
    ZF.flag("ASCR", 1);
    ZF.flag("SFSR", 1);
    ZF.flag("ENUE", 1);
    ZF.flag("TUPV", 1);
    ZF.flag("DMWW", 0);
    ZF.flag("DSWW", 0);
}


/* Test with interfaces in Zfitter class */
void test_ZFitterClass(ZFitter& ZF) {

    /* set flags */
    std::cout << "##### call ZFitter::info(0) #####" << std::endl << std::endl;
    setFlags(ZF);
    ZF.info(0); // output flag info

    /* calculate EW common blocks */
    std::cout << "##### call ZFitter::calcCommonBlocks() "
              << "with the flag PRNT=1 #####" << std::endl;
    ZF.calcCommonBlocks();

    /* sqrt(s) = mZ */
    double sqrt_s = ZF.getZMASS();
    double s = sqrt_s*sqrt_s;

    /* set cuts */
    std::cout << std::endl << "##### call ZFitter::setCuts() "
              << "and ZFitter::info(1) #####" << std::endl << std::endl;
    int indexFermion;
    for (indexFermion=0; indexFermion<11; indexFermion++) {
        ZF.setCuts(indexFermion, -1, 0., 0., 0.01*s, 0., 180., 0.25*s);
    }
    ZF.setCuts(11, 0, 15., 10., 0., 35.0, 145.0, 0.);
    ZF.info(1); // print cut info

    /* Atomic Patity Violation */
    std::cout << "##### call ZFitter::calcAPV() #####" << std::endl << std::endl;
    ZF.calcAPV();
    std::cout << "  C1U = " << ZF.C1U << " "
              << "  C1D = " << ZF.C1D << " "
              << "  C2U = " << ZF.C2U << " "
              << "  C2D = " << ZF.C2D << std::endl << std::endl;

    /* Cross sections, Asymmetries and Tau polarizations at \sqrt{s}=mZ */

    std::cout << "##### call ZFitter::calcXS_AFB() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        ZF.calcXS_AFB(indexFermion, sqrt_s);
        std::cout << "  " << indexFermion << "  " << ZF.XS[indexFermion] << "  ";
        if (indexFermion!=0&&indexFermion!=10) std::cout << ZF.AFB[indexFermion];
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Cross sections" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total and partial decay widths */
        double Gamma_Z = ZF.getCommonWidths(11)*0.001;
        double Gamma_e = ZF.getCommonWidths(1)*0.001;
        double Gamma_f = Gamma_e;
        if (indexFermion!=11) Gamma_f = ZF.getCommonWidths(indexFermion)*0.001;
        ZF.calcXS(indexFermion, sqrt_s, Gamma_Z, Gamma_e, Gamma_f);
        std::cout << "  " << indexFermion << "  " << ZF.XS[indexFermion] 
                  << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB_2() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.getCommonWidths(11)*0.001;
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
        std::cout << "  " << indexFermion << "  ";
        if (indexFermion!=0&&indexFermion!=10) {
            ZF.calcXS_AFB_2(indexFermion, sqrt_s, Gamma_Z, 0,
                            GVE, GAE, GVF, GAF);
            std::cout << ZF.XS[indexFermion] << "  " << ZF.AFB[indexFermion];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcXS_AFB_3() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.getCommonWidths(11)*0.001;
        /* get the effective vector and axial-vector couplings */
        double GVF, GAF;
        if (indexFermion<11) {
            GAF = sqrt(ZF.getCommonAROTFZ(indexFermion))*0.5;
            GVF = ZF.getCommonARVEFZ(indexFermion)*GAF;
        } else {
            GAF = sqrt(ZF.getCommonAROTFZ(1))*0.5;
            GVF = ZF.getCommonARVEFZ(1)*GAF;
        }
        std::cout << "  " << indexFermion << "  ";
        double GVF2 = GVF*GVF;
        double GAF2 = GAF*GAF;
        if (indexFermion==1||indexFermion==2||indexFermion==3||indexFermion==11) {
            ZF.calcXS_AFB_3(indexFermion, sqrt_s, Gamma_Z, 0, GVF2, GAF2);
            std::cout << ZF.XS[indexFermion] << "  " << ZF.AFB[indexFermion];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

   std::cout << "##### call ZFitter::calcXS_AFB_4() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Cross sections, FB asymmetries" << std::endl;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        /* get the total width */
        double Gamma_Z = ZF.getCommonWidths(11)*0.001;
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
        std::cout << "  " << indexFermion << "  ";
        if (indexFermion==1||indexFermion==2||indexFermion==3) {
            ZF.calcXS_AFB_4(indexFermion, sqrt_s, Gamma_Z, PFOUR, PVAE2, PVAF2);
            std::cout << ZF.XS[indexFermion] << "  " << ZF.AFB[indexFermion];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcTauPol() #####" << std::endl
              << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Tau polarization, Tau polarization asymmetry"
              << std::endl;
    ZF.calcTauPol(sqrt_s);
    std::cout << "  3  " << ZF.TAUPOL << "  " << ZF.TAUAFB << std::endl
              << std::endl;

    std::cout << "##### call ZFitter::calcTauPol_2() #####" << std::endl
              << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Tau polarization, Tau polarization asymmetry"
              << std::endl;
    /* get the total width */
    double Gamma_Z = ZF.getCommonWidths(11)*0.001;
    /* get the effective vector and axial-vector couplings */
    double GAE = sqrt(ZF.getCommonAROTFZ(1))*0.5;
    double GVE = ZF.getCommonARVEFZ(1) * GAE;
    double GAF = sqrt(ZF.getCommonAROTFZ(3))*0.5;
    double GVF = ZF.getCommonARVEFZ(3) * GAF;
    ZF.calcTauPol_2(sqrt_s, Gamma_Z, 0, GVE, GAE, GVF, GAF);
    std::cout << "  3  "<< ZF.TAUPOL << "  " << ZF.TAUAFB << std::endl
              << std::endl;

    std::cout << "##### call ZFitter::calcALR() #####" << std::endl << std::endl;
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Left-right asymmetries " << std::endl;
    double POL = 0.63;
    double ALRI;
    for (indexFermion=0; indexFermion<10; indexFermion++) {
        ZF.calcALR(indexFermion, sqrt_s, POL);
        ALRI = (ZF.XSMI[indexFermion] - ZF.XSPL[indexFermion])
                / (ZF.XSMI[indexFermion] + ZF.XSPL[indexFermion]) / POL;
        std::cout << "  " << indexFermion << "  ";
        if (indexFermion!=8) std::cout << ALRI;
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "##### call ZFitter::calcDXS() "
              << "with the flag ISPP=2 #####" << std::endl;
    ZF.flag("ISPP", 2); // default: 2
    std::cout << "  sqrt(s) = " << sqrt_s << std::endl << std::endl;
    std::cout << "  INDF, Differential cross sections" << std::endl;
    for (indexFermion = 0; indexFermion < 11; indexFermion++) {
        std::cout << "  " << indexFermion << "  ";
        for (int i=0; i<=6; i++) {
            double CSA = -1.0 + 2.0/6.0*(double)i;
            ZF.calcDXS(indexFermion, sqrt_s, CSA);
            std::cout << ZF.DXS[indexFermion] << "  ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

}


int main(int argc, char** argv) {

    //std::cout.precision(8);

    double mZ = 91.1876;
    double mt = 178.0;
    double mh = 100.0;
    double alphas = 0.117;
    double dalpha5h = 0.027572;
    double alphaem = 1.0/137.0359895;
    double vtb = 1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    /* call a constructor of ZFitter */
    ZFitter ZF(mZ, mt, mh, alphas, dalpha5h, alphaem, vtb,
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

