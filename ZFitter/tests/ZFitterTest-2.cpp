/* 
 * File:   ZFitterTest-2.cpp
 * Author: mishima
 *
 * Created on Feb 21, 2011, 5:00:52 AM
 */

/*
 *  Test for Table 6.1 in hep-ph/0507146
 *
 */

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "../src/ZFitter.h"


int main(int argc, char** argv) {
 
    /* set parameters */
    double mZ = 91.1875;
    double mt = 175.0;
    double mh = 150.0;
    double alphas = 0.118;
    double dalpha5h = 0.02758;
    double vtb = 1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    /* call a constructor of ZFitter with the above parameters (essential)
     * set flags and cuts to be their defalut values */
    ZFitter ZF(mZ, mt, mh, alphas, dalpha5h, vtb,
               mu_constituent, md_constituent);

    int indexFermion;

    /* sqrt(s) = mZ */
    double sqrt_s = ZF.getZMASS();
    double s = sqrt_s*sqrt_s;

    /* print input parameters (not essential) */
    ZF.printInputs();

    /* set flags (not necessary if using the default flags) */
    // default: AMT4=4, ALEM=3
    int AFBC = 1, SCAL = 0, SCRE = 0, AMT4 = 6, BORN = 0,
        BOXD = 1, CONV = 1, FINR = 1, FOT2 = 3, GAMS = 1,
        DIAG = 1, INTF = 1, BARB = 2, PART = 0, POWR = 1,
        PRNT = 0, ALEM = 2, QCDC = 3, VPOL = 1, WEAK = 1,
        FTJR = 1, EXPR = 0, EXPF = 0, HIGS = 0, AFMT = 3,
        CZAK = 1, PREC = 10,HIG2 = 0, ALE2 = 3, GFER = 2,
        ISPP = 2, FSRS = 1, MISC = 0, MISD = 1, IPFC = 5,
        IPSC = 0, IPTO = 3, FBHO = 0, FSPP = 0, FUNA = 0,
        ASCR = 1, SFSR = 1, ENUE = 1, TUPV = 1, DMWW = 0,
        DSWW = 0;
    int flags[46] = {AFBC, SCAL, SCRE, AMT4, BORN, BOXD, CONV, FINR, FOT2, GAMS,
                     DIAG, INTF, BARB, PART, POWR, PRNT, ALEM, QCDC, VPOL, WEAK,
                     FTJR, EXPR, EXPF, HIGS, AFMT, CZAK, PREC, HIG2, ALE2, GFER,
                     ISPP, FSRS, MISC, MISD, IPFC, IPSC, IPTO, FBHO, FSPP, FUNA,
                     ASCR, SFSR, ENUE, TUPV, DMWW, DSWW};
    ZF.setAllFlags(flags, 1);

    /* set cuts (not necessary if using the default cuts) */
    /*
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
     */

    /* calculate EW common blocks (necessary before computing observables) */
    ZF.calcCommonBlocks();
    std::cout << std::endl;

    /* calcualte and print EW precision obserbables */
    std::cout << "###  see Table 6.1 in hep-ph/0507146  ###" << std::endl;
    ZF.printPO();

    return (EXIT_SUCCESS);
}

