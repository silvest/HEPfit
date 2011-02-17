/* 
 * File:   ZFitterTest.cpp
 * Author: mishima
 *
 * Created on Feb 15, 2011, 11:25:35 PM
 */

#include <stdlib.h>
#include <iostream>
#include <cstring>
#include "../src/ZFitter.h"


void test1(ZFitter ZF) {
    ZF.test(0);
}

void test2(ZFitter ZF) {

    /* set flags */
    //ZF.flag("AFBC",0);
    //ZF.info(0); // for flag info

    /* calculate EW common blocks */
    ZF.calcCommonBlocks();

    /* set cuts */
    ZF.setCuts(11, 0, 15., 10., 0., 35., 145., 0.);
    ZF.info(1); // print cut info

    double sqrt_s = ZF.getZMASS();

    int indexFermion;
    for (indexFermion=0; indexFermion<12; indexFermion++) {
        ZF.calcXS_AFB(indexFermion, sqrt_s);
        std::cout << ZF.XS[indexFermion] << "  "
                  << ZF.AFB[indexFermion] << std::endl;
    }
    std::cout << std::endl;

}


int main(int argc, char** argv) {

    double mZ =91.1876;
    double mt = 178.0;
    double mh = 100.0;
    double alphas = 0.117;
    double dalpha5h = 0.027572;
    double alphaem = 1.0/137.0359895;
    double vtb =1.0;
    double mu_constituent = 0.1;
    double md_constituent = 0.1;

    /* call a constructor of ZFitter */
    ZFitter ZF(mZ, mt, mh, alphas, dalpha5h, alphaem, vtb,
               mu_constituent, md_constituent);

    if (argc>1) {
        if (strcmp(argv[1], "-t")==0) {
            /* test with Subroutine ZFTEST in ZFITTER */
            test1(ZF);
        } else {
            std::cout << "use the option -t for ZFTEST" << std::endl;
        }
    } else {
        /* compute EW precision observables */
        test2(ZF);
    }
    
    return (EXIT_SUCCESS);
}

