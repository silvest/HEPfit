/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include "../src/MFV.h"

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "MFVtest test 1" << std::endl;
    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
    MFV MFVM(VCKM_i, 0.15, .02, 1.5, 0.1, 170., 4.5,
        UPMNS_i, 0.0005, 0.1, 1.5, 0., 0., 0., 10.,1.);
  
    MFVM.setParameter(1.E4, 1.e4, 1.e4, 1.e2, 1.e2, 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0.);
    gslpp::matrix<gslpp::complex> pippo(MFVM.getMsu2());
    std::cout << MFVM.getRu()*pippo*MFVM.getRu().hconjugate() << std::endl;
    std::cout << pippo << std::endl;
    //std::cout << MFVM.getMsu2() << std::endl;
    //std::cout << MFVM.getRd() << std::endl;
    //std::cout << MFVM.getMsd2() << std::endl;
}

void test2() {
    std::cout << "MFVtest test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (MFVtest) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% MFVtest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (MFVtest)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (MFVtest)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (MFVtest)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (MFVtest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

