/* 
 * File:   StandardModelTest.cpp
 * Author: silvest
 *
 * Created on Nov 30, 2010, 2:58:59 PM
 */

#include <stdlib.h>
#include <iostream>
#include "../src/StandardModel.h"

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "StandardModelTest test 1" << std::endl;
    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
    StandardModel SM(VCKM_i, .003, .007, 1.5, .1, 174., 4.28,
            UPMNS_i, 0.0005, 0.1, 1.7, 0., 0., 0.,
            1.16639E-5, 0.119, 1.0/137.0360, 91.18760, 0.05907, 130.);
    std::cout << SM.getMt() << "  " << SM.getMc() << std::endl;
    std::cout << SM.getVCKM()(2,2) << std::endl;


    std::cout << "EWPO:" << std::endl;
    std::cout << "  mW = " << SM.mW() << " GeV" << std::endl;
    std::cout << "  sin^2theta = " << SM.sin2thw() << std::endl;
    SM.setMZ(99.);
    std::cout << "  mW = " << SM.mW() << " GeV" << std::endl;
    std::cout << "  sin^2theta = " << SM.sin2thw() << std::endl;
    std::cout << "  mW = " << SM.mW() << " GeV" << std::endl;
    std::cout << "  sin^2theta = " << SM.sin2thw() << std::endl;
    double xx[1000000];
    for(int i=0; i<10000000; i++)
        xx[i%1000000]=SM.mW();

}

void test2() {
    std::cout << "StandardModelTest test 2" << std::endl;
    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
    Parameters Par;
    Par.Set("VCKM", VCKM_i);
    Par.Set("mu", .003);
    Par.Set("md", .007);
    Par.Set("mc", 1.5);
    Par.Set("ms", .1);
    Par.Set("mt", 174.);
    Par.Set("mb", 4.28);
    Par.Set("UPMNS", UPMNS_i);
    Par.Set("me", .5e-4);
    Par.Set("mmu", .1);
    Par.Set("mtau", 1.7);
    Par.Set("mnu1", 0.);
    Par.Set("mnu2", 0.);
    Par.Set("mnu3", 0.);
    Par.Set("GF", 1.16639E-5);
    Par.Set("alsMz", 0.119);
    Par.Set("ale", 1.0/137.0360);
    Par.Set("mZ", 91.18760);
    Par.Set("dAle5Mz", 0.05907);
    Par.Set("mHl", 130.);
    StandardModel SM(Par);
    std::cout << SM.getMt() << "  " << SM.getMc() << std::endl;
    std::cout << SM.getVCKM()(2,2) << std::endl;

    std::cout << "EWPO:" << std::endl;
    std::cout << "  mW = " << SM.mW() << " GeV" << std::endl;
    std::cout << "  sin^2theta = " << SM.sin2thw() << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% StandardModelTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (StandardModelTest)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (StandardModelTest)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (StandardModelTest)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (StandardModelTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

