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
    StandardModel SM(VCKM_i, 0.15, .02, 1.5, 0.1, 170., 4.5,
        UPMNS_i, 0.0005, 0.1, 1.5, 0., 0., 0., 1.e-5);
    std::cout << SM.getMt() << "  " << SM.getMc() << std::endl;
    std::cout << SM.getVCKM()(2,2) << std::endl;
}

void test2() {
    std::cout << "StandardModelTest test 2" << std::endl;
    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
    Parameters Par;
    Par.Set("mu",.003);
    Par.Set("VCKM",VCKM_i);
    Par.Set("UPMNS",UPMNS_i);
    Par.Set("md",.007);
    Par.Set("mc",1.5);
    Par.Set("ms",.1);
    Par.Set("mt",174.);
    Par.Set("mb",4.28);
    Par.Set("me",.5e-3);
    Par.Set("mmu",.1);
    Par.Set("mtau",1.7);
    Par.Set("mnu1",1.e-10);
    Par.Set("mnu2",1.e-11);
    Par.Set("mnu3",1.e-12);
    StandardModel SM(Par);
    std::cout << SM.getMt() << "  " << SM.getMc() << std::endl;
    std::cout << SM.getVCKM()(2,2) << std::endl;
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

