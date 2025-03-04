/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <stdexcept>
#include "../src/StandardModel.h"
using namespace std;

//void test1() {
//    std::cout << "StandardModelTest test 1" << std::endl;
//    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
//    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
//
//    StandardModel SM(VCKM_i, .003, .007, .1, 1.5, 4.28, 174.,
//            UPMNS_i, 0.0005, 0.1, 1.7, 0., 0., 0.,
//            1.16639E-5, 0.119, 1.0/137.0360, 91.18760, 0.05907, 130.,
//            160.,4.5,1.5);
//    std::cout << SM.getMass(QCD::TOP) << std::endl;
//    std::cout << SM.getMass(StandardModel::ELECTRON) << std::endl;
//    std::cout << SM.getAlsM() << std::endl;
//
//    std::cout << "  als(4.5,0) = " << SM.als(4.5,0) << std::endl;
//    std::cout << "  lambda4 NLO = " << SM.lambda4(1) << std::endl;
//    std::cout << "  lambda4 LO = " << SM.lambda4(0) << std::endl;
//    double mbpole = SM.mbar2mp(4.5);
//    std::cout << "  mbpole = " << mbpole << std::endl;
//    std::cout << " back to mbar... " << SM.mp2mbar(mbpole) << std::endl;
//}

//void test2() {
//    std::cout << "StandardModelTest test 2" << std::endl;
//    gslpp::matrix<gslpp::complex> VCKM_i(3,3,1.);
//    gslpp::matrix<gslpp::complex> UPMNS_i(3,3,2.);
//    Parameters Par;
//    Par.Set("VCKM", VCKM_i);
//    Par.Set("mu", .003);
//    Par.Set("md", .007);
//    Par.Set("mc", 1.5);
//    Par.Set("ms", .1);
//    Par.Set("mt", 174.);
//    Par.Set("mb", 4.28);
//    Par.Set("UPMNS", UPMNS_i);
//    Par.Set("me", .5e-4);
//    Par.Set("mmu", .1);
//    Par.Set("mtau", 1.7);
//    Par.Set("mnu1", 0.);
//    Par.Set("mnu2", 0.);
//    Par.Set("mnu3", 0.);
//    Par.Set("GF", 1.16639E-5);
//    Par.Set("AlsMz", 0.119);
//    Par.Set("ale", 1.0/137.0360);
//    Par.Set("Mz", 91.18760);
//    Par.Set("dAle5Mz", 0.05907);
//    Par.Set("mHl", 130.);
//    StandardModel SM(Par);
//    std::cout << SM.getVCKM()(2,2) << std::endl;
//}

void test3() {
    WilsonCoefficient WC(8, NDR, NLO);
    WC.setMu(10.);
    std::cout << WC.getMu() << std::endl;
    WC.setMu(6.);
    std::cout << WC.getMu() << std::endl;
    WC.setMu(110.);
    std::cout << WC.getMu() << std::endl;
    WC.setMu(55.);
    std::cout << WC.getMu() << std::endl;
}

int main(int argc, char** argv) {

    try {
        //test1();

        //test2();

        test3();

        return EXIT_SUCCESS;
    } catch (const runtime_error& e) {
        cerr << e.what() << endl;
        return EXIT_FAILURE;
    }
}

