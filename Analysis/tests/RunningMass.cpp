/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include <QCD.h>

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "RunningMass test 1" << std::endl;
    QCD qcd = QCD();
    qcd.setParameter("AlsMz", 0.1181);
    qcd.setParameter("Mz", 91.1876);
    qcd.setParameter("mup", 0.0037);
    qcd.setParameter("mdown", 0.0063);
    qcd.setParameter("mstrange", 0.100);
    qcd.setParameter("mcharm", 1.3);
    qcd.setParameter("mbottom", 4.21);
    qcd.setParameter("mtop", 163.);
    qcd.setParameter("mut", 163.);
    qcd.setParameter("mub", 4.21);
    qcd.setParameter("muc", 1.3);

    std::cout << qcd.Mp2Mbar(173.34 + 0.76) << std::endl;
    std::cout << qcd.Mp2Mbar(173.34) << std::endl;
    std::cout << qcd.Mp2Mbar(173.34 - 0.76) << std::endl;
}

void test2() {
    std::cout << "RunningMass test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (RunningMass) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% RunningMass" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (RunningMass)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (RunningMass)" << std::endl;

//    std::cout << "%TEST_STARTED% test2 (RunningMass)\n" << std::endl;
//    test2();
//    std::cout << "%TEST_FINISHED% time=0 test2 (RunningMass)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

