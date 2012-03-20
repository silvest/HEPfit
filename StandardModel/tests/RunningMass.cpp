/* 
 * File:   newsimpletest.cpp
 * Author: marco
 *
 * Created on Aug 4, 2011, 1:13:13 PM
 */

#include <stdlib.h>
#include <iostream>
#include "QCD.h"

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "newsimpletest test 1" << std::endl;
    QCD qcd = QCD();
    qcd.SetQCDParameter("AlsMz", 0.119);
    qcd.SetQCDParameter("Mz", 91.1876);
    qcd.SetQCDParameter("mup", 0.0037);
    qcd.SetQCDParameter("mdown", 0.0063);
    qcd.SetQCDParameter("mstrange", 0.100);
    qcd.SetQCDParameter("mcharm", 1.3);
    qcd.SetQCDParameter("mbottom", 4.21);
    qcd.SetQCDParameter("mtop", 163.);
    qcd.SetQCDParameter("mut", 163.);
    qcd.SetQCDParameter("mub", 4.21);
    qcd.SetQCDParameter("muc", 1.3);

    std::cout << qcd.Mp2Mbar(173.3) << std::endl;
    std::cout << qcd.Mp2Mbar(174.2) << std::endl;
    std::cout << qcd.Mp2Mbar(172.4) << std::endl;
}

void test2() {
    std::cout << "newsimpletest test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (newsimpletest) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% newsimpletest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (newsimpletest)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (newsimpletest)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (newsimpletest)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (newsimpletest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

