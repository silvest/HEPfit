/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdlib.h>
#include <iostream>
#include "../src/InputParser.h"

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "InputParserTest test 1" << std::endl;
    std::vector<Observable> Obs;
    std::vector<ModelParameter> ModPars;
    InputParser::ReadParameters("test.conf",ModPars,Obs);
    for (std::vector<ModelParameter>::iterator it=ModPars.begin() ; it < ModPars.end(); it++ )
        std::cout << *it;
    for (std::vector<Observable>::iterator xt=Obs.begin() ; xt < Obs.end(); xt++ )
        std::cout << *xt;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% InputParserTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (InputParserTest)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (InputParserTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}
