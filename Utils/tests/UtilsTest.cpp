/* 
 * File:   UtilsTest.cpp
 * Author: marco
 *
 * Created on Feb 2, 2011, 5:32:16 PM
 */

#include <stdlib.h>
#include <iostream>
#include "../src/Parameters.h"
#include "../../gslpp/src/gslpp_complex.h"

/*
 * Simple C++ Test Suite
 */

void test1() {
    std::cout << "UtilsTest test 1" << std::endl;

  int i;
  for(i=0; i<Parameters::NumberOfTypes; i++)
    std::cout << Parameters::TypeList[i] << std::endl;

  Parameters P;
  int j;
  double d;
  std::string s;
  gslpp::complex z;


  P.Set("Mt",180.);
  P.Set("N",10);
  P.Set("Model","SM");
  P.Set("Mt",170.);
  P.Set("mu",10.+5.*gslpp::complex::i());

  P.Get("Mt",d);
  P.Get("N",i);
  P.Get("Model",s);
  P.Get("mu",z);
  std::cout << i << std::endl;
  std::cout << d << std::endl;
  std::cout << s << std::endl;
  std::cout << z << std::endl;
// errors
//  P.Set("Mt",1);
//  P.Get("Mt",j);
//  P.Get("Mq",j);
}

void test2() {
    std::cout << "UtilsTest test 2" << std::endl;
    std::cout << "%TEST_FAILED% time=0 testname=test2 (UtilsTest) message=error message sample" << std::endl;
}

int main(int argc, char** argv) {
    std::cout << "%SUITE_STARTING% UtilsTest" << std::endl;
    std::cout << "%SUITE_STARTED%" << std::endl;

    std::cout << "%TEST_STARTED% test1 (UtilsTest)" << std::endl;
    test1();
    std::cout << "%TEST_FINISHED% time=0 test1 (UtilsTest)" << std::endl;

    std::cout << "%TEST_STARTED% test2 (UtilsTest)\n" << std::endl;
    test2();
    std::cout << "%TEST_FINISHED% time=0 test2 (UtilsTest)" << std::endl;

    std::cout << "%SUITE_FINISHED% time=0" << std::endl;

    return (EXIT_SUCCESS);
}

