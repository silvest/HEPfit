/*
 * File:   EWSMEW2testrunner.cpp
 * Author: mishima
 */

#include <cppunit/BriefTestProgressListener.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestRunner.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
using namespace std;


int main() {

    /* Options for outputs */
    int prec_def = 8;
    cout.precision(prec_def);
    //int wd = prec_def + 7;
    cout.setf(ios::floatfield);
    //cout.setf(ios::scientific);
    
    try {
        
        ////////////////////////////////////////////////////////////////////            
        
        // Create the event manager and test controller
        CPPUNIT_NS::TestResult controller;
        
        // Add a listener that colllects test result
        CPPUNIT_NS::TestResultCollector result;
        controller.addListener(&result);
        
        // Add a listener that print dots as test run.
        CPPUNIT_NS::BriefTestProgressListener progress;
        controller.addListener(&progress);
        
        // Add the top suite to the test runner
        CPPUNIT_NS::TestRunner runner;
        runner.addTest(CPPUNIT_NS::TestFactoryRegistry::getRegistry().makeTest());
        runner.run(controller);
        
        // Print test in a compiler compatible format.
        CPPUNIT_NS::CompilerOutputter outputter(&result, CPPUNIT_NS::stdCOut());
        outputter.write();
        
        return result.wasSuccessful() ? 0 : 1;
    
        
        ////////////////////////////////////////////////////////////////////        

        return EXIT_SUCCESS;
    } catch (const char* c) {
        cerr << c << endl;
        return EXIT_FAILURE;
    }            
    
}
