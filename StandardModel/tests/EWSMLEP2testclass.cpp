/*
 * File:   EWSMLEP2testclass.cpp
 * Author: mishima
 */

#include "EWSMLEP2testclass.h"
#include "EWSMsetParameters.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMLEP2testclass);

EWSMLEP2testclass::EWSMLEP2testclass() {
}

EWSMLEP2testclass::~EWSMLEP2testclass() {
}

void EWSMLEP2testclass::setUp() {
    mySM = new StandardModel(true);
    setSMparameters(*mySM);
    mySM->InitializeModel();
    myLEP2 = new EWSMTwoFermionsLEP2(*mySM);

    Mw = mySM->Mw_tree();/* Tests are done with the tree-level Mw */
    Mw2 = Mw*Mw;
    Mz = mySM->getMz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = mySM->getMtpole();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 
}

void EWSMLEP2testclass::tearDown() {
    delete myLEP2; 
    delete mySM; 
}

void EWSMLEP2testclass::MwTEST() {
    double expected = 0.0;
    double result = 100.0;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}
