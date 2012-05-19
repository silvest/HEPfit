/*
 * File:   EWSMEW2testclass.cpp
 * Author: mishima
 */

#include "EWSMEW2testclass.h"
#include "EWSMsetParameters.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMEW2testclass);

EWSMEW2testclass::EWSMEW2testclass() {
}

EWSMEW2testclass::~EWSMEW2testclass() {
}

void EWSMEW2testclass::setUp() {
    mySM = new StandardModel();
    setSMparameters(*mySM);   
    myCache = new EWSMcache(*mySM, true);
    myEW2 = new EWSMTwoLoopEW(*myCache);

    Mw = myCache->Mw(mySM->Mw_tree());
    Mw2 = Mw*Mw;
    Mz = myCache->Mz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = myCache->Mt();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 
}

void EWSMEW2testclass::tearDown() {
    delete myEW2;        
    delete myCache;
    delete mySM;  
}

void EWSMEW2testclass::DeltaAlpha_l() {
    double ZFITTER = 0.00007761651950913550; /* ZFITTER result*/
    double result = myEW2->DeltaAlpha_l();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}







