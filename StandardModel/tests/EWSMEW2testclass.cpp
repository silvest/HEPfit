/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
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
    mySM->InitializeModel();
    setSMparameters(*mySM);   
    myCache = new EWSMcache(*mySM);
    myEW2 = new EWSMTwoLoopEW(*myCache);
    myEWSM = new EWSM(*mySM);

    myCache->setFlagDebug(true);

    Mw = myCache->Mw(mySM->Mw_tree());/* Tests are done with the tree-level Mw */
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
    //delete myEWSM;
    delete myEW2;        
    delete myCache;
    //delete mySM;  
}

void EWSMEW2testclass::DeltaAlpha_l() {
    double ZFITTER = 0.00007761651950913550; /* ZFITTER result*/
    double result = myEW2->DeltaAlpha_l(Mz2);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMEW2testclass::rho2() {
    double ZFITTER = -6.405572408621286; /* ZFITTER result*/
    double result = myEW2->rho_2();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMEW2testclass::tau2() {
    double ZFITTER = 1.645531097693745; /* ZFITTER result*/
    double result = myEW2->tau_2();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMEW2testclass::taub() {
    double ZFITTER = -0.005678047858754-0.000033678404245; /* ZFITTER result*/
    double result = myEWSM->taub();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

// This is O(\alpha) contribution, which has to be subtracted from rho_Z^b 
// and kappa_Z^b when adding the tau_b contribution.
void EWSMEW2testclass::CORBB() {
    double ZFITTER = 0.006397903768469; /* ZFITTER result*/
    double result = myCache->ale()/8.0/M_PI/sW2*pow(myCache->Mt()/Mw, 2.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}





