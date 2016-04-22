/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMQCD3testclass.h"
#include "EWSMsetParameters.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMQCD3testclass);

EWSMQCD3testclass::EWSMQCD3testclass() {
}

EWSMQCD3testclass::~EWSMQCD3testclass() {
}

void EWSMQCD3testclass::setUp() {
    mySM = new StandardModel();
    mySM->InitializeModel();
    setSMparameters(*mySM);   
    myCache = new EWSMcache(*mySM);
    myQCD3 = new EWSMThreeLoopQCD(*myCache);

    myCache->setFlagDebug(true);

    Mw = myCache->Mw(mySM->Mw_tree());/* Tests are done with the tree-level Mw */
    Mw2 = Mw*Mw;
    Mz = myCache->Mz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = myCache->Mt();
    
    AlsMt_ratio = 0.107443278759216/myCache->alsMt();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 
}

void EWSMQCD3testclass::tearDown() {
    delete myQCD3;
    delete myCache;
    //delete mySM;  
}

void EWSMQCD3testclass::DeltaRho() {
    double ZFITTER = -0.000204951302281; /* ZFITTER result*/
    double result = myQCD3->DeltaRho(Mw);
    result *= AlsMt_ratio*AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);     
}

void EWSMQCD3testclass::deltaKappa_rem_l_real() {
    double ZFITTER = -0.000025563165772; /* ZFITTER result*/
    double result = myQCD3->deltaKappa_rem_l(mySM->NEUTRINO_1, Mw).real();
    result *= AlsMt_ratio*AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);     
}

void EWSMQCD3testclass::deltaKappa_rem_q_real() {
    double ZFITTER = -0.000025563165772; /* ZFITTER result*/
    double result = myQCD3->deltaKappa_rem_q(mySM->UP, Mw).real();
    result *= AlsMt_ratio*AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);     
}

void EWSMQCD3testclass::deltaKappa_rem_l_imag() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myQCD3->deltaKappa_rem_l(mySM->NEUTRINO_1, Mw).imag();
    result *= AlsMt_ratio*AlsMt_ratio;
    //double delta = fabs(epsilon*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);     
}

void EWSMQCD3testclass::deltaKappa_rem_q_imag() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myQCD3->deltaKappa_rem_q(mySM->UP, Mw).imag();
    result *= AlsMt_ratio*AlsMt_ratio;
    //double delta = fabs(epsilon*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);     
}







