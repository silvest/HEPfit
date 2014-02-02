/*
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMEW3testclass.h"
#include "EWSMsetParameters.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMEW3testclass);

EWSMEW3testclass::EWSMEW3testclass() {
}

EWSMEW3testclass::~EWSMEW3testclass() {
}

void EWSMEW3testclass::setUp() {
    mySM = new StandardModel();
    mySM->InitializeModel();
    setSMparameters(*mySM);   
    myCache = new EWSMcache(*mySM);
    myEW3 = new EWSMThreeLoopEW(*myCache);

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

void EWSMEW3testclass::tearDown() {
    delete myEW3;        
    delete myCache;
    //delete mySM;  
}

void EWSMEW3testclass::DeltaAlpha_l() {
    double ZFITTER = 0.00000110911008841406; /* ZFITTER result*/
    double result = myEW3->DeltaAlpha_l(Mz2);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}



