/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMQCD2testclass.h"
#include "EWSMsetParameters.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMQCD2testclass);

EWSMQCD2testclass::EWSMQCD2testclass() {
}

EWSMQCD2testclass::~EWSMQCD2testclass() {
}

void EWSMQCD2testclass::setUp() {
    mySM = new StandardModel();
    mySM->InitializeModel();
    setSMparameters(*mySM);   
    myCache = new EWSMcache(*mySM);
    myQCD2 = new EWSMTwoLoopQCD(*myCache);

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

void EWSMQCD2testclass::tearDown() {
    delete myQCD2;        
    delete myCache;
    //delete mySM;  
}

void EWSMQCD2testclass::F1() {
    double ZFITTER = -1.027558816350545; /* ZFITTER result*/
    double result = myQCD2->F1(Mw*Mw/Mt/Mt, Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMQCD2testclass::F1_0() {
    double ZFITTER = -1.188052388163498; /* ZFITTER result*/
    double result = myQCD2->F1(0.0, Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMQCD2testclass::V1() {
    double ZFITTER = 0.289346402688953; /* ZFITTER result*/
    double result = myQCD2->V1(Mz*Mz/4.0/Mt/Mt);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMQCD2testclass::A1() {
    double ZFITTER = -6.747476785054744; /* ZFITTER result*/
    double result = myQCD2->A1(Mz*Mz/4.0/Mt/Mt);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMQCD2testclass::A1_0() {
    double ZFITTER = -6.897143619502220; /* ZFITTER result*/
    double result = myQCD2->A1(0.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMQCD2testclass::DeltaRho() {
    double ZFITTER = -0.000938665868188; /* ZFITTER result*/
    double result = myQCD2->DeltaRho(Mw);
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);     
}

void EWSMQCD2testclass::DeltaR_ud() {
    double ZFITTER = 0.000066523629089; /* ZFITTER result*/
    double result = myQCD2->DeltaR_ud(Mw);
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void EWSMQCD2testclass::DeltaR_tb() {
    double XTBQCD_real = 0.003850296748989; /* ZFITTER result*/
    double result = myQCD2->DeltaR_tb(Mw);
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XTBQCD_real, result, delta);     
}

void EWSMQCD2testclass::DeltaR_rem() {
    double ZFITTER = 0.000497564669985; /* ZFITTER result*/
    double result = myQCD2->DeltaR_rem(Mw);
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}
    
void EWSMQCD2testclass::DeltaRho_ud() {
    double ZFITTER = 0.000088177162158; /* ZFITTER result*/
    double result = myQCD2->DeltaRho_ud(Mw);
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}

void EWSMQCD2testclass::DeltaRho_tb() {
    double ZFITTER = -0.000899010923276; /* ZFITTER result*/
    double result = myQCD2->DeltaRho_tb(Mw);
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}

void EWSMQCD2testclass::DeltaKappa_ud_real() {
    double ZFITTER = -0.000091039011232; /* ZFITTER result*/
    double result = myQCD2->DeltaKappa_ud(Mw).real();
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}

void EWSMQCD2testclass::DeltaKappa_tb_real() {
    double ZFITTER = -0.003834690219081; /* ZFITTER result*/
    double result = myQCD2->DeltaKappa_tb(Mw).real();
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}
   
void EWSMQCD2testclass::deltaRho_rem_l_real() {
    double ZFITTER = 0.000216009269228; /* ZFITTER result*/
    double result = myQCD2->deltaRho_rem_l(mySM->NEUTRINO_1, Mw).real();
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}

void EWSMQCD2testclass::deltaRho_rem_q_real() {
    double ZFITTER = 0.000216009269228; /* ZFITTER result*/
    double result = myQCD2->deltaRho_rem_q(mySM->UP, Mw).real();
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}

void EWSMQCD2testclass::deltaKappa_rem_l_real() {
    double ZFITTER = -0.000530988904362; /* ZFITTER result*/
    double result = myQCD2->deltaKappa_rem_l(mySM->NEUTRINO_1, Mw).real();
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}

void EWSMQCD2testclass::deltaKappa_rem_q_real() {
    double ZFITTER = -0.000530988904362; /* ZFITTER result*/
    double result = myQCD2->deltaKappa_rem_q(mySM->UP, Mw).real();
    result *= AlsMt_ratio;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);        
}





