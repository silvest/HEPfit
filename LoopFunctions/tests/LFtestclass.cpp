/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include <clooptools.h>
#include "LFtestclass.h"
#include "ClausenFunctions.h"

CPPUNIT_TEST_SUITE_REGISTRATION(LFtestclass);

LFtestclass::LFtestclass() {
    /* Initialize LoppTools library */
    ltini();
}

LFtestclass::~LFtestclass() {
}

void LFtestclass::setUp() {
    myPL = new Polylogarithms();
    myClausen = new ClausenFunctions();
    myPV = new PVfunctions(true);
    myLT = new LoopToolsWrapper();

    /* test */
    //delete myLT;
    //myLT = new LoopTools();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-10; 

    Mz = 91.1875;
    Mw = 80.360848365211552;
    Mt = 175.0;
    mh = 126.0;
    Mz2 = Mz*Mz;
    Mw2 = Mw*Mw;
    Mt2 = Mt*Mt;
    mh2 = mh*mh;
}

void LFtestclass::tearDown() {
    delete myPL;
    delete myClausen;
}

void LFtestclass::Li2_m12_re() {
    double expect = -4.650655935392819;
    double result = myPL->Li2(-12.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_m12_im() {
    double expect = 0.0;
    double result = myPL->Li2(-12.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_0_re() {
    double expect = 0.0;
    double result = myPL->Li2(0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_0_im() {
    double expect = 0.0;
    double result = myPL->Li2(0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_01234_re() {
    double expect = 0.1274314216364471;
    double result = myPL->Li2(0.1234).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_01234_im() {
    double expect = 0.0;
    double result = myPL->Li2(0.1234).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_12_re() {
    double expect = 0.1173506750161737;
    double result = myPL->Li2(12.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li2_12_im() {
    double expect = -7.806564475830408;
    double result = myPL->Li2(12.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_m52131231311() {
    double expect = -1903.63788375121;
    double result = myPL->Li3(-5213123131.1);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_m12() {
    double expect = -6.72727772529119;
    double result = myPL->Li3(-12.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_0() {
    double expect = 0.0;
    double result = myPL->Li3(0.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_00005() {
    double expect = 0.000500031254630606;
    double result = myPL->Li3(0.0005);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_0132() {
    double expect = 0.134268275005437;
    double result = myPL->Li3(0.132);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_043456() {
    double expect = 0.461929179415877;
    double result = myPL->Li3(0.43456);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_05() {
    double expect = 0.53721319360804;
    double result = myPL->Li3(0.5);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_05000001() {
    double expect = 0.537213310056148;
    double result = myPL->Li3(0.5000001);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_0752() {
    double expect = 0.847036513150411;
    double result = myPL->Li3(0.752);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_09999999() {
    double expect = 1.20205673866627;
    double result = myPL->Li3(0.9999999);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::Li3_1() {
    double expect = 1.20205690315959;
    double result = myPL->Li3(1.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

 void LFtestclass::Cl3_0002() {
     /* Mathematica: SetPrecision[Re[PolyLog[3,Exp[I 0.002]]],15] */
    double expect = 1.20204147394334;
    double result = myClausen->Cl3(0.002);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
 void LFtestclass::Cl3_1() {
    double expect = 0.448573007280017;
    double result = myClausen->Cl3(1.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
 void LFtestclass::Cl3_Pi() {
    double expect = -0.90154267736970;
    double result = myClausen->Cl3(M_PI);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
void LFtestclass::A0_Mw2_Mz2() {
    double expect = - myLT->PV_A0(Mw2, Mz2);
    double result = myPV->A0(Mw2, Mz2);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);     
} 
 
void LFtestclass::B0_Mw2_Mz2_Mw2_Mw2_real() {
    //double expect = 0.2487465514341856;
    double expect = myLT->PV_B0(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B0(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
void LFtestclass::B0_Mw2_Mz2_Mw2_Mw2_imag() {
    //double expect = 0.0;
    double expect = myLT->PV_B0(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B0(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
void LFtestclass::B1_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B1(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B1(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
void LFtestclass::B1_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B1(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B1(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}

void LFtestclass::B11_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B11(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B11(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B11_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B11(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B11(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_Mw2_Mz2_real() {
    double expect = myLT->PV_B00(Mw2, Mz2, Mw2, Mz2).real();
    double result = myPV->B00(Mw2, Mz2, Mw2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_Mw2_Mz2_imag() {
    double expect = myLT->PV_B00(Mw2, Mz2, Mw2, Mz2).imag();
    double result = myPV->B00(Mw2, Mz2, Mw2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B00(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B00(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B00(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B00(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_0_Mw2_real() {
    double expect = myLT->PV_B00(Mw2, Mz2, 0.0, Mw2).real();
    double result = myPV->B00(Mw2, Mz2, 0.0, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_0_Mw2_imag() {
    double expect = myLT->PV_B00(Mw2, Mz2, 0.0, Mw2).imag();
    double result = myPV->B00(Mw2, Mz2, 0.0, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_Mw2_0_real() {
    double expect = myLT->PV_B00(Mw2, Mz2, Mw2, 0.0).real();
    double result = myPV->B00(Mw2, Mz2, Mw2, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_Mw2_0_imag() {
    double expect = myLT->PV_B00(Mw2, Mz2, Mw2, 0.0).imag();
    double result = myPV->B00(Mw2, Mz2, Mw2, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_0_0_real() {
    double expect = myLT->PV_B00(Mw2, Mz2, 0.0, 0.0).real();
    double result = myPV->B00(Mw2, Mz2, 0.0, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_Mz2_0_0_imag() {
    double expect = myLT->PV_B00(Mw2, Mz2, 0.0, 0.0).imag();
    double result = myPV->B00(Mw2, Mz2, 0.0, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_Mw2_Mz2_real() {
    double expect = myLT->PV_B00(Mw2, 0.0, Mw2, Mz2).real();
    double result = myPV->B00(Mw2, 0.0, Mw2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_Mw2_Mz2_imag() {
    double expect = myLT->PV_B00(Mw2, 0.0, Mw2, Mz2).imag();
    double result = myPV->B00(Mw2, 0.0, Mw2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_Mw2_Mw2_real() {
    double expect = myLT->PV_B00(Mw2, 0.0, Mw2, Mw2).real();
    double result = myPV->B00(Mw2, 0.0, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_Mw2_Mw2_imag() {
    double expect = myLT->PV_B00(Mw2, 0.0, Mw2, Mw2).imag();
    double result = myPV->B00(Mw2, 0.0, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_0_Mw2_real() {
    double expect = myLT->PV_B00(Mw2, 0.0, 0.0, Mw2).real();
    double result = myPV->B00(Mw2, 0.0, 0.0, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_0_Mw2_imag() {
    double expect = myLT->PV_B00(Mw2, 0.0, 0.0, Mw2).imag();
    double result = myPV->B00(Mw2, 0.0, 0.0, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_Mw2_0_real() {
    double expect = myLT->PV_B00(Mw2, 0.0, Mw2, 0.0).real();
    double result = myPV->B00(Mw2, 0.0, Mw2, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_Mw2_0_imag() {
    double expect = myLT->PV_B00(Mw2, 0.0, Mw2, 0.0).imag();
    double result = myPV->B00(Mw2, 0.0, Mw2, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_0_0_real() {
    double expect = myLT->PV_B00(Mw2, 0.0, 0.0, 0.0).real();
    double result = myPV->B00(Mw2, 0.0, 0.0, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00_Mw2_0_0_0_imag() {
    double expect = myLT->PV_B00(Mw2, 0.0, 0.0, 0.0).imag();
    double result = myPV->B00(Mw2, 0.0, 0.0, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_Mz2_Mw2_Mz2_real() {
    double expect = myLT->PV_B0p(Mw2, Mz2, Mw2, Mz2).real();
    double result = myPV->B0p(Mw2, Mz2, Mw2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_Mz2_Mw2_Mz2_imag() {
    double expect = myLT->PV_B0p(Mw2, Mz2, Mw2, Mz2).imag();
    double result = myPV->B0p(Mw2, Mz2, Mw2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B0p(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B0p(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B0p(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B0p(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_Mw2_Mz2_real() {
    double expect = myLT->PV_B0p(Mw2, 0.0, Mw2, Mz2).real();
    double result = myPV->B0p(Mw2, 0.0, Mw2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_Mw2_Mz2_imag() {
    double expect = myLT->PV_B0p(Mw2, 0.0, Mw2, Mz2).imag();
    double result = myPV->B0p(Mw2, 0.0, Mw2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_Mw2_Mw2_real() {
    double expect = myLT->PV_B0p(Mw2, 0.0, Mw2, Mw2).real();
    double result = myPV->B0p(Mw2, 0.0, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_Mw2_Mw2_imag() {
    double expect = myLT->PV_B0p(Mw2, 0.0, Mw2, Mw2).imag();
    double result = myPV->B0p(Mw2, 0.0, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_0_Mw2_real() {
    double expect = myLT->PV_B0p(Mw2, 0.0, 0.0, Mw2).real();
    double result = myPV->B0p(Mw2, 0.0, 0.0, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_0_Mw2_imag() {
    double expect = myLT->PV_B0p(Mw2, 0.0, 0.0, Mw2).imag();
    double result = myPV->B0p(Mw2, 0.0, 0.0, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_Mw2_0_real() {
    double expect = myLT->PV_B0p(Mw2, 0.0, Mw2, 0.0).real();
    double result = myPV->B0p(Mw2, 0.0, Mw2, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_Mw2_0_imag() {
    double expect = myLT->PV_B0p(Mw2, 0.0, Mw2, 0.0).imag();
    double result = myPV->B0p(Mw2, 0.0, Mw2, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_0_0_real() {
    double expect = myLT->PV_B0p(Mw2, 0.0, 0.0, 0.0).real();
    double result = myPV->B0p(Mw2, 0.0, 0.0, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B0p_Mw2_0_0_0_imag() {
    double expect = myLT->PV_B0p(Mw2, 0.0, 0.0, 0.0).imag();
    double result = myPV->B0p(Mw2, 0.0, 0.0, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B1p_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B1p(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B1p(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B1p_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B1p(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B1p(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B11p_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B11p(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B11p(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B11p_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B11p(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B11p(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_Mz2_Mw2_Mz2_real() {
    double expect = myLT->PV_B00p(Mw2, Mz2, Mw2, Mz2).real();
    double result = myPV->B00p(Mw2, Mz2, Mw2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_Mz2_Mw2_Mz2_imag() {
    double expect = myLT->PV_B00p(Mw2, Mz2, Mw2, Mz2).imag();
    double result = myPV->B00p(Mw2, Mz2, Mw2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_Mz2_Mw2_Mw2_real() {
    double expect = myLT->PV_B00p(Mw2, Mz2, Mw2, Mw2).real();
    double result = myPV->B00p(Mw2, Mz2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_Mz2_Mw2_Mw2_imag() {
    double expect = myLT->PV_B00p(Mw2, Mz2, Mw2, Mw2).imag();
    double result = myPV->B00p(Mw2, Mz2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_Mw2_Mz2_real() {
    double expect = myLT->PV_B00p(Mw2, 0.0, Mw2, Mz2).real();
    double result = myPV->B00p(Mw2, 0.0, Mw2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_Mw2_Mz2_imag() {
    double expect = myLT->PV_B00p(Mw2, 0.0, Mw2, Mz2).imag();
    double result = myPV->B00p(Mw2, 0.0, Mw2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_Mw2_Mw2_real() {
    double expect = myLT->PV_B00p(Mw2, 0.0, Mw2, Mw2).real();
    double result = myPV->B00p(Mw2, 0.0, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_Mw2_Mw2_imag() {
    double expect = myLT->PV_B00p(Mw2, 0.0, Mw2, Mw2).imag();
    double result = myPV->B00p(Mw2, 0.0, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_0_Mw2_real() {
    double expect = myLT->PV_B00p(Mw2, 0.0, 0.0, Mw2).real();
    double result = myPV->B00p(Mw2, 0.0, 0.0, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_0_Mw2_imag() {
    double expect = myLT->PV_B00p(Mw2, 0.0, 0.0, Mw2).imag();
    double result = myPV->B00p(Mw2, 0.0, 0.0, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_Mw2_0_real() {
    double expect = myLT->PV_B00p(Mw2, 0.0, Mw2, 0.0).real();
    double result = myPV->B00p(Mw2, 0.0, Mw2, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_Mw2_0_imag() {
    double expect = myLT->PV_B00p(Mw2, 0.0, Mw2, 0.0).imag();
    double result = myPV->B00p(Mw2, 0.0, Mw2, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_0_0_real() {
    double expect = myLT->PV_B00p(Mw2, 0.0, 0.0, 0.0).real();
    double result = myPV->B00p(Mw2, 0.0, 0.0, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::B00p_Mw2_0_0_0_imag() {
    double expect = myLT->PV_B00p(Mw2, 0.0, 0.0, 0.0).imag();
    double result = myPV->B00p(Mw2, 0.0, 0.0, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_Mz2_Mw2_Mz2_Mw2_real() {
    //double expect = -7.999904668009068e-05 * (-1.0);
    double expect = - myLT->PV_C0(Mz2, Mw2, Mz2, Mw2).real();
    double result = myPV->C0(Mz2, Mw2, Mz2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::C0_Mz2_Mw2_Mz2_Mw2_imag() {
    //double expect = 0.0;
    double expect = - myLT->PV_C0(Mz2, Mw2, Mz2, Mw2).imag();
    double result = myPV->C0(Mz2, Mw2, Mz2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::C0_0_Mt2_Mz2_Mw2_real() {
    double expect = - myLT->PV_C0(0.0, Mt2, Mz2, Mw2).real();
    double result = myPV->C0(0.0, Mt2, Mz2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mt2_Mz2_Mw2_imag() {
    double expect = - myLT->PV_C0(0.0, Mt2, Mz2, Mw2).imag();
    double result = myPV->C0(0.0, Mt2, Mz2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mt2_Mw2_Mw2_real() {
    double expect = - myLT->PV_C0(0.0, Mt2, Mw2, Mw2).real();
    double result = myPV->C0(0.0, Mt2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mt2_Mw2_Mw2_imag() {
    double expect = - myLT->PV_C0(0.0, Mt2, Mw2, Mw2).imag();
    double result = myPV->C0(0.0, Mt2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mw2_Mt2_Mw2_real() {
    double expect = - myLT->PV_C0(0.0, Mw2, Mt2, Mw2).real();
    double result = myPV->C0(0.0, Mw2, Mt2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mw2_Mt2_Mw2_imag() {
    double expect = - myLT->PV_C0(0.0, Mw2, Mt2, Mw2).imag();
    double result = myPV->C0(0.0, Mw2, Mt2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mw2_Mw2_Mt2_real() {
    double expect = - myLT->PV_C0(0.0, Mw2, Mw2, Mt2).real();
    double result = myPV->C0(0.0, Mw2, Mw2, Mt2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mw2_Mw2_Mt2_imag() {
    double expect = - myLT->PV_C0(0.0, Mw2, Mw2, Mt2).imag();
    double result = myPV->C0(0.0, Mw2, Mw2, Mt2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mw2_Mw2_Mw2_real() {
    double expect = - myLT->PV_C0(0.0, Mw2, Mw2, Mw2).real();
    double result = myPV->C0(0.0, Mw2, Mw2, Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::C0_0_Mw2_Mw2_Mw2_imag() {
    double expect = - myLT->PV_C0(0.0, Mw2, Mw2, Mw2).imag();
    double result = myPV->C0(0.0, Mw2, Mw2, Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_s_t_Mz2_0_Mz2_0_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 4.020516092470443e-08;
    double expect = myLT->PV_D0(s, t, Mz2, 0.0, Mz2, 0.0).real();
    double result = myPV->D0(s, t, Mz2, 0.0, Mz2, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::D0_s_t_Mz2_0_Mz2_0_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 6.27514908604246e-08;
    double expect = myLT->PV_D0(s, t, Mz2, 0.0, Mz2, 0.0).imag();
    double result = myPV->D0(s, t, Mz2, 0.0, Mz2, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}  

void LFtestclass::D0_s_t_Mz2_0_Mz2_Mt2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 5.5248026028194e-09;
    double expect = myLT->PV_D0(s, t, Mz2, 0.0, Mz2, Mt2).real();
    double result = myPV->D0(s, t, Mz2, 0.0, Mz2, Mt2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::D0_s_t_Mz2_0_Mz2_Mt2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 3.49203127363618e-09;
    double expect = myLT->PV_D0(s, t, Mz2, 0.0, Mz2, Mt2).imag();
    double result = myPV->D0(s, t, Mz2, 0.0, Mz2, Mt2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::D0_0_0_Mt2_Mz2_Mw2_mh2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mz2, Mw2, mh2).real();
    double result = myPV->D0(0.0, 0.0, Mt2, Mz2, Mw2, mh2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mz2_Mw2_mh2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mz2, Mw2, mh2).imag();
    double result = myPV->D0(0.0, 0.0, Mt2, Mz2, Mw2, mh2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mt2_Mw2_mh2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mt2, Mw2, mh2).real();
    double result = myPV->D0(0.0, 0.0, Mt2, Mt2, Mw2, mh2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mt2_Mw2_mh2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mt2, Mw2, mh2).imag();
    double result = myPV->D0(0.0, 0.0, Mt2, Mt2, Mw2, mh2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mz2_Mt2_mh2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mz2, Mt2, mh2).real();
    double result = myPV->D0(0.0, 0.0, Mt2, Mz2, Mt2, mh2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mz2_Mt2_mh2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mz2, Mt2, mh2).imag();
    double result = myPV->D0(0.0, 0.0, Mt2, Mz2, Mt2, mh2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mz2_Mw2_Mt2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mz2, Mw2, Mt2).real();
    double result = myPV->D0(0.0, 0.0, Mt2, Mz2, Mw2, Mt2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mt2_Mz2_Mw2_Mt2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mt2, Mz2, Mw2, Mt2).imag();
    double result = myPV->D0(0.0, 0.0, Mt2, Mz2, Mw2, Mt2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mt2_Mt2_mh2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mt2, Mt2, mh2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, Mt2, Mt2, mh2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mt2_Mt2_mh2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mt2, Mt2, mh2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, Mt2, Mt2, mh2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mt2_Mw2_Mt2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mt2, Mw2, Mt2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, Mt2, Mw2, Mt2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mt2_Mw2_Mt2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mt2, Mw2, Mt2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, Mt2, Mw2, Mt2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mw2_Mt2_Mt2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mw2, Mt2, Mt2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, Mw2, Mt2, Mt2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mw2_Mt2_Mt2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mw2, Mt2, Mt2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, Mw2, Mt2, Mt2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mz2_Mz2_mh2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mz2, Mz2, mh2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, Mz2, Mz2, mh2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mz2_Mz2_mh2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mz2, Mz2, mh2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, Mz2, Mz2, mh2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mz2_mh2_Mz2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mz2, mh2, Mz2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, Mz2, mh2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mz2_mh2_Mz2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mz2, mh2, Mz2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, Mz2, mh2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_mh2_Mz2_Mz2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, mh2, Mz2, Mz2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, mh2, Mz2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_mh2_Mz2_Mz2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, mh2, Mz2, Mz2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, mh2, Mz2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_mh2_Mz2_Mz2_Mz2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, mh2, Mz2, Mz2, Mz2).real();
    double result = myPV->D0(0.0, 0.0, mh2, Mz2, Mz2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_mh2_Mz2_Mz2_Mz2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, mh2, Mz2, Mz2, Mz2).imag();
    double result = myPV->D0(0.0, 0.0, mh2, Mz2, Mz2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mz2_Mz2_Mz2_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mz2, Mz2, Mz2).real();
    double result = myPV->D0(0.0, 0.0, Mz2, Mz2, Mz2, Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

void LFtestclass::D0_0_0_Mz2_Mz2_Mz2_Mz2_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    double expect = myLT->PV_D0(0.0, 0.0, Mz2, Mz2, Mz2, Mz2).imag();
    double result = myPV->D0(0.0, 0.0, Mz2, Mz2, Mz2, Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);
}

