/*
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <cmath>
#include "LFtestclass.h"
#include "ClausenFunctions.h"


CPPUNIT_TEST_SUITE_REGISTRATION(LFtestclass);

LFtestclass::LFtestclass() {
}

LFtestclass::~LFtestclass() {
}

void LFtestclass::setUp() {
    myPL = new Polylogarithms();
    myClausen = new ClausenFunctions();
    myPV = new PVfunctions();
    myLT = new LoopTools();

    /* test */
    //delete myLT;
    //myLT = new LoopTools();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-10; 

    Mz = 91.1875;
    Mw = 80.360848365211552;
    mH = 150.0;
    Mt = 175.0;
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
 
void LFtestclass::A0_Mw_Mz() {
    double expect = myLT->PV_A0(Mw, Mz);
    double result = myPV->A0(Mw, Mz);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);     
} 
 
void LFtestclass::B0_Mw_Mz2_Mw_Mw_real() {
    //double expect = 0.2487465514341856;
    double expect = myLT->PV_B0(Mw, Mz*Mz, Mw, Mw).real();
    double result = myPV->B0(Mw, Mz*Mz, Mw, Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
void LFtestclass::B0_Mw_Mz2_Mw_Mw_imag() {
    //double expect = 0.0;
    double expect = myLT->PV_B0(Mw, Mz*Mz, Mw, Mw).imag();
    double result = myPV->B0(Mw, Mz*Mz, Mw, Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);      
}
 
void LFtestclass::C0_Mz2_Mw_Mz_Mw_real() {
    //double expect = -7.999904668009068e-05 * (-1.0);
    double expect = myLT->PV_C0(Mz*Mz, Mw, Mz, Mw).real();
    double result = myPV->C0(Mz*Mz, Mw, Mz, Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::C0_Mz2_Mw_Mz_Mw_imag() {
    //double expect = 0.0;
    double expect = myLT->PV_C0(Mz*Mz, Mw, Mz, Mw).imag();
    double result = myPV->C0(Mz*Mz, Mw, Mz, Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::D0_s_t_Mz_0_Mz_0_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 4.020516092470443e-08;
    double expect = myLT->PV_D0(s, t, Mz, 0.0, Mz, 0.0).real();
    double result = myPV->D0(s, t, Mz, 0.0, Mz, 0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::D0_s_t_Mz_0_Mz_0_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 6.27514908604246e-08;
    double expect = myLT->PV_D0(s, t, Mz, 0.0, Mz, 0.0).imag();
    double result = myPV->D0(s, t, Mz, 0.0, Mz, 0.0).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}  

void LFtestclass::D0_s_t_Mz_0_Mz_Mt_real() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 5.5248026028194e-09;
    double expect = myLT->PV_D0(s, t, Mz, 0.0, Mz, Mt).real();
    double result = myPV->D0(s, t, Mz, 0.0, Mz, Mt).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}

void LFtestclass::D0_s_t_Mz_0_Mz_Mt_imag() {
    double s = 200.0*200.0, t = 50.0*50.0;
    //double expect = 3.49203127363618e-09;
    double expect = myLT->PV_D0(s, t, Mz, 0.0, Mz, Mt).imag();
    double result = myPV->D0(s, t, Mz, 0.0, Mz, Mt).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expect, result, delta);    
}



