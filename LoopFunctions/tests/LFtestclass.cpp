/*
 * File:   LFtestclass.cpp
 * Author: mishima
 */

#include "LFtestclass.h"


CPPUNIT_TEST_SUITE_REGISTRATION(LFtestclass);

LFtestclass::LFtestclass() {
}

LFtestclass::~LFtestclass() {
}

void LFtestclass::setUp() {
    myPL = new Polylogarithms;

    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-10; 

}

void LFtestclass::tearDown() {
    delete myPL;
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

 