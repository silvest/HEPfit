/*
 * File:   LF_CppUnitTestClass.cpp
 * Author: mishima
 *
 * Created on Aug 5, 2011, 12:36:51 AM
 */

#include "LF_CppUnitTestClass.h"


CPPUNIT_TEST_SUITE_REGISTRATION(LF_CppUnitTestClass);

LF_CppUnitTestClass::LF_CppUnitTestClass() {
}

LF_CppUnitTestClass::~LF_CppUnitTestClass() {
}

void LF_CppUnitTestClass::setUp() {
}

void LF_CppUnitTestClass::tearDown() {
}

void LF_CppUnitTestClass::testMethod() {
    CPPUNIT_ASSERT(true);
}

void LF_CppUnitTestClass::testFailedMethod() {
    CPPUNIT_ASSERT(false);
}

