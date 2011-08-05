/*
 * File:   LF_CppUnitTestClass.h
 * Author: mishima
 *
 * Created on Aug 5, 2011, 12:36:51 AM
 */

#ifndef LF_CPPUNITTESTCLASS_H
#define	LF_CPPUNITTESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

class LF_CppUnitTestClass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(LF_CppUnitTestClass);

    CPPUNIT_TEST(testMethod);
    CPPUNIT_TEST(testFailedMethod);

    CPPUNIT_TEST_SUITE_END();

public:
    LF_CppUnitTestClass();
    virtual ~LF_CppUnitTestClass();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();
};

#endif	/* LF_CPPUNITTESTCLASS_H */

