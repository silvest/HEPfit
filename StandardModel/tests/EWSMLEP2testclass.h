/*
 * File:   EWSMLEP2testclass.h
 * Author: mishima
 */

#ifndef EWSMLEP2TESTCLASS_H
#define	EWSMLEP2TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "EWSMTwoFermionsLEP2.h"
using namespace std;


class EWSMLEP2testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMLEP2testclass);
    CPPUNIT_TEST(sqrtsTEST);
    CPPUNIT_TEST(MwTEST);
    
    
    CPPUNIT_TEST_SUITE_END();

public:
    EWSMLEP2testclass();
    virtual ~EWSMLEP2testclass();
    void setUp();
    void tearDown();

    void setModelParameters(StandardModel& Model_i);
    
private:
    StandardModel* mySM;
    EWSMTwoFermionsLEP2* myLEP2;
    
    double epsilon;
    double sqrt_s;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void sqrtsTEST();    
    void MwTEST();    
    
};

#endif	/* EWSMLEP2TESTCLASS_H */

