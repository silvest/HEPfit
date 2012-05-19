/*
 * File:   EWSMEW2testclass.h
 * Author: mishima
 */

#ifndef EWSMEW2TESTCLASS_H
#define	EWSMEW2TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "EWSMTwoLoopEW.h"
using namespace std;


class EWSMEW2testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMEW2testclass);
    CPPUNIT_TEST(DeltaAlpha_l);

    CPPUNIT_TEST_SUITE_END();

public:
    EWSMEW2testclass();
    virtual ~EWSMEW2testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcache* myCache;
    EWSMTwoLoopEW* myEW2;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void DeltaAlpha_l();
    
};

#endif	/* EWSMEW2TESTCLASS_H */

