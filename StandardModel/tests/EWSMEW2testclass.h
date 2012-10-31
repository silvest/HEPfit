/*
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
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
#include "EWSM.h"
using namespace std;


class EWSMEW2testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMEW2testclass);
    CPPUNIT_TEST(DeltaAlpha_l);
    CPPUNIT_TEST(rho2);
    CPPUNIT_TEST(tau2);    
    CPPUNIT_TEST(taub);
    CPPUNIT_TEST(CORBB);    
    
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
    EWSM* myEWSM;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void DeltaAlpha_l();
    void rho2();
    void tau2();
    void taub();
    void CORBB();


};

#endif	/* EWSMEW2TESTCLASS_H */

