/*
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMQCD3TESTCLASS_H
#define	EWSMQCD3TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "EWSMThreeLoopQCD.h"
using namespace std;


class EWSMQCD3testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMQCD3testclass);
    CPPUNIT_TEST(DeltaRho);
    CPPUNIT_TEST(deltaKappa_rem_l_real);
    CPPUNIT_TEST(deltaKappa_rem_q_real);    
    CPPUNIT_TEST(deltaKappa_rem_l_imag);
    CPPUNIT_TEST(deltaKappa_rem_q_imag);

    
    CPPUNIT_TEST_SUITE_END();

public:
    EWSMQCD3testclass();
    virtual ~EWSMQCD3testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcache* myCache;
    EWSMThreeLoopQCD* myQCD3;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    double AlsMt_ratio; /* The ratio of alpha_s(Mt) in ZFITTER to the current program */

    //void DeltaAlpha_t();
    void DeltaRho();
    void deltaKappa_rem_l_real();
    void deltaKappa_rem_q_real();    
    void deltaKappa_rem_l_imag();
    void deltaKappa_rem_q_imag();    
    
};

#endif	/* EWSMQCD3TESTCLASS_H */

