/*
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMEW3TESTCLASS_H
#define	EWSMEW3TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "EWSMThreeLoopEW.h"
using namespace std;


class EWSMEW3testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMEW3testclass);
    CPPUNIT_TEST(DeltaAlpha_l);

    CPPUNIT_TEST_SUITE_END();

public:
    EWSMEW3testclass();
    virtual ~EWSMEW3testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcache* myCache;
    EWSMThreeLoopEW* myEW3;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void DeltaAlpha_l();
    
};

#endif	/* EWSMEW3TESTCLASS_H */

