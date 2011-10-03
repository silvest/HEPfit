/*
 * File:   EW3testclass.h
 * Author: mishima
 */

#ifndef EW3TESTCLASS_H
#define	EW3TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "ThreeLoopEW.h"
using namespace std;


class EW3testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EW3testclass);
    CPPUNIT_TEST(DeltaAlpha_l);

    CPPUNIT_TEST_SUITE_END();

public:
    EW3testclass();
    virtual ~EW3testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcommon* myEWSMC;
    ThreeLoopEW* myEW3;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void setSMparameters(StandardModel& SM_i);
    
    void DeltaAlpha_l();
    
};

#endif	/* EW3TESTCLASS_H */

