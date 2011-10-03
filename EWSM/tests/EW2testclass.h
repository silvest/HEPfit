/*
 * File:   EW2testclass.h
 * Author: mishima
 */

#ifndef EW2TESTCLASS_H
#define	EW2TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "TwoLoopEW.h"
using namespace std;


class EW2testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EW2testclass);
    CPPUNIT_TEST(DeltaAlpha_l);

    CPPUNIT_TEST_SUITE_END();

public:
    EW2testclass();
    virtual ~EW2testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcommon* myEWSMC;
    TwoLoopEW* myEW2;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void setSMparameters(StandardModel& SM_i);
    
    void DeltaAlpha_l();
    
};

#endif	/* EW2TESTCLASS_H */

