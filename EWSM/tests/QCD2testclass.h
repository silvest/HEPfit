/*
 * File:   QCD2testclass.h
 * Author: mishima
 */

#ifndef QCD2TESTCLASS_H
#define	QCD2TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "TwoLoopQCD.h"
using namespace std;


class QCD2testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(QCD2testclass);
    CPPUNIT_TEST(F1);
    CPPUNIT_TEST(F1_0);
    CPPUNIT_TEST(V1);
    CPPUNIT_TEST(A1);
    CPPUNIT_TEST(A1_0);

    CPPUNIT_TEST_SUITE_END();

public:
    QCD2testclass();
    virtual ~QCD2testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcommon* myEWSMC;
    TwoLoopQCD* myQCD2;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void setSMparameters(StandardModel& SM_i);

    void F1();
    void F1_0();
    void V1();
    void A1();
    void A1_0();
    
    
    
    
};

#endif	/* QCD2TESTCLASS_H */

