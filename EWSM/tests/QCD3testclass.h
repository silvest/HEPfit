/*
 * File:   QCD3testclass.h
 * Author: mishima
 */

#ifndef QCD3TESTCLASS_H
#define	QCD3TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "ThreeLoopQCD.h"
using namespace std;


class QCD3testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(QCD3testclass);
    CPPUNIT_TEST(DeltaRho);
    CPPUNIT_TEST(deltaKappa_rem_l_real);
    CPPUNIT_TEST(deltaKappa_rem_q_real);    
    CPPUNIT_TEST(deltaKappa_rem_l_imag);
    CPPUNIT_TEST(deltaKappa_rem_q_imag);

    
    CPPUNIT_TEST_SUITE_END();

public:
    QCD3testclass();
    virtual ~QCD3testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcommon* myEWSMC;
    ThreeLoopQCD* myQCD3;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    double AlsMt_ratio; /* The ratio of alpha_s(Mt) in ZFITTER to the current program */
                         
    void setSMparameters(StandardModel& SM_i);

    //void DeltaAlpha_t();
    void DeltaRho();
    void deltaKappa_rem_l_real();
    void deltaKappa_rem_q_real();    
    void deltaKappa_rem_l_imag();
    void deltaKappa_rem_q_imag();    
    
};

#endif	/* QCD3TESTCLASS_H */

