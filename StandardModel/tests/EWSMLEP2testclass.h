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
    CPPUNIT_TEST(GammaZTEST);
    CPPUNIT_TEST(chi_Z_real);
    CPPUNIT_TEST(chi_Z_imag);
    CPPUNIT_TEST(G1_mu_noWeak);    
    CPPUNIT_TEST(G2_mu_noWeak);    
    CPPUNIT_TEST(G3_mu_noWeak);    
    CPPUNIT_TEST(F_za_real);    
    CPPUNIT_TEST(F_za_imag);    
    
    CPPUNIT_TEST_SUITE_END();

public:
    EWSMLEP2testclass();
    virtual ~EWSMLEP2testclass();
    void setUp();
    void tearDown();

    void setModelParameters(StandardModel& Model_i);
    
private:
    StandardModel* SM;
    EWSMTwoFermionsLEP2* myLEP2;    
    
    double epsilon;
    double sqrt_s, Mw, GammaZ;
    double s, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void sqrtsTEST();
    void MwTEST();
    void GammaZTEST();
    void chi_Z_real();
    void chi_Z_imag();
    void G1_mu_noWeak();
    void G2_mu_noWeak();
    void G3_mu_noWeak();
    void F_za_real();
    void F_za_imag();
    
};

#endif	/* EWSMLEP2TESTCLASS_H */

