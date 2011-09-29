/*
 * File:   testclass.h
 * Author: mishima
 */

#ifndef TESTCLASS_H
#define	TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "OneLoopEW.h"
using namespace std;


class testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(testclass);
    CPPUNIT_TEST(SigmaWW_bos_Mw_0);
    CPPUNIT_TEST(SigmaWW_fer_Mw_0);
    CPPUNIT_TEST(SigmaWW_bos_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaWW_bos_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaWW_fer_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaWW_fer_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaZZ_bos_Mw_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_bos_Mw_Mz2_imag);
    CPPUNIT_TEST(SigmaZZ_fer_Mw_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_fer_Mw_Mz2_imag);
    CPPUNIT_TEST(PiZgamma_bos_Mw_Mz2_real);
    CPPUNIT_TEST(PiZgamma_bos_Mw_Mz2_imag);
    CPPUNIT_TEST(PiZgamma_fer_Mw_Mz2_real);
    CPPUNIT_TEST(PiZgamma_fer_Mw_Mz2_imag);
    CPPUNIT_TEST(PiGammaGamma_fer_Mw_0);    
    CPPUNIT_TEST(B0_diff_Mw2_Mz_Mw);
    CPPUNIT_TEST(B1_diff_Mw2_Mz_Mw);
    CPPUNIT_TEST(B21_diff_Mw2_Mz_Mw);
    CPPUNIT_TEST(Bf_diff_Mw2_Mz_Mw);
    CPPUNIT_TEST(B0p_diff_Mw2_0_Mw);
    CPPUNIT_TEST(SigmaWW_bos_diff_0_real);
    CPPUNIT_TEST(SigmaWW_fer_diff_0_real);    
    CPPUNIT_TEST(SigmaWW_bos_diff_Mw2_real);
    CPPUNIT_TEST(SigmaWW_fer_diff_Mw2_real);
    CPPUNIT_TEST(SigmaZZ_bos_diff_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_fer_diff_Mz2_real);
    CPPUNIT_TEST(PiGammaGamma_fer_diff_0_real);
    CPPUNIT_TEST(PiGammaGamma_bos_diff_Mz2_real);
    CPPUNIT_TEST(PiGammaGamma_fer_diff_Mz2_real);
    CPPUNIT_TEST(PiZgamma_bos_diff_Mz2_real);
    CPPUNIT_TEST(PiZgamma_fer_diff_Mz2_real);    
    CPPUNIT_TEST(DeltaRhobar_bos_Mw);
    CPPUNIT_TEST(DeltaRhobar_fer_Mw);
    CPPUNIT_TEST(DeltaRhobarW_bos_Mw);
    CPPUNIT_TEST(DeltaRhobarW_fer_Mw);
    CPPUNIT_TEST(DeltaRhobar_bos_Mz);
    CPPUNIT_TEST(DeltaRhobar_fer_Mz);
    CPPUNIT_TEST(TEST_DeltaRhobar_bos_Mw);
    CPPUNIT_TEST(TEST_DeltaRhobarW_bos_Mw);
    CPPUNIT_TEST(SigmaPrimeWW_bos_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaPrimeWW_bos_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaPrimeWW_fer_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaPrimeWW_fer_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaPrimeZZ_bos_Mw_Mz2_real); 
    CPPUNIT_TEST(SigmaPrimeZZ_bos_Mw_Mz2_imag); 
    CPPUNIT_TEST(SigmaPrimeZZ_fer_Mw_Mz2_real); 
    CPPUNIT_TEST(SigmaPrimeZZ_fer_Mw_Mz2_imag);
    CPPUNIT_TEST(C0_Mz2_Mt_Mw_Mt_real);
    CPPUNIT_TEST(C0_Mz2_Mt_Mw_Mt_imag);    
    CPPUNIT_TEST(C0_Mz2_0_Mw_0_real);
    CPPUNIT_TEST(C0_Mz2_0_Mw_0_imag);
    CPPUNIT_TEST(C0_Mz2_Mw_Mt_Mw_real);
    CPPUNIT_TEST(C0_Mz2_Mw_Mt_Mw_imag);
    CPPUNIT_TEST(C0_Mz2_Mw_0_Mw_real);
    CPPUNIT_TEST(C0_Mz2_Mw_0_Mw_imag);

    
    
    CPPUNIT_TEST_SUITE_END();

public:
    testclass();
    virtual ~testclass();
    void setUp();
    void tearDown();

private:
    void testMethod();
    void testFailedMethod();   
    
    StandardModel* mySM;
    EWSMcommon* myEWSMC;
    OneLoopEW* myOLEW;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void setSMparameters(StandardModel& SM_i);
    
    void SigmaWW_bos_Mw_0();
    void SigmaWW_fer_Mw_0();    
    void SigmaWW_bos_Mw_Mw2_real();
    void SigmaWW_bos_Mw_Mw2_imag();
    void SigmaWW_fer_Mw_Mw2_real();     
    void SigmaWW_fer_Mw_Mw2_imag();     
    void SigmaZZ_bos_Mw_Mz2_real();
    void SigmaZZ_bos_Mw_Mz2_imag();
    void SigmaZZ_fer_Mw_Mz2_real();
    void SigmaZZ_fer_Mw_Mz2_imag();
    void PiZgamma_bos_Mw_Mz2_real();
    void PiZgamma_bos_Mw_Mz2_imag();
    void PiZgamma_fer_Mw_Mz2_real();
    void PiZgamma_fer_Mw_Mz2_imag();
    void PiGammaGamma_fer_Mw_0();
    
    void B0_diff_Mw2_Mz_Mw();
    void B1_diff_Mw2_Mz_Mw();    
    void B21_diff_Mw2_Mz_Mw();
    void Bf_diff_Mw2_Mz_Mw();
    void B0p_diff_Mw2_0_Mw();
    
    void SigmaWW_bos_diff_0_real();    
    void SigmaWW_fer_diff_0_real();
    void SigmaWW_bos_diff_Mw2_real();
    void SigmaWW_fer_diff_Mw2_real();
    void SigmaZZ_bos_diff_Mz2_real();
    void SigmaZZ_fer_diff_Mz2_real();
    void PiGammaGamma_fer_diff_0_real();
    void PiGammaGamma_bos_diff_Mz2_real();
    void PiGammaGamma_fer_diff_Mz2_real();
    void PiZgamma_bos_diff_Mz2_real();
    void PiZgamma_fer_diff_Mz2_real();
    
    void DeltaRhobar_bos_Mw();
    void DeltaRhobar_fer_Mw();
    void DeltaRhobarW_bos_Mw();
    void DeltaRhobarW_fer_Mw();

    void DeltaRhobar_bos_Mz();
    void DeltaRhobar_fer_Mz();    
    
    void TEST_DeltaRhobar_bos_Mw();
    void TEST_DeltaRhobarW_bos_Mw();    
    
    void SigmaPrimeWW_bos_Mw_Mw2_real();
    void SigmaPrimeWW_bos_Mw_Mw2_imag();
    void SigmaPrimeWW_fer_Mw_Mw2_real();
    void SigmaPrimeWW_fer_Mw_Mw2_imag();
    void SigmaPrimeZZ_bos_Mw_Mz2_real(); 
    void SigmaPrimeZZ_bos_Mw_Mz2_imag(); 
    void SigmaPrimeZZ_fer_Mw_Mz2_real(); 
    void SigmaPrimeZZ_fer_Mw_Mz2_imag();     
    
    void C0_Mz2_Mt_Mw_Mt_real();
    void C0_Mz2_Mt_Mw_Mt_imag();
    void C0_Mz2_0_Mw_0_real();
    void C0_Mz2_0_Mw_0_imag();
    void C0_Mz2_Mw_Mt_Mw_real();
    void C0_Mz2_Mw_Mt_Mw_imag();
    void C0_Mz2_Mw_0_Mw_real();
    void C0_Mz2_Mw_0_Mw_imag();

    
    
};

#endif	/* TESTCLASS_H */

