/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMEW1HVTESTCLASS_H
#define	EWSMEW1HVTESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "EWSMOneLoopEW_HV.h"
using namespace std;


class EWSMEW1HVtestclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMEW1HVtestclass);
    CPPUNIT_TEST(F_Hollik_0_Mw_Mw_real);
    CPPUNIT_TEST(F_Hollik_0_Mw_Mw_imag);
    CPPUNIT_TEST(Fprime_Hollik_0_Mw_Mw_real);
    CPPUNIT_TEST(Fprime_Hollik_0_Mw_Mw_imag);
    CPPUNIT_TEST(SigmaWW_bos_Mz_Mz2_real);
    CPPUNIT_TEST(SigmaWW_bos_Mz_Mz2_imag);
    CPPUNIT_TEST(SigmaZZ_bos_Mz_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_bos_Mz_Mz2_imag);
    CPPUNIT_TEST(SigmaZgamma_bos_Mz_Mz2_real);
    CPPUNIT_TEST(SigmaZgamma_bos_Mz_Mz2_imag);
    CPPUNIT_TEST(SigmaZgamma_bos_Mz_0_real);
    CPPUNIT_TEST(SigmaZgamma_bos_Mz_0_imag);
    CPPUNIT_TEST(SigmaGammaGamma_bos_Mz_Mz2_real);
    CPPUNIT_TEST(SigmaGammaGamma_bos_Mz_Mz2_imag);
    CPPUNIT_TEST(PiGammaGamma_bos_Mz_Mz2_real);
    CPPUNIT_TEST(PiGammaGamma_bos_Mz_Mz2_imag);
    CPPUNIT_TEST(PiGammaGamma_bos_Mz_0_real);
    CPPUNIT_TEST(PiGammaGamma_bos_Mz_0_imag);
    
    CPPUNIT_TEST_SUITE_END();

public:
    EWSMEW1HVtestclass();
    virtual ~EWSMEW1HVtestclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMOneLoopEW_HV* myEW1HV;
    
    double epsilon;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    void F_Hollik_0_Mw_Mw_real();
    void F_Hollik_0_Mw_Mw_imag();
    void Fprime_Hollik_0_Mw_Mw_real();
    void Fprime_Hollik_0_Mw_Mw_imag();
    void SigmaWW_bos_Mz_Mz2_real();
    void SigmaWW_bos_Mz_Mz2_imag();
    void SigmaZZ_bos_Mz_Mz2_real();
    void SigmaZZ_bos_Mz_Mz2_imag();
    void SigmaZgamma_bos_Mz_Mz2_real();
    void SigmaZgamma_bos_Mz_Mz2_imag();
    void SigmaZgamma_bos_Mz_0_real();
    void SigmaZgamma_bos_Mz_0_imag();
    void SigmaGammaGamma_bos_Mz_Mz2_real();
    void SigmaGammaGamma_bos_Mz_Mz2_imag();
    void PiGammaGamma_bos_Mz_Mz2_real();
    void PiGammaGamma_bos_Mz_Mz2_imag();
    void PiGammaGamma_bos_Mz_0_real();
    void PiGammaGamma_bos_Mz_0_imag();
    
};

#endif	/* EWSMEW1HVTESTCLASS_H */

