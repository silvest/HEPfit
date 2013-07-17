/*
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LFTESTCLASS_H
#define	LFTESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "Polylogarithms.h"
#include "ClausenFunctions.h"
#include "PVfunctions.h"
//#include "LoopTools.h"
#include "../src/LoopTools.h"
using namespace std;


class LFtestclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(LFtestclass);
    CPPUNIT_TEST(Li2_m12_re);
    CPPUNIT_TEST(Li2_m12_im);
    CPPUNIT_TEST(Li2_0_re);
    CPPUNIT_TEST(Li2_0_im);
    CPPUNIT_TEST(Li2_01234_re);
    CPPUNIT_TEST(Li2_01234_im);
    CPPUNIT_TEST(Li2_12_re);
    CPPUNIT_TEST(Li2_12_im);
    CPPUNIT_TEST(Li3_m52131231311);
    CPPUNIT_TEST(Li3_m12);
    CPPUNIT_TEST(Li3_0);
    CPPUNIT_TEST(Li3_00005);
    CPPUNIT_TEST(Li3_0132);
    CPPUNIT_TEST(Li3_043456);
    CPPUNIT_TEST(Li3_05);
    CPPUNIT_TEST(Li3_05000001);
    CPPUNIT_TEST(Li3_0752);
    CPPUNIT_TEST(Li3_09999999);    
    CPPUNIT_TEST(Li3_1);    
    CPPUNIT_TEST(Cl3_0002);    
    CPPUNIT_TEST(Cl3_1);    
    CPPUNIT_TEST(Cl3_Pi);  
    CPPUNIT_TEST(A0_Mw_Mz);
    CPPUNIT_TEST(B0_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B0_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B1_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B1_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B21_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B21_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B22_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B22_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B0p_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B0p_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B1p_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B1p_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B21p_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B21p_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(B22p_Mw_Mz2_Mw_Mw_real);
    CPPUNIT_TEST(B22p_Mw_Mz2_Mw_Mw_imag);
    CPPUNIT_TEST(C0_Mz2_Mw_Mz_Mw_real);
    CPPUNIT_TEST(C0_Mz2_Mw_Mz_Mw_imag);    
    CPPUNIT_TEST(D0_s_t_Mz_0_Mz_0_real);
    CPPUNIT_TEST(D0_s_t_Mz_0_Mz_0_imag); 
    CPPUNIT_TEST(D0_s_t_Mz_0_Mz_Mt_real);
    CPPUNIT_TEST(D0_s_t_Mz_0_Mz_Mt_imag); 
    
    CPPUNIT_TEST_SUITE_END();

public:
    LFtestclass();
    virtual ~LFtestclass();
    void setUp();
    void tearDown();

private:
    Polylogarithms *myPL;
    ClausenFunctions *myClausen;
    PVfunctions *myPV;
    LoopTools *myLT;
    double epsilon;
    double Mz, Mw, mH, Mt;

    void Li2_m12_re();
    void Li2_m12_im();
    void Li2_0_re();
    void Li2_0_im();
    void Li2_01234_re();
    void Li2_01234_im();
    void Li2_12_re();
    void Li2_12_im();
    
    void Li3_m52131231311();
    void Li3_m12();
    void Li3_0();
    void Li3_00005();
    void Li3_0132();
    void Li3_043456();
    void Li3_05();
    void Li3_05000001();
    void Li3_0752();
    void Li3_09999999();
    void Li3_1();

    void Cl3_0002();
    void Cl3_1();
    void Cl3_Pi();    
    
    void A0_Mw_Mz();
    
    void B0_Mw_Mz2_Mw_Mw_real();
    void B0_Mw_Mz2_Mw_Mw_imag();
    void B1_Mw_Mz2_Mw_Mw_real();
    void B1_Mw_Mz2_Mw_Mw_imag();
    void B21_Mw_Mz2_Mw_Mw_real();
    void B21_Mw_Mz2_Mw_Mw_imag();
    void B22_Mw_Mz2_Mw_Mw_real();
    void B22_Mw_Mz2_Mw_Mw_imag();

    void B0p_Mw_Mz2_Mw_Mw_real();
    void B0p_Mw_Mz2_Mw_Mw_imag();
    void B1p_Mw_Mz2_Mw_Mw_real();
    void B1p_Mw_Mz2_Mw_Mw_imag();
    void B21p_Mw_Mz2_Mw_Mw_real();
    void B21p_Mw_Mz2_Mw_Mw_imag();
    void B22p_Mw_Mz2_Mw_Mw_real();
    void B22p_Mw_Mz2_Mw_Mw_imag();

    void C0_Mz2_Mw_Mz_Mw_real();
    void C0_Mz2_Mw_Mz_Mw_imag();    
    
    void D0_s_t_Mz_0_Mz_0_real();
    void D0_s_t_Mz_0_Mz_0_imag();    
    void D0_s_t_Mz_0_Mz_Mt_real();
    void D0_s_t_Mz_0_Mz_Mt_imag();     
    
};

#endif	/* LFTESTCLASS_H */

