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
#include "LoopToolsWrapper.h"
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
    CPPUNIT_TEST(A0_Mw2_Mz2);
    CPPUNIT_TEST(B0_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B0_Mw2_Mz2_Mw2_Mw2_imag);
    CPPUNIT_TEST(B1_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B1_Mw2_Mz2_Mw2_Mw2_imag);
    CPPUNIT_TEST(B21_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B21_Mw2_Mz2_Mw2_Mw2_imag);
    //
    CPPUNIT_TEST(B22_Mw2_Mz2_Mw2_Mz2_real);
    CPPUNIT_TEST(B22_Mw2_Mz2_Mw2_Mz2_imag);
    CPPUNIT_TEST(B22_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B22_Mw2_Mz2_Mw2_Mw2_imag);
    CPPUNIT_TEST(B22_Mw2_Mz2_0_Mw2_real);
    CPPUNIT_TEST(B22_Mw2_Mz2_0_Mw2_imag);
    CPPUNIT_TEST(B22_Mw2_Mz2_Mw2_0_real);
    CPPUNIT_TEST(B22_Mw2_Mz2_Mw2_0_imag);
    CPPUNIT_TEST(B22_Mw2_Mz2_0_0_real);
    CPPUNIT_TEST(B22_Mw2_Mz2_0_0_imag);
    CPPUNIT_TEST(B22_Mw2_0_Mw2_Mz2_real);
    CPPUNIT_TEST(B22_Mw2_0_Mw2_Mz2_imag);
    CPPUNIT_TEST(B22_Mw2_0_Mw2_Mw2_real);
    CPPUNIT_TEST(B22_Mw2_0_Mw2_Mw2_imag);
    CPPUNIT_TEST(B22_Mw2_0_0_Mw2_real);
    CPPUNIT_TEST(B22_Mw2_0_0_Mw2_imag);
    CPPUNIT_TEST(B22_Mw2_0_Mw2_0_real);
    CPPUNIT_TEST(B22_Mw2_0_Mw2_0_imag);
    CPPUNIT_TEST(B22_Mw2_0_0_0_real);
    CPPUNIT_TEST(B22_Mw2_0_0_0_imag);
    //
    CPPUNIT_TEST(B0p_Mw2_Mz2_Mw2_Mz2_real);
    CPPUNIT_TEST(B0p_Mw2_Mz2_Mw2_Mz2_imag);
    CPPUNIT_TEST(B0p_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B0p_Mw2_Mz2_Mw2_Mw2_imag);
    CPPUNIT_TEST(B0p_Mw2_0_Mw2_Mz2_real);
    CPPUNIT_TEST(B0p_Mw2_0_Mw2_Mz2_imag);
    CPPUNIT_TEST(B0p_Mw2_0_Mw2_Mw2_real);
    CPPUNIT_TEST(B0p_Mw2_0_Mw2_Mw2_imag);
    CPPUNIT_TEST(B0p_Mw2_0_0_Mw2_real);
    CPPUNIT_TEST(B0p_Mw2_0_0_Mw2_imag);
    CPPUNIT_TEST(B0p_Mw2_0_Mw2_0_real);
    CPPUNIT_TEST(B0p_Mw2_0_Mw2_0_imag);
    //CPPUNIT_TEST(B0p_Mw2_0_0_0_real);
    //CPPUNIT_TEST(B0p_Mw2_0_0_0_imag);
    //
    CPPUNIT_TEST(B1p_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B1p_Mw2_Mz2_Mw2_Mw2_imag);
    CPPUNIT_TEST(B21p_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B21p_Mw2_Mz2_Mw2_Mw2_imag);

    CPPUNIT_TEST(B22p_Mw2_Mz2_Mw2_Mz2_real);
    CPPUNIT_TEST(B22p_Mw2_Mz2_Mw2_Mz2_imag);
    CPPUNIT_TEST(B22p_Mw2_Mz2_Mw2_Mw2_real);
    CPPUNIT_TEST(B22p_Mw2_Mz2_Mw2_Mw2_imag);
    CPPUNIT_TEST(B22p_Mw2_0_Mw2_Mz2_real);
    CPPUNIT_TEST(B22p_Mw2_0_Mw2_Mz2_imag);
    CPPUNIT_TEST(B22p_Mw2_0_Mw2_Mw2_real);
    CPPUNIT_TEST(B22p_Mw2_0_Mw2_Mw2_imag);
    CPPUNIT_TEST(B22p_Mw2_0_0_Mw2_real);
    CPPUNIT_TEST(B22p_Mw2_0_0_Mw2_imag);
    CPPUNIT_TEST(B22p_Mw2_0_Mw2_0_real);
    CPPUNIT_TEST(B22p_Mw2_0_Mw2_0_imag);
    //CPPUNIT_TEST(B22p_Mw2_0_0_0_real);
    //CPPUNIT_TEST(B22p_Mw2_0_0_0_imag);

    CPPUNIT_TEST(C0_Mz2_Mw2_Mz2_Mw2_real);
    CPPUNIT_TEST(C0_Mz2_Mw2_Mz2_Mw2_imag);
    CPPUNIT_TEST(C0_0_Mt2_Mz2_Mw2_real);
    CPPUNIT_TEST(C0_0_Mt2_Mz2_Mw2_imag);
    CPPUNIT_TEST(C0_0_Mt2_Mw2_Mw2_real);
    CPPUNIT_TEST(C0_0_Mt2_Mw2_Mw2_imag);
    CPPUNIT_TEST(C0_0_Mw2_Mt2_Mw2_real);
    CPPUNIT_TEST(C0_0_Mw2_Mt2_Mw2_imag);
    CPPUNIT_TEST(C0_0_Mw2_Mw2_Mt2_real);
    CPPUNIT_TEST(C0_0_Mw2_Mw2_Mt2_imag);
    CPPUNIT_TEST(C0_0_Mw2_Mw2_Mw2_real);
    CPPUNIT_TEST(C0_0_Mw2_Mw2_Mw2_imag);

    //CPPUNIT_TEST(D0_s_t_Mz2_0_Mz2_0_real);
    //CPPUNIT_TEST(D0_s_t_Mz2_0_Mz2_0_imag);
    //CPPUNIT_TEST(D0_s_t_Mz2_0_Mz2_Mt2_real);
    //CPPUNIT_TEST(D0_s_t_Mz2_0_Mz2_Mt2_imag);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mz2_Mw2_mh2_real);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mz2_Mw2_mh2_imag);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mt2_Mw2_mh2_real);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mt2_Mw2_mh2_imag);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mz2_Mt2_mh2_real);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mz2_Mt2_mh2_imag);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mz2_Mw2_Mt2_real);
    //CPPUNIT_TEST(D0_0_0_Mt2_Mz2_Mw2_Mt2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mt2_Mt2_mh2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mt2_Mt2_mh2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mt2_Mw2_Mt2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mt2_Mw2_Mt2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mw2_Mt2_Mt2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mw2_Mt2_Mt2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mz2_Mz2_mh2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mz2_Mz2_mh2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mz2_mh2_Mz2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mz2_mh2_Mz2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_mh2_Mz2_Mz2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_mh2_Mz2_Mz2_imag);
    //CPPUNIT_TEST(D0_0_0_mh2_Mz2_Mz2_Mz2_real);
    //CPPUNIT_TEST(D0_0_0_mh2_Mz2_Mz2_Mz2_imag);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mz2_Mz2_Mz2_real);
    //CPPUNIT_TEST(D0_0_0_Mz2_Mz2_Mz2_Mz2_imag);
    
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
    LoopToolsWrapper *myLT;
    double epsilon;
    double Mz, Mw, Mt, mh, Mz2, Mw2, Mt2, mh2;

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
    
    void A0_Mw2_Mz2();
    
    void B0_Mw2_Mz2_Mw2_Mw2_real();
    void B0_Mw2_Mz2_Mw2_Mw2_imag();
    void B1_Mw2_Mz2_Mw2_Mw2_real();
    void B1_Mw2_Mz2_Mw2_Mw2_imag();
    void B21_Mw2_Mz2_Mw2_Mw2_real();
    void B21_Mw2_Mz2_Mw2_Mw2_imag();

    void B22_Mw2_Mz2_Mw2_Mz2_real();
    void B22_Mw2_Mz2_Mw2_Mz2_imag();
    void B22_Mw2_Mz2_Mw2_Mw2_real();
    void B22_Mw2_Mz2_Mw2_Mw2_imag();
    void B22_Mw2_Mz2_0_Mw2_real();
    void B22_Mw2_Mz2_0_Mw2_imag();
    void B22_Mw2_Mz2_Mw2_0_real();
    void B22_Mw2_Mz2_Mw2_0_imag();
    void B22_Mw2_Mz2_0_0_real();
    void B22_Mw2_Mz2_0_0_imag();
    void B22_Mw2_0_Mw2_Mz2_real();
    void B22_Mw2_0_Mw2_Mz2_imag();
    void B22_Mw2_0_Mw2_Mw2_real();
    void B22_Mw2_0_Mw2_Mw2_imag();
    void B22_Mw2_0_0_Mw2_real();
    void B22_Mw2_0_0_Mw2_imag();
    void B22_Mw2_0_Mw2_0_real();
    void B22_Mw2_0_Mw2_0_imag();
    void B22_Mw2_0_0_0_real();
    void B22_Mw2_0_0_0_imag();

    void B0p_Mw2_Mz2_Mw2_Mz2_real();
    void B0p_Mw2_Mz2_Mw2_Mz2_imag();
    void B0p_Mw2_Mz2_Mw2_Mw2_real();
    void B0p_Mw2_Mz2_Mw2_Mw2_imag();
    void B0p_Mw2_0_Mw2_Mz2_real();
    void B0p_Mw2_0_Mw2_Mz2_imag();
    void B0p_Mw2_0_Mw2_Mw2_real();
    void B0p_Mw2_0_Mw2_Mw2_imag();
    void B0p_Mw2_0_0_Mw2_real();
    void B0p_Mw2_0_0_Mw2_imag();
    void B0p_Mw2_0_Mw2_0_real();
    void B0p_Mw2_0_Mw2_0_imag();
    void B0p_Mw2_0_0_0_real();
    void B0p_Mw2_0_0_0_imag();

    void B1p_Mw2_Mz2_Mw2_Mw2_real();
    void B1p_Mw2_Mz2_Mw2_Mw2_imag();
    void B21p_Mw2_Mz2_Mw2_Mw2_real();
    void B21p_Mw2_Mz2_Mw2_Mw2_imag();

    void B22p_Mw2_Mz2_Mw2_Mz2_real();
    void B22p_Mw2_Mz2_Mw2_Mz2_imag();
    void B22p_Mw2_Mz2_Mw2_Mw2_real();
    void B22p_Mw2_Mz2_Mw2_Mw2_imag();
    void B22p_Mw2_0_Mw2_Mz2_real();
    void B22p_Mw2_0_Mw2_Mz2_imag();
    void B22p_Mw2_0_Mw2_Mw2_real();
    void B22p_Mw2_0_Mw2_Mw2_imag();
    void B22p_Mw2_0_0_Mw2_real();
    void B22p_Mw2_0_0_Mw2_imag();
    void B22p_Mw2_0_Mw2_0_real();
    void B22p_Mw2_0_Mw2_0_imag();
    void B22p_Mw2_0_0_0_real();
    void B22p_Mw2_0_0_0_imag();

    void C0_Mz2_Mw2_Mz2_Mw2_real();
    void C0_Mz2_Mw2_Mz2_Mw2_imag();
    void C0_0_Mt2_Mz2_Mw2_real();
    void C0_0_Mt2_Mz2_Mw2_imag();
    void C0_0_Mt2_Mw2_Mw2_real();
    void C0_0_Mt2_Mw2_Mw2_imag();
    void C0_0_Mw2_Mt2_Mw2_real();
    void C0_0_Mw2_Mt2_Mw2_imag();
    void C0_0_Mw2_Mw2_Mt2_real();
    void C0_0_Mw2_Mw2_Mt2_imag();
    void C0_0_Mw2_Mw2_Mw2_real();
    void C0_0_Mw2_Mw2_Mw2_imag();
    
    void D0_s_t_Mz2_0_Mz2_0_real();
    void D0_s_t_Mz2_0_Mz2_0_imag();
    void D0_s_t_Mz2_0_Mz2_Mt2_real();
    void D0_s_t_Mz2_0_Mz2_Mt2_imag();
    void D0_0_0_Mt2_Mz2_Mw2_mh2_real();
    void D0_0_0_Mt2_Mz2_Mw2_mh2_imag();
    void D0_0_0_Mt2_Mt2_Mw2_mh2_real();
    void D0_0_0_Mt2_Mt2_Mw2_mh2_imag();
    void D0_0_0_Mt2_Mz2_Mt2_mh2_real();
    void D0_0_0_Mt2_Mz2_Mt2_mh2_imag();
    void D0_0_0_Mt2_Mz2_Mw2_Mt2_real();
    void D0_0_0_Mt2_Mz2_Mw2_Mt2_imag();
    void D0_0_0_Mz2_Mt2_Mt2_mh2_real();
    void D0_0_0_Mz2_Mt2_Mt2_mh2_imag();
    void D0_0_0_Mz2_Mt2_Mw2_Mt2_real();
    void D0_0_0_Mz2_Mt2_Mw2_Mt2_imag();
    void D0_0_0_Mz2_Mw2_Mt2_Mt2_real();
    void D0_0_0_Mz2_Mw2_Mt2_Mt2_imag();
    void D0_0_0_Mz2_Mz2_Mz2_mh2_real();
    void D0_0_0_Mz2_Mz2_Mz2_mh2_imag();
    void D0_0_0_Mz2_Mz2_mh2_Mz2_real();
    void D0_0_0_Mz2_Mz2_mh2_Mz2_imag();
    void D0_0_0_Mz2_mh2_Mz2_Mz2_real();
    void D0_0_0_Mz2_mh2_Mz2_Mz2_imag();
    void D0_0_0_mh2_Mz2_Mz2_Mz2_real();
    void D0_0_0_mh2_Mz2_Mz2_Mz2_imag();
    void D0_0_0_Mz2_Mz2_Mz2_Mz2_real();
    void D0_0_0_Mz2_Mz2_Mz2_Mz2_imag();

};

#endif	/* LFTESTCLASS_H */

