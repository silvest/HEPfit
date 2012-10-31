/*
 * Copyright (C) 2012 SUSYfit Collaboration
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
    
    CPPUNIT_TEST_SUITE_END();

public:
    LFtestclass();
    virtual ~LFtestclass();
    void setUp();
    void tearDown();

private:
    Polylogarithms *myPL;
    ClausenFunctions *myClausen;
    double epsilon;
    
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
    
};

#endif	/* LFTESTCLASS_H */

