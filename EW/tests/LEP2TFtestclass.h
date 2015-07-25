/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef LEP2TFTESTCLASS_H
#define	LEP2TFTESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "LEP2TwoFermions.h"
using namespace std;

const double GeVminus2_to_nb = 389379.338;

class LEP2TFtestclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(LEP2TFtestclass);
    CPPUNIT_TEST(sqrtsTEST);
    CPPUNIT_TEST(MwTEST);
    CPPUNIT_TEST(GammaZTEST);
    CPPUNIT_TEST(sigma_mu);   

    CPPUNIT_TEST_SUITE_END();

public:
    LEP2TFtestclass();
    virtual ~LEP2TFtestclass();
    void setUp();
    void tearDown();

    void setModelParameters(StandardModel& Model_i);
    
private:
    StandardModel* SM;
    LEP2TwoFermions* myLEP2TF;
    
    double epsilon;
    double sqrt_s, Mw, GammaZ;
    double s, Mw2, Mz, Mz2, cW2, sW2, Mt;
    
    static const bool noRCs[5], withRCs[5];
    
    void sqrtsTEST();
    void MwTEST();
    void GammaZTEST();
    void sigma_mu();

};

#endif	/* LEP2TFTESTCLASS_H */

