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
    CPPUNIT_TEST(B_WW_d_0_hat_real);
    CPPUNIT_TEST(B_WW_d_0_hat_imag);
    CPPUNIT_TEST(rhoef_NonUnitary_MU_real);
    CPPUNIT_TEST(rhoef_NonUnitary_MU_imag);
    CPPUNIT_TEST(kappae_NonUnitary_MU_real);
    CPPUNIT_TEST(kappae_NonUnitary_MU_imag);
    CPPUNIT_TEST(kappaf_NonUnitary_MU_real);
    CPPUNIT_TEST(kappaf_NonUnitary_MU_imag);
    CPPUNIT_TEST(kappaef_NonUnitary_MU_real);
    CPPUNIT_TEST(kappaef_NonUnitary_MU_imag);
    CPPUNIT_TEST(rhoef_NonUnitary_UP_real);
    CPPUNIT_TEST(rhoef_NonUnitary_UP_imag);
    CPPUNIT_TEST(kappae_NonUnitary_UP_real);
    CPPUNIT_TEST(kappae_NonUnitary_UP_imag);
    CPPUNIT_TEST(kappaf_NonUnitary_UP_real);
    CPPUNIT_TEST(kappaf_NonUnitary_UP_imag);
    CPPUNIT_TEST(kappaef_NonUnitary_UP_real);
    CPPUNIT_TEST(kappaef_NonUnitary_UP_imag);
    CPPUNIT_TEST(Delta_rhoef_TOP_NonUnitary_real);
    CPPUNIT_TEST(Delta_rhoef_TOP_NonUnitary_imag);
    CPPUNIT_TEST(Delta_kappae_TOP_NonUnitary_real);
    CPPUNIT_TEST(Delta_kappae_TOP_NonUnitary_imag);
    CPPUNIT_TEST(Delta_kappaf_TOP_NonUnitary_real);
    CPPUNIT_TEST(Delta_kappaf_TOP_NonUnitary_imag);
    CPPUNIT_TEST(Delta_kappaef_TOP_NonUnitary_real);
    CPPUNIT_TEST(Delta_kappaef_TOP_NonUnitary_imag);
    CPPUNIT_TEST(B_WW_0_real);
    CPPUNIT_TEST(B_WW_0_imag);
    CPPUNIT_TEST(B_ZZ_0_real);
    CPPUNIT_TEST(B_ZZ_0_imag);
    CPPUNIT_TEST(Delta_rho_ef_WW_real);
    CPPUNIT_TEST(Delta_rho_ef_WW_imag);
    CPPUNIT_TEST(Delta_rho_ef_ZZ_real);
    CPPUNIT_TEST(Delta_rho_ef_ZZ_imag);
    CPPUNIT_TEST(Delta_rho_ef_WW_charm_real);
    CPPUNIT_TEST(Delta_rho_ef_WW_charm_imag);
    CPPUNIT_TEST(Delta_rho_ef_ZZ_charm_real);
    CPPUNIT_TEST(Delta_rho_ef_ZZ_charm_imag);
    CPPUNIT_TEST(G1_mu_Box_TEST);
    CPPUNIT_TEST(G2_mu_Box_TEST);
    CPPUNIT_TEST(G3_mu_Box_TEST);
    
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
    EWSMTwoFermionsLEP2* myLEP2_NU;
    
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
    void B_WW_d_0_hat_real();
    void B_WW_d_0_hat_imag();
    
    void rhoef_NonUnitary_MU_real();
    void rhoef_NonUnitary_MU_imag();
    void kappae_NonUnitary_MU_real();
    void kappae_NonUnitary_MU_imag();
    void kappaf_NonUnitary_MU_real();
    void kappaf_NonUnitary_MU_imag();
    void kappaef_NonUnitary_MU_real();
    void kappaef_NonUnitary_MU_imag();

    void rhoef_NonUnitary_UP_real();
    void rhoef_NonUnitary_UP_imag();
    void kappae_NonUnitary_UP_real();
    void kappae_NonUnitary_UP_imag();
    void kappaf_NonUnitary_UP_real();
    void kappaf_NonUnitary_UP_imag();
    void kappaef_NonUnitary_UP_real();
    void kappaef_NonUnitary_UP_imag();
    
    void Delta_rhoef_TOP_NonUnitary_real();
    void Delta_rhoef_TOP_NonUnitary_imag();
    void Delta_kappae_TOP_NonUnitary_real();
    void Delta_kappae_TOP_NonUnitary_imag();
    void Delta_kappaf_TOP_NonUnitary_real();
    void Delta_kappaf_TOP_NonUnitary_imag();
    void Delta_kappaef_TOP_NonUnitary_real();
    void Delta_kappaef_TOP_NonUnitary_imag();

    void B_WW_0_real();
    void B_WW_0_imag();
    void B_ZZ_0_real();
    void B_ZZ_0_imag();    
    void Delta_rho_ef_WW_real();
    void Delta_rho_ef_WW_imag();
    void Delta_rho_ef_ZZ_real();
    void Delta_rho_ef_ZZ_imag();

    void Delta_rho_ef_WW_charm_real();
    void Delta_rho_ef_WW_charm_imag();
    void Delta_rho_ef_ZZ_charm_real();
    void Delta_rho_ef_ZZ_charm_imag();
    
    void G1_mu_Box_TEST();
    void G2_mu_Box_TEST();
    void G3_mu_Box_TEST();
    
};

#endif	/* EWSMLEP2TESTCLASS_H */

