/* 
 * Copyright (C) 2012-2014 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <complex>

#include "EWSMLEP2testclass.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMLEP2testclass);

EWSMLEP2testclass::EWSMLEP2testclass() {
}

EWSMLEP2testclass::~EWSMLEP2testclass() {
}

void EWSMLEP2testclass::setUp() {
    SM = new StandardModel();
    SM->InitializeModel();
    setModelParameters(*SM);
    myCache = new EWSMcache(*SM);
    myLEP2 = new EWSMTwoFermionsLEP2(*SM, *myCache);
    myLEP2_NU = new EWSMTwoFermionsLEP2(*SM, *myCache);

    myCache->setFlagDebug(true);
    myLEP2->setBDebug(true);
    myLEP2_NU->setBDebug(true);

    sqrt_s = 200.0;
    Mw = 80.360848365211552;
    GammaZ = 2.494980275134754;
    // alpha(s)=0.007821563447202 or 0.007821563459046
    
    s = sqrt_s*sqrt_s;
    Mw2 = Mw*Mw;
    Mz = SM->getMz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = SM->getMtpole();

    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 
}

void EWSMLEP2testclass::tearDown() {
    delete myLEP2;
    delete myCache;
    //delete SM; 
}

/*  Parameters for StandardModel class  */
void EWSMLEP2testclass::setModelParameters(StandardModel& Model_i) {
    std::map<std::string, double> Parameters;
    
    // 17+5 parameters defined in StandardModel
    //Parameters["GF"] = 1.16637E-5;
    Parameters["mneutrino_1"] = 0.0;
    Parameters["mneutrino_2"] = 0.0;
    Parameters["mneutrino_3"] = 0.0;
    //Parameters["melectron"] = 0.54857990943e-3;
    //Parameters["mmu"] = 0.1134289256;
    //Parameters["mtau"] = 1.77682;
    Parameters["lambda"] = 0.2253;
    Parameters["A"] = 0.808;
    Parameters["rhob"] = 0.132;
    Parameters["etab"] = 0.341;
    //Parameters["ale"] = 1.0/137.035999679;
    //Parameters["dAle5Mz"] = 0.02758;
    //Parameters["mHl"] = 130.0;
    Parameters["muw"] = 80.0;
    //Parameters["mub"] = 4.2;
    //Parameters["muc"] = 1.4;
    //
    Parameters["phiEpsK"] = 0.0;
    Parameters["DeltaMK"] = 0.0;
    Parameters["KbarEpsK"] = 0.0;
    Parameters["Dmk"] = 0.0;
    Parameters["SM_M12D" ] = 0.0;
    
    // 26+16+1 parameters defined in QCD    
    //Parameters["AlsMz"] = 0.1184;
    //Parameters["Mz"] = 91.1876;
    //Parameters["mup"] = 0.003;
    //Parameters["mdown"] = 0.007;
    //Parameters["mcharm"] = 1.5;
    //Parameters["mstrange"] = 0.1;
    //Parameters["mtop"] = 174.0;
    //Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 175.0;
    Parameters["mub"] = 4.2;
    Parameters["muc"] = 1.4;
    Parameters["MBd"] = 0.0;
    Parameters["MBs"] = 0.0;
    Parameters["MBp"] = 0.0;
    Parameters["MK0"] = 0.0;
    Parameters["MKp"] = 0.0;
    Parameters["FBs"] = 0.0;
    Parameters["FBsoFBd"] = 0.0;
    Parameters["BBsoBBd"] = 0.0;
    Parameters["BBs1"] = 0.0;
    Parameters["BBs2"] = 0.0;  
    Parameters["BBs3"] = 0.0;
    Parameters["BBs4"] = 0.0;
    Parameters["BBs5"] = 0.0;
    Parameters["BBsscale"] = 0.0;
    Parameters["BBsscheme"] = 0.0;
    //
    Parameters["MD"] = 0.0;
    Parameters["FD"] = 0.0;
    Parameters["BD1"] = 0.0;
    Parameters["BD2"] = 0.0;
    Parameters["BD3"] = 0.0;
    Parameters["BD4"] = 0.0;
    Parameters["BD5"] = 0.0;
    Parameters["BDscale"] = 0.0;
    Parameters["BDscheme"] = 0.0;
    Parameters["BK1"] = 0.0;
    Parameters["BK2"] = 0.0;
    Parameters["BK3"] = 0.0;
    Parameters["BK4"] = 0.0;
    Parameters["BK5"] = 0.0;
    Parameters["BKscale"] = 0.0;
    Parameters["BKscheme"] = 0.0; 
    Parameters["FK"] = 0.0; 
    
    /** To make comparisons with ZFITTER codes **/
    Parameters["GF"] = 1.16637E-5;  // for GFER=2
    Parameters["melectron"] = 0.51099907e-3;
    Parameters["mmu"] = 0.105658389;
    Parameters["mtau"] = 1.77705;
    Parameters["mup"] = 0.062;
    Parameters["mdown"] = 0.083;
    //Parameters["mcharm"] = 1.50; // In ZFitter, 1.5 is the pole mass, but we input it to loop functions by hand for comparisons. 
    Parameters["mstrange"] = 0.215;
    //Parameters["mbottom"] = 4.70; // In ZFitter, 4.7 is the pole mass, but we input it to loop functions by hand for comparisons. 
    Parameters["ale"] = 1.0/137.0359895;

    Parameters["Mz"] = 91.1875;
    Parameters["AlsMz"] = 0.118;
    Parameters["dAle5Mz"] = 0.02758;
    Parameters["mHl"] = 125.7;        
    Parameters["mtop"] = 173.18;
    
    /* mb(Mz)= 0.56381685 and mc(Mz)= 2.8194352 in ZFitter */
    /* In case where bDebug is true, the RG running effects of the quark masses are neglected. */
    Parameters["mcharm"] = 0.56381685;
    Parameters["mbottom"] = 2.8194352;
    
    //Model_i.Init(Parameters);
    Model_i.Update(Parameters);
}

void EWSMLEP2testclass::sqrtsTEST() {
    double expected = 200.0;
    double result = sqrt_s;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::MwTEST() {
    double expected = 80.360848365211552;
    double result = Mw;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::GammaZTEST() {
    double expected = 2.494980275134754;
    double result = GammaZ;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::chi_Z_real() {
    double expected = 0.4714952963;
    double result = myLEP2->chi_Z(s,Mw,GammaZ).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::chi_Z_imag() {
    double expected = -0.0162861206;
    double result = myLEP2->chi_Z(s,Mw,GammaZ).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::G1_mu_noWeak() {
    double expected = 1.387618480024821;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double mfp = 0.0;
    double result = myLEP2->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::G2_mu_noWeak() {
    double expected = 1.162518560532168;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double mfp = 0.0;
    double result = myLEP2->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::G3_mu_noWeak() {
    double expected = 1.021139815000306;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double mfp = 0.0;
    double result = myLEP2->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, mfp, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::F_za_real() {
    double expected = 1.2063205379;
    double result = myLEP2->F_za_0(s, Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::F_za_imag() {
    double expected = 5.3999103856;
    double result = myLEP2->F_za_0(s, Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMLEP2testclass::B_WW_d_0_hat_real() {
    double t = -0.4, u = -0.3;
    double expected = myLEP2->B_WW_d_0_hat_TEST(s,t,u,Mw).real();
    double result = myLEP2->B_WW_d_0_hat(s,t,u,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::B_WW_d_0_hat_imag() {
    double t = -0.4, u = -0.3;
    double expected = myLEP2->B_WW_d_0_hat_TEST(s,t,u,Mw).imag();
    double result = myLEP2->B_WW_d_0_hat(s,t,u,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMLEP2testclass::rhoef_NonUnitary_MU_real() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);       
}

void EWSMLEP2testclass::rhoef_NonUnitary_MU_imag() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);           
}

void EWSMLEP2testclass::kappae_NonUnitary_MU_real() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappae_NonUnitary_MU_imag() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaf_NonUnitary_MU_real() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaf_NonUnitary_MU_imag() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaef_NonUnitary_MU_real() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaef_NonUnitary_MU_imag() {
    double t = -0.4;
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double expected = myLEP2_NU->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::rhoef_NonUnitary_UP_real() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);       
}

void EWSMLEP2testclass::rhoef_NonUnitary_UP_imag() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->rho_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);           
}

void EWSMLEP2testclass::kappae_NonUnitary_UP_real() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappae_NonUnitary_UP_imag() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->kappa_e(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaf_NonUnitary_UP_real() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaf_NonUnitary_UP_imag() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->kappa_f(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaef_NonUnitary_UP_real() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double result = myLEP2->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::kappaef_NonUnitary_UP_imag() {
    double t = -0.4;
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double expected = myLEP2_NU->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double result = myLEP2->kappa_ef(s,t,Mw,I3f,Qf,mf,0.0,true,true,false).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::Delta_rhoef_TOP_NonUnitary_real() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_rho_ef_TOP(s,t,u,Mw,true).real();
    double result = myLEP2->Delta_rho_ef_TOP(s,t,u,Mw,true).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);        
}

void EWSMLEP2testclass::Delta_rhoef_TOP_NonUnitary_imag() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_rho_ef_TOP(s,t,u,Mw,true).imag();
    double result = myLEP2->Delta_rho_ef_TOP(s,t,u,Mw,true).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);             
}

void EWSMLEP2testclass::Delta_kappae_TOP_NonUnitary_real() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_kappa_e_TOP(s,t,u,Mw,true).real();
    double result = myLEP2->Delta_kappa_e_TOP(s,t,u,Mw,true).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::Delta_kappae_TOP_NonUnitary_imag() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_kappa_e_TOP(s,t,u,Mw,true).imag();
    double result = myLEP2->Delta_kappa_e_TOP(s,t,u,Mw,true).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::Delta_kappaf_TOP_NonUnitary_real() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_kappa_f_TOP(s,t,u,Mw,true).real();
    double result = myLEP2->Delta_kappa_f_TOP(s,t,u,Mw,true).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::Delta_kappaf_TOP_NonUnitary_imag() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_kappa_f_TOP(s,t,u,Mw,true).imag();
    double result = myLEP2->Delta_kappa_f_TOP(s,t,u,Mw,true).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::Delta_kappaef_TOP_NonUnitary_real() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_kappa_ef_TOP(s,t,u,Mw,true).real();
    double result = myLEP2->Delta_kappa_ef_TOP(s,t,u,Mw,true).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::Delta_kappaef_TOP_NonUnitary_imag() {
    double mf = SM->getQuarks(SM->BOTTOM).getMass();
    double t = -0.4, u = -s - t + 2.0*mf*mf;
    double expected = myLEP2_NU->Delta_kappa_ef_TOP(s,t,u,Mw,true).imag();
    double result = myLEP2->Delta_kappa_ef_TOP(s,t,u,Mw,true).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::B_WW_0_real() {
    //double I3f = SM->getLeptons(SM->MU).getIsospin();
    //double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0, u = - s - t;
    double mu = Mw;
    double expected = -0.00010760320088880713; // ??????
    //double result = myLEP2->B_WW_d_0_hat(s, t, u, Mw).real();
    double result = myLEP2->B_WW_d_0(mu, s, t, u, Mw).real();
    //double result = myLEP2_NU->B_WW_d_0_hat(s, t, u, Mw).real(); //!! TEST !!
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMLEP2testclass::B_WW_0_imag() {
    //double I3f = SM->getLeptons(SM->MU).getIsospin();
    //double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0, u = - s - t;
    double mu = Mw;
    double expected = 0.00048542322208538363; // ??????
    //double result = myLEP2->B_WW_d_0_hat(s, t, u, Mw).imag();
    double result = myLEP2->B_WW_d_0(mu, s, t, u, Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);        
}

void EWSMLEP2testclass::B_ZZ_0_real() {
    double t = - s/2.0, u = - s - t;
    double mu = Mw;
    double expected = 0.00036346250937252961;
    double result = myLEP2->B_ZZ_0(mu, s, t, u).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMLEP2testclass::B_ZZ_0_imag() {
    double t = - s/2.0, u = - s - t;
    double mu = Mw;
    double expected = 0.00040616916658261799;
    double result = myLEP2->B_ZZ_0(mu, s, t, u).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);        
}

void EWSMLEP2testclass::Delta_rho_ef_WW_real() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double t = - s/2.0, u = - s - t;
    //double mu = Mw;
    double expected = -0.0068839947; // ??????
    //double result = myLEP2->Delta_rho_ef_WW_hat(s,t,u,Mw,I3f).real();
    double result = myLEP2_NU->Delta_rho_ef_WW_hat(s,t,u,Mw,I3f).real(); //!! TEST !!
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);     
}

void EWSMLEP2testclass::Delta_rho_ef_WW_imag() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double t = - s/2.0, u = - s - t;
    //double mu = Mw;
    double expected = 0.0310553110; // ??????
    //double result = myLEP2->Delta_rho_ef_WW_hat(s,t,u,Mw,I3f).imag();
    double result = myLEP2_NU->Delta_rho_ef_WW_hat(s,t,u,Mw,I3f).imag(); //!! TEST !!
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMLEP2testclass::Delta_rho_ef_ZZ_real() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double t = - s/2.0, u = - s - t;
    double mu = Mw;
    double expected = 0.0024644708;
    double result = myLEP2->Delta_rho_ef_ZZ(mu,s,t,u,Mw,I3f,Qf).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);      
}

void EWSMLEP2testclass::Delta_rho_ef_ZZ_imag() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double t = - s/2.0, u = - s - t;
    double mu = Mw;
    double expected = 0.0027540448;
    double result = myLEP2->Delta_rho_ef_ZZ(mu,s,t,u,Mw,I3f,Qf).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);     
}

void EWSMLEP2testclass::Delta_rho_ef_WW_charm_real() {
    double I3f = SM->getQuarks(SM->CHARM).getIsospin();
    double snew = 205.0*205.0;
    double t = - snew/2.0, u = - snew - t;
    //double mu = Mw;
    double expected = -0.0190989409; // ??????
    //double result = myLEP2->Delta_rho_ef_WW_hat(snew,t,u,Mw,I3f).real();
    double result = myLEP2_NU->Delta_rho_ef_WW_hat(snew,t,u,Mw,I3f).real(); //!! TEST !!
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);     
}

void EWSMLEP2testclass::Delta_rho_ef_WW_charm_imag() {
    double I3f = SM->getQuarks(SM->CHARM).getIsospin();
    double snew = 205.0*205.0;
    double t = - snew/2.0, u = - snew - t;
    //double mu = Mw;
    double expected = 0.0028803190; // ??????
    //double result = myLEP2->Delta_rho_ef_WW_hat(snew,t,u,Mw,I3f).imag();
    double result = myLEP2_NU->Delta_rho_ef_WW_hat(snew,t,u,Mw,I3f).imag(); //!! TEST !!
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);     
}

void EWSMLEP2testclass::Delta_rho_ef_ZZ_charm_real() {
    double I3f = SM->getQuarks(SM->CHARM).getIsospin();
    double Qf = SM->getQuarks(SM->CHARM).getCharge();
    double snew = 205.0*205.0;
    double t = - snew/2.0, u = - snew - t;
    double mu = Mw;
    double expected = -0.0024895438;
    double result = myLEP2->Delta_rho_ef_ZZ(mu,snew,t,u,Mw,I3f,Qf).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);     
}

void EWSMLEP2testclass::Delta_rho_ef_ZZ_charm_imag() {
    double I3f = SM->getQuarks(SM->CHARM).getIsospin();
    double Qf = SM->getQuarks(SM->CHARM).getCharge();
    double snew = 205.0*205.0;
    double t = - snew/2.0, u = - snew - t;
    double mu = Mw;
    double expected = -0.0033604175;
    double result = myLEP2->Delta_rho_ef_ZZ(mu,snew,t,u,Mw,I3f,Qf).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);   
}

void EWSMLEP2testclass::G1_mu_WWBox_TEST() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_1(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true, false);
    double result = myLEP2->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_1_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, true, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G2_mu_WWBox_TEST() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_2(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true, false);
    double result = myLEP2->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_2_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, true, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G3_mu_WWBox_TEST() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_3(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true, false);
    double result = myLEP2->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_3_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, true, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}
void EWSMLEP2testclass::G1_UP_WWBox_TEST() {
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_1(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true, false);
    double result = myLEP2->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_1_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, true, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G2_UP_WWBox_TEST() {
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_2(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true, false);
    double result = myLEP2->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_2_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, true, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G3_UP_WWBox_TEST() {
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_3(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true, false);
    double result = myLEP2->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_3_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, true, false);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G1_mu_ZZBox_TEST() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_1(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, false, true);
    double result = myLEP2->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_1_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G2_mu_ZZBox_TEST() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_2(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, false, true);
    double result = myLEP2->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_2_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G3_mu_ZZBox_TEST() {
    double I3f = SM->getLeptons(SM->MU).getIsospin();
    double Qf = SM->getLeptons(SM->MU).getCharge();
    double mf = SM->getLeptons(SM->MU).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_3(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, false, true);
    double result = myLEP2->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_3_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}
void EWSMLEP2testclass::G1_UP_ZZBox_TEST() {
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_1(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, false, true);
    double result = myLEP2->G_1_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_1_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G2_UP_ZZBox_TEST() {
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_2(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, false, true);
    double result = myLEP2->G_2_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_2_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}

void EWSMLEP2testclass::G3_UP_ZZBox_TEST() {
    double I3f = SM->getQuarks(SM->UP).getIsospin();
    double Qf = SM->getQuarks(SM->UP).getCharge();
    double mf = SM->getQuarks(SM->UP).getMass();
    double t = - s/2.0;
    double expected = myLEP2->G_3(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, false, true);
    double result = myLEP2->G_3_noBox(s, Mw, GammaZ, I3f, Qf, mf, 0.0, false)
                    + myLEP2->G_3_box(s, t, Mw, GammaZ, I3f, Qf, mf, 0.0, false, true);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);
}







