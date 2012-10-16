/*
 * File:   EWSMEW1HVtestclass.cpp
 * Author: mishima
 */

#include "EWSMEW1HVtestclass.h"
#include "EWSMsetParameters.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMEW1HVtestclass);

EWSMEW1HVtestclass::EWSMEW1HVtestclass() {
}

EWSMEW1HVtestclass::~EWSMEW1HVtestclass() {
}

void EWSMEW1HVtestclass::setUp() {
    mySM = new StandardModel(true);
    mySM->InitializeModel();
    setSMparameters(*mySM);   
    myEW1HV = new EWSMOneLoopEW_HV(*mySM);

    Mw = mySM->Mw_tree();/* Tests are done with the tree-level Mw */
    Mw2 = Mw*Mw;
    Mz = mySM->getMz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = mySM->getMtpole();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 
}

void EWSMEW1HVtestclass::tearDown() {
    delete myEW1HV; 
    delete mySM;  
}

void EWSMEW1HVtestclass::F_Hollik_0_Mw_Mw_real() {
    double expected = 0.0;
    double result = myEW1HV->F_Hollik(0.0, Mw, Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMEW1HVtestclass::F_Hollik_0_Mw_Mw_imag() {
    double expected = 0.0;
    double result = myEW1HV->F_Hollik(0.0, Mw, Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);      
}

void EWSMEW1HVtestclass::Fprime_Hollik_0_Mw_Mw_real() {
    double expected = 1.0/6.0/Mw2; 
    double result = myEW1HV->Fprime_Hollik(Mw, 0.0, Mw, Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMEW1HVtestclass::Fprime_Hollik_0_Mw_Mw_imag() {
    double expected = 0.0; 
    double result = myEW1HV->Fprime_Hollik(Mw, 0.0, Mw, Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);  
}

void EWSMEW1HVtestclass::SigmaWW_bos_Mz_Mz2_real() {
    double expected = myEW1HV->SigmaWW_bos(Mz,Mz2,Mw).real();
    double result = myEW1HV->SigmaWW_bos_Hollik(Mz,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaWW_bos_Mz_Mz2_imag() {
    double expected = myEW1HV->SigmaWW_bos(Mz,Mz2,Mw).imag();
    double result = myEW1HV->SigmaWW_bos_Hollik(Mz,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaZZ_bos_Mz_Mz2_real() {
    double expected = myEW1HV->SigmaZZ_bos(Mz,Mz2,Mw).real();
    double result = myEW1HV->SigmaZZ_bos_Hollik(Mz,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaZZ_bos_Mz_Mz2_imag() {
    double expected = myEW1HV->SigmaZZ_bos(Mz,Mz2,Mw).imag();
    double result = myEW1HV->SigmaZZ_bos_Hollik(Mz,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaZgamma_bos_Mz_Mz2_real() {
    double expected = myEW1HV->SigmaZgamma_bos(Mz,Mz2,Mw).real();
    double result = myEW1HV->SigmaZgamma_bos_Hollik(Mz,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaZgamma_bos_Mz_Mz2_imag() {
    double expected = myEW1HV->SigmaZgamma_bos(Mz,Mz2,Mw).imag();
    double result = myEW1HV->SigmaZgamma_bos_Hollik(Mz,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta); 
}

void EWSMEW1HVtestclass::SigmaZgamma_bos_Mz_0_real() {
    double expected = myEW1HV->SigmaZgamma_bos(Mz,0.0,Mw).real();
    double result = myEW1HV->SigmaZgamma_bos_Hollik(Mz,0.0,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaZgamma_bos_Mz_0_imag() {
    double expected = myEW1HV->SigmaZgamma_bos(Mz,0.0,Mw).imag();
    double result = myEW1HV->SigmaZgamma_bos_Hollik(Mz,0.0,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta); 
}

void EWSMEW1HVtestclass::SigmaGammaGamma_bos_Mz_Mz2_real() {
    double expected = myEW1HV->SigmaGammaGamma_bos(Mz,Mz2,Mw).real();
    double result = myEW1HV->SigmaGammaGamma_bos_Hollik(Mz,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::SigmaGammaGamma_bos_Mz_Mz2_imag() {
    double expected = myEW1HV->SigmaGammaGamma_bos(Mz,Mz2,Mw).imag();
    double result = myEW1HV->SigmaGammaGamma_bos_Hollik(Mz,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::PiGammaGamma_bos_Mz_Mz2_real() {
    double expected = myEW1HV->PiGammaGamma_bos(Mz,Mz2,Mw).real();
    double result = myEW1HV->PiGammaGamma_bos_Hollik(Mz,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::PiGammaGamma_bos_Mz_Mz2_imag() {
    double expected = myEW1HV->PiGammaGamma_bos(Mz,Mz2,Mw).imag();
    double result = myEW1HV->PiGammaGamma_bos_Hollik(Mz,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::PiGammaGamma_bos_Mz_0_real() {
    double expected = myEW1HV->PiGammaGamma_bos(Mz,0.0,Mw).real();
    double result = myEW1HV->PiGammaGamma_bos_Hollik(Mz,0.0,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}

void EWSMEW1HVtestclass::PiGammaGamma_bos_Mz_0_imag() {
    double expected = myEW1HV->PiGammaGamma_bos(Mz,0.0,Mw).imag();
    double result = myEW1HV->PiGammaGamma_bos_Hollik(Mz,0.0,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, result, delta);    
}





