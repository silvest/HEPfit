/*
 * File:   EWSMEW1testclass.cpp
 * Author: mishima
 */

#include "EWSMEW1testclass.h"
#include "EWSMsetParameters.h"

CPPUNIT_TEST_SUITE_REGISTRATION(EWSMEW1testclass);

EWSMEW1testclass::EWSMEW1testclass() {
}

EWSMEW1testclass::~EWSMEW1testclass() {
}

void EWSMEW1testclass::setUp() {
    mySM = new StandardModel(true);
    mySM->InitializeModel();
    setSMparameters(*mySM);   
    myCache = new EWSMcache(*mySM);
    myEW1 = new EWSMOneLoopEW(*myCache);

    Mw = myCache->Mw(mySM->Mw_tree());/* Tests are done with the tree-level Mw */
    Mw2 = Mw*Mw;
    Mz = myCache->Mz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = myCache->Mt();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-10; 
    epsilon_Li2 = 1.0e-7; /* for quantities with the dilogarithm */

}

void EWSMEW1testclass::tearDown() {
    delete myEW1;        
    delete myCache;
    delete mySM;    
}

void EWSMEW1testclass::DeltaAlpha_l() {
    double ZFITTER = 0.031418982830747; /* ZFITTER result*/
    double result = myEW1->DeltaAlpha_l(Mz2);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void EWSMEW1testclass::SigmaWW_bos_Mw_0() {
    double W0 = -1.702845692408391; /* ZFITTER result*/
    double result = myEW1->SigmaWW_bos(Mw,0.0,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0*Mw2, result, delta);
}

void EWSMEW1testclass::SigmaWW_fer_Mw_0() {
    double W0F = 7.306899137576438; /* ZFITTER result*/
    double result = myEW1->SigmaWW_fer(Mw,0.0,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0F*Mw2, result, delta);
}

void EWSMEW1testclass::SigmaWW_bos_Mw_Mw2_real() {
    double XMM1 = -4.212516913677337; /* ZFITTER result*/
    double result = myEW1->SigmaWW_bos(Mw,Mw2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1*Mw2, result, delta);
}

void EWSMEW1testclass::SigmaWW_bos_Mw_Mw2_imag() {
    double XMM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myEW1->SigmaWW_bos(Mw,Mw2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1*Mw2, result, delta);
}

void EWSMEW1testclass::SigmaWW_fer_Mw_Mw2_real() {
    double XMM1F = 11.188122745796212; /* ZFITTER result*/
    double result = myEW1->SigmaWW_fer(Mw,Mw2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1F*Mw2, result, delta);   
}

void EWSMEW1testclass::SigmaWW_fer_Mw_Mw2_imag() {
    double XMM1F = 9.422358620516286; /* ZFITTER result*/
    double result = myEW1->SigmaWW_fer(Mw,Mw2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1F*Mw2, result, delta); 
}

void EWSMEW1testclass::SigmaZZ_bos_Mw_Mz2_real() {
    double XZM1 = -3.320247142827423; /* ZFITTER result*/
    double result = myEW1->SigmaZZ_bos(Mw,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1*Mw2, result, delta);  
}

void EWSMEW1testclass::SigmaZZ_bos_Mw_Mz2_imag() {
    double XZM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myEW1->SigmaZZ_bos(Mw,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1*Mw2, result, delta);     
}

void EWSMEW1testclass::SigmaZZ_fer_Mw_Mz2_real() {
    double XZM1F = 14.431690239155724; /* ZFITTER result*/
    double result = myEW1->SigmaZZ_fer(Mw,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1F*Mw2, result, delta);  
}

void EWSMEW1testclass::SigmaZZ_fer_Mw_Mz2_imag() {
    double XZM1F = 9.893711064659454; /* ZFITTER result*/
    double result = myEW1->SigmaZZ_fer(Mw,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1F*Mw2, result, delta); 
}

void EWSMEW1testclass::PiZgamma_bos_Mw_Mz2_real() {
    double XAMM1 = 1.999074211569152; /* ZFITTER result*/
    double result = myEW1->PiZgamma_bos(Mw,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1, result, delta);    
}

void EWSMEW1testclass::PiZgamma_bos_Mw_Mz2_imag() {
    double XAMM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myEW1->PiZgamma_bos(Mw,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1, result, delta); 
}

void EWSMEW1testclass::PiZgamma_fer_Mw_Mz2_real() {
    double XAMM1F = - 1.641562584085382; /* ZFITTER result*/
    double result = myEW1->PiZgamma_fer(Mw,Mz2,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1F, result, delta);    
}

void EWSMEW1testclass::PiZgamma_fer_Mw_Mz2_imag() {
    double XAMM1F = - 4.547472366444150; /* ZFITTER result*/
    double result = myEW1->PiZgamma_fer(Mw,Mz2,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1F, result, delta); 
}

void EWSMEW1testclass::PiGammaGamma_fer_Mw_0() {
    double ZFITTER = 110.687763755483232; /* ZFITTER result*/
    double result = myEW1->PiGammaGamma_fer(Mw,0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::B0_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B0(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->B0(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void EWSMEW1testclass::B1_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B1(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->B1(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = - log(cW2)/2.0;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void EWSMEW1testclass::B21_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B21(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->B21(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = log(cW2)/3.0;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void EWSMEW1testclass::Bf_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->Bf(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->Bf(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = - log(cW2)/3.0;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void EWSMEW1testclass::B0p_diff_Mw2_0_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B0p(Mw, Mw2, 0.0, Mw).real();
    double result_Mz = myPV->B0p(Mz, Mw2, 0.0, Mw).real();
    double MZtoMW = - log(cW2)/2.0/Mw2;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void EWSMEW1testclass::SigmaWW_bos_diff_0_real() {
    double result_Mw = myEW1->SigmaWW_bos(Mw,0.0,Mw).real();
    double result_Mz = myEW1->SigmaWW_bos(Mz,0.0,Mw).real();
    double MZtoMW = Mw2*(3.0/4.0/cW2 - 3.0/2.0 + 3.0/4.0*cW2 + 3.0/4.0*cW2*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::SigmaWW_fer_diff_0_real() {
    double result_Mw = myEW1->SigmaWW_fer(Mw,0.0,Mw).real();
    double result_Mz = myEW1->SigmaWW_fer(Mz,0.0,Mw).real();
    double MZtoMW = 0.0;
    for (int i=0; i<6; i++) {
        MZtoMW += - 1.0/2.0*pow(myCache->ml((StandardModel::lepton)i),2.0)*log(cW2);
        MZtoMW += - 3.0/2.0*pow(myCache->mq((StandardModel::quark)i,Mz),2.0)*log(cW2);
    }
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::SigmaWW_bos_diff_Mw2_real() {
    double result_Mw = myEW1->SigmaWW_bos(Mw,Mw2,Mw).real();
    double result_Mz = myEW1->SigmaWW_bos(Mz,Mw2,Mw).real();
    double MZtoMW = Mw2*(3.0/4.0/cW2 - 17.0/3.0)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::SigmaWW_fer_diff_Mw2_real() {
    double result_Mw = myEW1->SigmaWW_fer(Mw,Mw2,Mw).real();
    double result_Mz = myEW1->SigmaWW_fer(Mz,Mw2,Mw).real();
    double MZtoMW = 24.0/6.0*Mw2*log(cW2);
    for (int i=0; i<6; i++) {
        MZtoMW += - 1.0/2.0*pow(myCache->ml((StandardModel::lepton)i),2.0)*log(cW2);
        MZtoMW += - 3.0/2.0*pow(myCache->mq((StandardModel::quark)i,Mz),2.0)*log(cW2);
    }
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::SigmaZZ_bos_diff_Mz2_real() {
    double result_Mw = myEW1->SigmaZZ_bos(Mw,Mz2,Mw).real();
    double result_Mz = myEW1->SigmaZZ_bos(Mz,Mz2,Mw).real();
    double MZtoMW = Mw2*(11.0/12.0/cW2 + 7.0/6.0 - 7.0*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::SigmaZZ_fer_diff_Mz2_real() {
    double result_Mw = myEW1->SigmaZZ_fer(Mw,Mz2,Mw).real();
    double result_Mz = myEW1->SigmaZZ_fer(Mz,Mz2,Mw).real();
    double MZtoMW = Mw2*(24.0/6.0/cW2 - 2.0*sW2/3.0/cW2*12.0 
                         + 4.0*sW2*sW2/3.0/cW2*8.0)*log(cW2);
    for (int i=0; i<6; i++) {
        MZtoMW += - 1.0/2.0*pow(myCache->ml((StandardModel::lepton)i),2.0)*log(cW2);
        MZtoMW += - 3.0/2.0*pow(myCache->mq((StandardModel::quark)i,Mz),2.0)*log(cW2);
    }
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::PiGammaGamma_fer_diff_0_real() {
    double result_Mw = myEW1->PiGammaGamma_fer(Mw,0.0).real();
    double result_Mz = myEW1->PiGammaGamma_fer(Mz,0.0).real();
    double MZtoMW = 4.0/3.0*8.0*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::PiGammaGamma_bos_diff_Mz2_real() {
    double result_Mw = myEW1->PiGammaGamma_bos(Mw,Mz2,Mw).real();
    double result_Mz = myEW1->PiGammaGamma_bos(Mz,Mz2,Mw).real();
    double MZtoMW = 1.0/cW2*(1.0/12.0/cW2 + 7.0/6.0 - 7.0*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::PiGammaGamma_fer_diff_Mz2_real() {
    double result_Mw = myEW1->PiGammaGamma_fer(Mw,Mz2).real();
    double result_Mz = myEW1->PiGammaGamma_fer(Mz,Mz2).real();
    double MZtoMW = 4.0/3.0*8.0*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::PiZgamma_bos_diff_Mz2_real() {
    double result_Mw = myEW1->PiZgamma_bos(Mw,Mz2,Mw).real();
    double result_Mz = myEW1->PiZgamma_bos(Mz,Mz2,Mw).real();
    double MZtoMW = (1.0/12.0/cW2 + 7.0/6.0 - 7.0*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::PiZgamma_fer_diff_Mz2_real() {
    double result_Mw = myEW1->PiZgamma_fer(Mw,Mz2,Mw).real();
    double result_Mz = myEW1->PiZgamma_fer(Mz,Mz2,Mw).real();
    double MZtoMW = (12.0/3.0 - 4.0/3.0*sW2*8.0)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void EWSMEW1testclass::DeltaRhobar_bos_Mw() {
    double ZFITTER = -0.892269770849914; /* ZFITTER result*/    
    double result = (myEW1->SigmaWW_bos(Mw,Mw2,Mw).real()
                     - myEW1->SigmaZZ_bos(Mw,Mz2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void EWSMEW1testclass::DeltaRhobar_fer_Mw() {
    double ZFITTER = -3.243567493359512; /* ZFITTER result*/    
    double result = (myEW1->SigmaWW_fer(Mw,Mw2,Mw).real()
                     - myEW1->SigmaZZ_fer(Mw,Mz2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::DeltaRhobarW_bos_Mw() {
    double ZFITTER = 2.509671221268946; /* ZFITTER result*/    
    double result = (myEW1->SigmaWW_bos(Mw,0.0,Mw).real()
                     - myEW1->SigmaWW_bos(Mw,Mw2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void EWSMEW1testclass::DeltaRhobarW_fer_Mw() {
    double ZFITTER = -3.881223608219774; /* ZFITTER result*/    
    double result = (myEW1->SigmaWW_fer(Mw,0.0,Mw).real()
                     - myEW1->SigmaWW_fer(Mw,Mw2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}
    
void EWSMEW1testclass::DeltaRhobar_bos_Mz() {
    double ZFITTER = -1.257094198330098; /* ZFITTER result*/    
    double result = (myEW1->SigmaWW_bos(Mz,Mw2,Mw).real()
                     - myEW1->SigmaZZ_bos(Mz,Mz2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void EWSMEW1testclass::DeltaRhobar_fer_Mz() {
    double ZFITTER = -3.132030503808709; /* ZFITTER result*/    
    double result = (myEW1->SigmaWW_fer(Mz,Mw2,Mw).real()
                     - myEW1->SigmaZZ_fer(Mz,Mz2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::TEST_DeltaRhobar_bos_Mw() {
    double test_result = myEW1->TEST_DeltaRhobar_bos(Mw);
    double result = (myEW1->SigmaWW_bos(Mw,Mw2,Mw).real()
                     - myEW1->SigmaZZ_bos(Mw,Mz2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_result, result, delta);    
}

void EWSMEW1testclass::TEST_DeltaRhobarW_bos_Mw() {
    double test_result = myEW1->TEST_DeltaRhobarW_bos(Mw);
    double result = (myEW1->SigmaWW_bos(Mw,0.0,Mw).real()
                     - myEW1->SigmaWW_bos(Mw,Mw2,Mw).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_result, result, delta);    
}

void EWSMEW1testclass::SCALE() {
    double SCALE = -0.002574576533187; /* ZFITTER result*/    
    double result = - mySM->getAle()/4.0/M_PI/sW2
                      * (myEW1->SigmaWW_bos(Mw,Mw2,Mw).real()
                         + myEW1->SigmaWW_fer(Mw,Mw2,Mw).real()
                         - myEW1->SigmaZZ_bos(Mw,Mz2,Mw).real()
                         - myEW1->SigmaZZ_fer(Mw,Mz2,Mw).real())/Mw2
                    - myEW1->DeltaRho(Mw);
    result *= cW2/sW2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(SCALE, result, delta);         
}
 
void EWSMEW1testclass::DeltaRho() {
    double ZFITTER = 0.012013822484982; /* ZFITTER result*/  
    double result = myEW1->DeltaRho(Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
} 

void EWSMEW1testclass::DeltaR_rem() {
    double DRREMN = 0.011412292032499; /* ZFITTER result*/   
    double result = myEW1->DeltaR_rem(Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DRREMN, result, delta);
}
    
void EWSMEW1testclass::SigmaPrimeWW_bos_Mw_Mw2_real() {
    double XWFM1 = 2.066448743517876; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_WW_bos_Mw2(Mw,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1, result, delta);
}

void EWSMEW1testclass::SigmaPrimeWW_bos_Mw_Mw2_imag() {
    double XWFM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_WW_bos_Mw2(Mw,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1, result, delta);
}

void EWSMEW1testclass::SigmaPrimeWW_fer_Mw_Mw2_real() {
    double XWFM1F = -0.966097280796035; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_WW_fer_Mw2(Mw,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F, result, delta);
}

void EWSMEW1testclass::SigmaPrimeWW_fer_Mw_Mw2_imag() {
    double XWFM1F = -9.424777230871838; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_WW_fer_Mw2(Mw,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F, result, delta);
}

void EWSMEW1testclass::SigmaPrimeZZ_bos_Mw_Mz2_real() {
    double XWFM1 = 3.109040605851439; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_ZZ_bos_Mz2(Mw,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1*cW2, result, delta);
}

void EWSMEW1testclass::SigmaPrimeZZ_bos_Mw_Mz2_imag() {
    double XWFM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_ZZ_bos_Mz2(Mw,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1*cW2, result, delta);
}

void EWSMEW1testclass::SigmaPrimeZZ_fer_Mw_Mz2_real() {
    double XWFM1F = -0.482507863333779; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_ZZ_fer_Mz2(Mw,Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F*cW2, result, delta);
}

void EWSMEW1testclass::SigmaPrimeZZ_fer_Mw_Mz2_imag() {
    double XWFM1F = -9.911978176545077; /* ZFITTER result*/
    double result = myEW1->SigmaPrime_ZZ_fer_Mz2(Mw,Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F*cW2, result, delta);
}

void EWSMEW1testclass::C0_Mz2_Mt_Mw_Mt_real() {
    double XS3T = 0.000024702795478; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mt,Mw,Mt).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T, result, delta);    
}

void EWSMEW1testclass::C0_Mz2_Mt_Mw_Mt_imag() {
    double XS3T = 0.000000000000000; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mt,Mw,Mt).imag();
    delete myPV;
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T, result, delta);    
}

void EWSMEW1testclass::C0_Mz2_0_Mw_0_real() {
    double XS3T0 = 0.000097066726252; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,0.0,Mw,0.0).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T0, result, delta);    
}

void EWSMEW1testclass::C0_Mz2_0_Mw_0_imag() {
    double XS3T0 = 0.000309606030132; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,0.0,Mw,0.0).imag();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T0, result, delta); 
}

void EWSMEW1testclass::C0_Mz2_Mw_Mt_Mw_real() {
    double XS3W = 0.000044262165905; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,Mt,Mw).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W, result, delta);    
}

void EWSMEW1testclass::C0_Mz2_Mw_Mt_Mw_imag() {
    double XS3W = 0.000000000000000; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,Mt,Mw).imag();
    delete myPV;
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W, result, delta);    
}

void EWSMEW1testclass::C0_Mz2_Mw_0_Mw_real() {
    double XS3W0 = 0.000172249280944; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,0.0,Mw).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W0, result, delta);    
}

void EWSMEW1testclass::C0_Mz2_Mw_0_Mw_imag() {
    double XS3W0 = 0.000000000000000; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,0.0,Mw).imag();
    delete myPV;
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W0, result, delta);    
}

void EWSMEW1testclass::FZa_0_real() {
    double V1ZZ = 1.079736267392905; /* ZFITTER result*/
    double result = myEW1->FZa_0(Mz*Mz,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1ZZ, result, delta);    
}

void EWSMEW1testclass::FZa_0_imag() {
    double V1ZIM = 1.712725454479850; /* ZFITTER result*/
    double result = myEW1->FZa_0(Mz*Mz,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1ZIM, result, delta);    
}
 
void EWSMEW1testclass::FWa_0_real() {
    double V1ZW = 1.175174895095040; /* ZFITTER result*/
    double result = myEW1->FWa_0(Mz*Mz,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1ZW, result, delta);    
}

void EWSMEW1testclass::FWa_0_imag() {
    double V1WIM = 2.082775415576450; /* ZFITTER result*/
    double result = myEW1->FWa_0(Mz*Mz,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1WIM, result, delta);    
}

void EWSMEW1testclass::FWn_0_real() {
    double V2ZWW = -1.359807910427424; /* ZFITTER result*/
    double result = myEW1->FWn_0(Mz*Mz,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V2ZWW, result, delta);    
}

void EWSMEW1testclass::FWn_0_imag() {
    double ZFITTER = 0.000000000000000; /* ZFITTER result*/
    double result = myEW1->FWn_0(Mz*Mz,Mw).imag();
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void EWSMEW1testclass::FWa_t_real() {
    double WWv11 = -0.695682293850053; /* ZFITTER result*/
    double result = myEW1->FWa_t(Mz*Mz,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv11, result, delta);    
}

void EWSMEW1testclass::FbarWa_t_real() {
    double WWv12 = 3.496220092867632; /* ZFITTER result*/
    double result = myEW1->FbarWa_t(Mz*Mz,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv12, result, delta);    
}

void EWSMEW1testclass::FWn_t_real() {
    double WWv2 = -1.460637714980599; /* ZFITTER result*/
    double result = myEW1->FWn_t(Mz*Mz,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv2, result, delta);      
}

void EWSMEW1testclass::TEST_FWn_0_real() {
    double V2ZWW = -1.359807910427424; /* ZFITTER result*/
    double result = myEW1->TEST_FWn(Mz*Mz, 0.0, Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V2ZWW, result, delta);    
}

void EWSMEW1testclass::TEST_FWn_t_real() {
    double WWv2 = -1.460637714980599; /* ZFITTER result*/
    double result = myEW1->TEST_FWn(Mz*Mz, myCache->Mt(), Mw).real()
                    - myEW1->TEST_FWn(Mz*Mz, 0.0, Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv2, result, delta);    
}

void EWSMEW1testclass::FW_diff_bb_dd_real() {
    double VTB = -2.649420523952066; /* ZFITTER result*/
    double result = myEW1->FW_q(Mz*Mz, mySM->BOTTOM, Mw).real()
                    - myEW1->FW_q(Mz*Mz, mySM->DOWN, Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(VTB, result, delta);    
}

void EWSMEW1testclass::f_convert() {
    double RENORM = 1.000000000000000; /* ZFITTER result*/    
    double result = myCache->f_AlphaToGF(Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(RENORM, result, delta);  
}

void EWSMEW1testclass::Xt_GF() {
    double TOPX2 = 0.003198951884234; /* ZFITTER result*/
    double result = myCache->Xt_GF();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(TOPX2, result, delta);    
}

void EWSMEW1testclass::DeltaRho_CORRHO() {
    double CORRHO = -0.012013822484982; /* ZFITTER result*/  
    double result = myEW1->DeltaRho(Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-CORRHO, result, delta);
} 

void EWSMEW1testclass::DeltaRho_G() {
    double DROBLO = 0.012013822484982; /* ZFITTER result*/  
    double result = myEW1->DeltaRho(Mw)*myCache->f_AlphaToGF(Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DROBLO, result, delta);
} 

void EWSMEW1testclass::DeltaRbar_rem() {
    double DRREMD = 0.004147533855948; /* ZFITTER result*/   
    double result = myEW1->DeltaRbar_rem(Mw);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DRREMD, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_NEUTRINO_1_real() {
    double ZFITTER = -0.002782064884822; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->NEUTRINO_1,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_NEUTRINO_2_real() {
    double ZFITTER = -0.002782064884822; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->NEUTRINO_2,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_NEUTRINO_3_real() {
    double ZFITTER = -0.002782064884822; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->NEUTRINO_3,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_ELECTRON_real() {
    double ZFITTER = -0.005521418517950; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->ELECTRON,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_MU_real() {
    double ZFITTER = -0.005521418517950; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->MU,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_TAU_real() {
    double ZFITTER = -0.005521418517950; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->TAU,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_UP_real() {
    double ZFITTER = -0.004833425293458; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->UP,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_DOWN_real() {
    double ZFITTER = -0.003920307415749; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->DOWN,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_CHARM_real() {
    double ZFITTER = -0.004833425293458; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->CHARM,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_STRANGE_real() {
    double ZFITTER = -0.003920307415749; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->STRANGE,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_NEUTRINO_1_imag() {
    double ZFITTER = -0.000306748406150; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->NEUTRINO_1,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_NEUTRINO_2_imag() {
    double ZFITTER = -0.000306748406150; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->NEUTRINO_2,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_NEUTRINO_3_imag() {
    double ZFITTER = -0.000306748406150; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->NEUTRINO_3,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_ELECTRON_imag() {
    double ZFITTER = -0.004905987392756; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->ELECTRON,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_MU_imag() {
    double ZFITTER = -0.004905987392756; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->MU,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_TAU_imag() {
    double ZFITTER = -0.004905987392756; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_l(mySM->TAU,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_UP_imag() {
    double ZFITTER = -0.003730010419090; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->UP,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_DOWN_imag() {
    double ZFITTER = -0.002196930756888; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->DOWN,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_CHARM_imag() {
    double ZFITTER = -0.003730010419090; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->CHARM,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_STRANGE_imag() {
    double ZFITTER = -0.002196930756888; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->STRANGE,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_NEUTRINO_1_real() {
    double ZFITTER = -0.000632677138537; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->NEUTRINO_1,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_NEUTRINO_2_real() {
    double ZFITTER = -0.000632677138537; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->NEUTRINO_2,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_NEUTRINO_3_real() {
    double ZFITTER = -0.000632677138537; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->NEUTRINO_3,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_ELECTRON_real() {
    double ZFITTER = 0.000905843167939; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->ELECTRON,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_MU_real() {
    double ZFITTER = 0.000905843167939; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->MU,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_TAU_real() {
    double ZFITTER = 0.000905843167939; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->TAU,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_UP_real() {
    double ZFITTER = 0.000468044616853; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->UP,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_DOWN_real() {
    double ZFITTER = -0.000044795485306; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->DOWN,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_CHARM_real() {
    double ZFITTER = 0.000468044616853; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->CHARM,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_STRANGE_real() {
    double ZFITTER = -0.000044795485306; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->STRANGE,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_NEUTRINO_1_imag() {
    double ZFITTER = 0.012600622680476; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->NEUTRINO_1,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_NEUTRINO_2_imag() {
    double ZFITTER = 0.012600622680476; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->NEUTRINO_2,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_NEUTRINO_3_imag() {
    double ZFITTER = 0.012600622680476; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->NEUTRINO_3,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_ELECTRON_imag() {
    double ZFITTER = 0.015168069190180; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->ELECTRON,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_MU_imag() {
    double ZFITTER = 0.015168069190180; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->MU,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_TAU_imag() {
    double ZFITTER = 0.015168069190180; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_l(mySM->TAU,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_UP_imag() {
    double ZFITTER = 0.014431287916457; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->UP,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_DOWN_imag() {
    double ZFITTER = 0.013575472413222; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->DOWN,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_CHARM_imag() {
    double ZFITTER = 0.014431287916457; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->CHARM,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_STRANGE_imag() {
    double ZFITTER = 0.013575472413222; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->STRANGE,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_BOTTOM_real() {
    double ZFITTER = -0.018424186901021; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->BOTTOM,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaRho_rem_BOTTOM_imag() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myEW1->deltaRho_rem_q(mySM->BOTTOM,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_BOTTOM_real() {
    double ZFITTER = 0.007207144257329; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->BOTTOM,Mw).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::deltaKappa_rem_BOTTOM_imag() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myEW1->deltaKappa_rem_q(mySM->BOTTOM,Mw).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void EWSMEW1testclass::rho_GammaW_leptons() {
    double ROW = 0.996736176315242; /* ZFITTER result*/
    double result = myEW1->rho_GammaW_l(mySM->NEUTRINO_1, mySM->ELECTRON,Mw);
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ROW, result, delta);   
}

void EWSMEW1testclass::rho_GammaW_quarks() {
    double ROW = 0.996281153730537; /* ZFITTER result*/
    double result = myEW1->rho_GammaW_q(mySM->UP, mySM->DOWN,Mw);
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ROW, result, delta);   
}



