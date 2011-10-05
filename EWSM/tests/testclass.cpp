/*
 * File:   testclass.cpp
 * Author: mishima
 */

#include "testclass.h"


CPPUNIT_TEST_SUITE_REGISTRATION(testclass);

testclass::testclass() {
}

testclass::~testclass() {
}

void testclass::setUp() {
    mySM = new StandardModel();
    testclass::setSMparameters(*mySM);   
    
    myEWSMC = new EWSMcommon(*mySM);
    myEWSMC->SetConstants();
    myEWSMC->Compute(mySM->Mw_tree());
    
    myOLEW = new OneLoopEW(*myEWSMC);

    Mw = myEWSMC->GetMw();
    Mw2 = Mw*Mw;
    Mz = myEWSMC->GetSM().getMz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = myEWSMC->GetSM().getQuarks(StandardModel::TOP).getMass();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-10; 
    epsilon_Li2 = 1.0e-7; /* for quantities with the dilogarithm */

}

void testclass::tearDown() {
    delete myOLEW;        
    delete myEWSMC;
    delete mySM;    
}

void testclass::setSMparameters(StandardModel& SM_i) {
    std::map<std::string, double> Parameters;
    // 17 parameters defined in StandardModel
    Parameters["GF"] = 1.16637E-5;
    Parameters["mneutrino_1"] = 0.0;
    Parameters["mneutrino_2"] = 0.0;
    Parameters["mneutrino_3"] = 0.0;
    Parameters["melectron"] = 0.54857990943e-3;
    Parameters["mmu"] = 0.1134289256;
    Parameters["mtau"] = 1.77682;
    Parameters["lambda"] = 0.2253;
    Parameters["A"] = 0.808;
    Parameters["rhob"] = 0.132;
    Parameters["etab"] = 0.341;
    Parameters["ale"] = 1.0/137.035999679;
    Parameters["dAle5Mz"] = 0.02758;
    Parameters["mHl"] = 130.0;
    Parameters["muw"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
    // 26 parameters defined in QCD    
    Parameters["AlsMz"] = 0.1184;
    Parameters["Mz"] = 91.1876;
    Parameters["mup"] = 0.003;
    Parameters["mdown"] = 0.007;
    Parameters["mcharm"] = 1.5;
    Parameters["mstrange"] = 0.1;
    Parameters["mtop"] = 174.0;
    Parameters["mbottom"] = 4.28;
    Parameters["mut"] = 0.0;
    Parameters["mub"] = 0.0;
    Parameters["muc"] = 0.0;
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
 
    /** Test for alpha_lep **/
    //Parameters["melectron"] = 0.00051099907;
    //Parameters["mmu"] = 0.105658389;
    //Parameters["mtau"] = 1.777;    
    //Parameters["AlsMz"] = 1.0/137.0359895;
    //Parameters["Mz"] = 91.187;    
    
    
    /** To make comparisons with ZFITTER codes **/
    Parameters["GF"] = 1.16637E-5;  // for GFER=2
    Parameters["melectron"] = 0.51099907e-3;
    Parameters["mmu"] = 0.105658389;
    Parameters["mtau"] = 1.77705;
    Parameters["mup"] = 0.062;
    Parameters["mdown"] = 0.083;
    Parameters["mcharm"] = 1.50;
    Parameters["mstrange"] = 0.215;
    Parameters["mbottom"] = 4.70;
    Parameters["ale"] = 1.0/137.0359895;


    /* TEST for Table 6.1 in hep-ph/0507146*/
    /* flags: AMT4=6, ALEM=2 */
    /* mcMz = 0.55381685, mbMz = 2.8194352 */
    Parameters["Mz"] = 91.1875;
    Parameters["mtop"] = 175.0;    
    Parameters["mHl"] = 150.0;
    Parameters["AlsMz"] = 0.118;
    Parameters["dAle5Mz"] = 0.02758;    
    
    
    SM_i.Init(Parameters);
}

void testclass::DeltaAlpha_l() {
    double ZFITTER = 0.031418982830747; /* ZFITTER result*/
    double result = myOLEW->DeltaAlpha_l();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void testclass::SigmaWW_bos_Mw_0() {
    double W0 = -1.702845692408391; /* ZFITTER result*/
    double result = myOLEW->SigmaWW_bos(Mw,0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0*Mw2, result, delta);
}

void testclass::SigmaWW_fer_Mw_0() {
    double W0F = 7.306899137576438; /* ZFITTER result*/
    double result = myOLEW->SigmaWW_fer(Mw,0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(W0F*Mw2, result, delta);
}

void testclass::SigmaWW_bos_Mw_Mw2_real() {
    double XMM1 = -4.212516913677337; /* ZFITTER result*/
    double result = myOLEW->SigmaWW_bos(Mw,Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1*Mw2, result, delta);
}

void testclass::SigmaWW_bos_Mw_Mw2_imag() {
    double XMM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myOLEW->SigmaWW_bos(Mw,Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1*Mw2, result, delta);
}

void testclass::SigmaWW_fer_Mw_Mw2_real() {
    double XMM1F = 11.188122745796212; /* ZFITTER result*/
    double result = myOLEW->SigmaWW_fer(Mw,Mw2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1F*Mw2, result, delta);   
}

void testclass::SigmaWW_fer_Mw_Mw2_imag() {
    double XMM1F = 9.422358620516286; /* ZFITTER result*/
    double result = myOLEW->SigmaWW_fer(Mw,Mw2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XMM1F*Mw2, result, delta); 
}

void testclass::SigmaZZ_bos_Mw_Mz2_real() {
    double XZM1 = -3.320247142827423; /* ZFITTER result*/
    double result = myOLEW->SigmaZZ_bos(Mw,Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1*Mw2, result, delta);  
}

void testclass::SigmaZZ_bos_Mw_Mz2_imag() {
    double XZM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myOLEW->SigmaZZ_bos(Mw,Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1*Mw2, result, delta);     
}

void testclass::SigmaZZ_fer_Mw_Mz2_real() {
    double XZM1F = 14.431690239155724; /* ZFITTER result*/
    double result = myOLEW->SigmaZZ_fer(Mw,Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1F*Mw2, result, delta);  
}

void testclass::SigmaZZ_fer_Mw_Mz2_imag() {
    double XZM1F = 9.893711064659454; /* ZFITTER result*/
    double result = myOLEW->SigmaZZ_fer(Mw,Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XZM1F*Mw2, result, delta); 
}

void testclass::PiZgamma_bos_Mw_Mz2_real() {
    double XAMM1 = -1.999074211569152; /* ZFITTER result*/
    double result = myOLEW->PiZgamma_bos(Mw,Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1, result, delta);    
}

void testclass::PiZgamma_bos_Mw_Mz2_imag() {
    double XAMM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myOLEW->PiZgamma_bos(Mw,Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1, result, delta); 
}

void testclass::PiZgamma_fer_Mw_Mz2_real() {
    double XAMM1F = 1.641562584085382; /* ZFITTER result*/
    double result = myOLEW->PiZgamma_fer(Mw,Mz2).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1F, result, delta);    
}

void testclass::PiZgamma_fer_Mw_Mz2_imag() {
    double XAMM1F = 4.547472366444150; /* ZFITTER result*/
    double result = myOLEW->PiZgamma_fer(Mw,Mz2).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XAMM1F, result, delta); 
}

void testclass::PiGammaGamma_fer_Mw_0() {
    double ZFITTER = 110.687763755483232; /* ZFITTER result*/
    double result = myOLEW->PiGammaGamma_fer(Mw,0.0).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::B0_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B0(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->B0(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void testclass::B1_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B1(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->B1(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = - log(cW2)/2.0;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void testclass::B21_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B21(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->B21(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = log(cW2)/3.0;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void testclass::Bf_diff_Mw2_Mz_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->Bf(Mw, Mw2, Mz, Mw).real();
    double result_Mz = myPV->Bf(Mz, Mw2, Mz, Mw).real();
    double MZtoMW = - log(cW2)/3.0;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void testclass::B0p_diff_Mw2_0_Mw() {
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result_Mw = myPV->B0p(Mw, Mw2, 0.0, Mw).real();
    double result_Mz = myPV->B0p(Mz, Mw2, 0.0, Mw).real();
    double MZtoMW = - log(cW2)/2.0/Mw2;
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);   
    delete myPV;
}

void testclass::SigmaWW_bos_diff_0_real() {
    double result_Mw = myOLEW->SigmaWW_bos(Mw,0.0).real();
    double result_Mz = myOLEW->SigmaWW_bos(Mz,0.0).real();
    double MZtoMW = Mw2*(3.0/4.0/cW2 - 3.0/2.0 + 3.0/4.0*cW2 + 3.0/4.0*cW2*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::SigmaWW_fer_diff_0_real() {
    double result_Mw = myOLEW->SigmaWW_fer(Mw,0.0).real();
    double result_Mz = myOLEW->SigmaWW_fer(Mz,0.0).real();
    double MZtoMW = 0.0;
    for (int i=0; i<6; i++) {
        MZtoMW += - 1.0/2.0*pow(mySM->getLeptons(i).getMass(),2.0)*log(cW2);
        MZtoMW += - 3.0/2.0*pow(mySM->getQuarks(i).getMass(),2.0)*log(cW2);
    }
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::SigmaWW_bos_diff_Mw2_real() {
    double result_Mw = myOLEW->SigmaWW_bos(Mw,Mw2).real();
    double result_Mz = myOLEW->SigmaWW_bos(Mz,Mw2).real();
    double MZtoMW = Mw2*(3.0/4.0/cW2 - 17.0/3.0)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::SigmaWW_fer_diff_Mw2_real() {
    double result_Mw = myOLEW->SigmaWW_fer(Mw,Mw2).real();
    double result_Mz = myOLEW->SigmaWW_fer(Mz,Mw2).real();
    double MZtoMW = 24.0/6.0*Mw2*log(cW2);
    for (int i=0; i<6; i++) {
        MZtoMW += - 1.0/2.0*pow(mySM->getLeptons(i).getMass(),2.0)*log(cW2);
        MZtoMW += - 3.0/2.0*pow(mySM->getQuarks(i).getMass(),2.0)*log(cW2);
    }
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::SigmaZZ_bos_diff_Mz2_real() {
    double result_Mw = myOLEW->SigmaZZ_bos(Mw,Mz2).real();
    double result_Mz = myOLEW->SigmaZZ_bos(Mz,Mz2).real();
    double MZtoMW = Mw2*(11.0/12.0/cW2 + 7.0/6.0 - 7.0*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::SigmaZZ_fer_diff_Mz2_real() {
    double result_Mw = myOLEW->SigmaZZ_fer(Mw,Mz2).real();
    double result_Mz = myOLEW->SigmaZZ_fer(Mz,Mz2).real();
    double MZtoMW = Mw2*(24.0/6.0/cW2 - 2.0*sW2/3.0/cW2*12.0 
                         + 4.0*sW2*sW2/3.0/cW2*8.0)*log(cW2);
    for (int i=0; i<6; i++) {
        MZtoMW += - 1.0/2.0*pow(mySM->getLeptons(i).getMass(),2.0)*log(cW2);
        MZtoMW += - 3.0/2.0*pow(mySM->getQuarks(i).getMass(),2.0)*log(cW2);
    }
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::PiGammaGamma_fer_diff_0_real() {
    double result_Mw = myOLEW->PiGammaGamma_fer(Mw,0.0).real();
    double result_Mz = myOLEW->PiGammaGamma_fer(Mz,0.0).real();
    double MZtoMW = 4.0/3.0*8.0*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::PiGammaGamma_bos_diff_Mz2_real() {
    double result_Mw = myOLEW->PiGammaGamma_bos(Mw,Mz2).real();
    double result_Mz = myOLEW->PiGammaGamma_bos(Mz,Mz2).real();
    double MZtoMW = 1.0/cW2*(1.0/12.0/cW2 + 7.0/6.0 - 7.0*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::PiGammaGamma_fer_diff_Mz2_real() {
    double result_Mw = myOLEW->PiGammaGamma_fer(Mw,Mz2).real();
    double result_Mz = myOLEW->PiGammaGamma_fer(Mz,Mz2).real();
    double MZtoMW = 4.0/3.0*8.0*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::PiZgamma_bos_diff_Mz2_real() {
    double result_Mw = myOLEW->PiZgamma_bos(Mw,Mz2).real();
    double result_Mz = myOLEW->PiZgamma_bos(Mz,Mz2).real();
    double MZtoMW = - (1.0/12.0/cW2 + 7.0/6.0 - 7.0*cW2)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::PiZgamma_fer_diff_Mz2_real() {
    double result_Mw = myOLEW->PiZgamma_fer(Mw,Mz2).real();
    double result_Mz = myOLEW->PiZgamma_fer(Mz,Mz2).real();
    double MZtoMW = - (12.0/3.0 - 4.0/3.0*sW2*8.0)*log(cW2);
    double delta = fabs(epsilon*result_Mw);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(result_Mw, result_Mz + MZtoMW, delta);
}

void testclass::DeltaRhobar_bos_Mw() {
    double ZFITTER = -0.892269770849914; /* ZFITTER result*/    
    double result = (myOLEW->SigmaWW_bos(Mw,Mw2).real()
                     - myOLEW->SigmaZZ_bos(Mw,Mz2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void testclass::DeltaRhobar_fer_Mw() {
    double ZFITTER = -3.243567493359512; /* ZFITTER result*/    
    double result = (myOLEW->SigmaWW_fer(Mw,Mw2).real()
                     - myOLEW->SigmaZZ_fer(Mw,Mz2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::DeltaRhobarW_bos_Mw() {
    double ZFITTER = 2.509671221268946; /* ZFITTER result*/    
    double result = (myOLEW->SigmaWW_bos(Mw,0.0).real()
                     - myOLEW->SigmaWW_bos(Mw,Mw2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void testclass::DeltaRhobarW_fer_Mw() {
    double ZFITTER = -3.881223608219774; /* ZFITTER result*/    
    double result = (myOLEW->SigmaWW_fer(Mw,0.0).real()
                     - myOLEW->SigmaWW_fer(Mw,Mw2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}
    
void testclass::DeltaRhobar_bos_Mz() {
    double ZFITTER = -1.257094198330098; /* ZFITTER result*/    
    double result = (myOLEW->SigmaWW_bos(Mz,Mw2).real()
                     - myOLEW->SigmaZZ_bos(Mz,Mz2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void testclass::DeltaRhobar_fer_Mz() {
    double ZFITTER = -3.132030503808709; /* ZFITTER result*/    
    double result = (myOLEW->SigmaWW_fer(Mz,Mw2).real()
                     - myOLEW->SigmaZZ_fer(Mz,Mz2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::TEST_DeltaRhobar_bos_Mw() {
    double test_result = myOLEW->TEST_DeltaRhobar_bos();
    double result = (myOLEW->SigmaWW_bos(Mw,Mw2).real()
                     - myOLEW->SigmaZZ_bos(Mw,Mz2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_result, result, delta);    
}

void testclass::TEST_DeltaRhobarW_bos_Mw() {
    double test_result = myOLEW->TEST_DeltaRhobarW_bos();
    double result = (myOLEW->SigmaWW_bos(Mw,0.0).real()
                     - myOLEW->SigmaWW_bos(Mw,Mw2).real())/Mw2;
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(test_result, result, delta);    
}

void testclass::SCALE() {
    double SCALE = -0.002574576533187; /* ZFITTER result*/    
    double result = - myEWSMC->GetSM().getAle()/4.0/M_PI/myEWSMC->GetSW2()
                      * (myOLEW->SigmaWW_bos(Mw,Mw2).real()
                         + myOLEW->SigmaWW_fer(Mw,Mw2).real()
                         - myOLEW->SigmaZZ_bos(Mw,Mz2).real()
                         - myOLEW->SigmaZZ_fer(Mw,Mz2).real())/Mw2
                    - myOLEW->DeltaRho();
    result *= myEWSMC->GetCW2()/myEWSMC->GetSW2();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(SCALE, result, delta);         
}
 
void testclass::DeltaRho() {
    double ZFITTER = 0.012013822484982; /* ZFITTER result*/  
    double result = myOLEW->DeltaRho();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
} 

void testclass::DeltaR_rem() {
    double DRREMN = 0.011412292032499; /* ZFITTER result*/   
    double result = myOLEW->DeltaR_rem();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DRREMN, result, delta);
}
    
void testclass::SigmaPrimeWW_bos_Mw_Mw2_real() {
    double XWFM1 = 2.066448743517876; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_WW_bos_Mw2(Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1, result, delta);
}

void testclass::SigmaPrimeWW_bos_Mw_Mw2_imag() {
    double XWFM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_WW_bos_Mw2(Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1, result, delta);
}

void testclass::SigmaPrimeWW_fer_Mw_Mw2_real() {
    double XWFM1F = -0.966097280796035; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_WW_fer_Mw2(Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F, result, delta);
}

void testclass::SigmaPrimeWW_fer_Mw_Mw2_imag() {
    double XWFM1F = -9.424777230871838; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_WW_fer_Mw2(Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F, result, delta);
}

void testclass::SigmaPrimeZZ_bos_Mw_Mz2_real() {
    double XWFM1 = 3.109040605851439; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_ZZ_bos_Mz2(Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1*cW2, result, delta);
}

void testclass::SigmaPrimeZZ_bos_Mw_Mz2_imag() {
    double XWFM1 = 0.000000000000000; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_ZZ_bos_Mz2(Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1*cW2, result, delta);
}

void testclass::SigmaPrimeZZ_fer_Mw_Mz2_real() {
    double XWFM1F = -0.482507863333779; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_ZZ_fer_Mz2(Mw).real();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F*cW2, result, delta);
}

void testclass::SigmaPrimeZZ_fer_Mw_Mz2_imag() {
    double XWFM1F = -9.911978176545077; /* ZFITTER result*/
    double result = myOLEW->SigmaPrime_ZZ_fer_Mz2(Mw).imag();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-XWFM1F*cW2, result, delta);
}

void testclass::C0_Mz2_Mt_Mw_Mt_real() {
    double XS3T = 0.000024702795478; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mt,Mw,Mt).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T, result, delta);    
}

void testclass::C0_Mz2_Mt_Mw_Mt_imag() {
    double XS3T = 0.000000000000000; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mt,Mw,Mt).imag();
    delete myPV;
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T, result, delta);    
}

void testclass::C0_Mz2_0_Mw_0_real() {
    double XS3T0 = 0.000097066726252; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,0.0,Mw,0.0).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T0, result, delta);    
}

void testclass::C0_Mz2_0_Mw_0_imag() {
    double XS3T0 = 0.000309606030132; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,0.0,Mw,0.0).imag();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3T0, result, delta); 
}

void testclass::C0_Mz2_Mw_Mt_Mw_real() {
    double XS3W = 0.000044262165905; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,Mt,Mw).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W, result, delta);    
}

void testclass::C0_Mz2_Mw_Mt_Mw_imag() {
    double XS3W = 0.000000000000000; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,Mt,Mw).imag();
    delete myPV;
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W, result, delta);    
}

void testclass::C0_Mz2_Mw_0_Mw_real() {
    double XS3W0 = 0.000172249280944; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,0.0,Mw).real();
    delete myPV;
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W0, result, delta);    
}

void testclass::C0_Mz2_Mw_0_Mw_imag() {
    double XS3W0 = 0.000000000000000; /* ZFITTER result*/
    PVfunctions* myPV;
    myPV = new PVfunctions();
    double result = myPV->C0(Mz2,Mw,0.0,Mw).imag();
    delete myPV;
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(XS3W0, result, delta);    
}

void testclass::FZa_0_real() {
    double V1ZZ = 1.079736267392905; /* ZFITTER result*/
    double result = myOLEW->FZa_0(Mz*Mz).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1ZZ, result, delta);    
}

void testclass::FZa_0_imag() {
    double V1ZIM = 1.712725454479850; /* ZFITTER result*/
    double result = myOLEW->FZa_0(Mz*Mz).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1ZIM, result, delta);    
}
 
void testclass::FWa_0_real() {
    double V1ZW = 1.175174895095040; /* ZFITTER result*/
    double result = myOLEW->FWa_0(Mz*Mz).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1ZW, result, delta);    
}

void testclass::FWa_0_imag() {
    double V1WIM = 2.082775415576450; /* ZFITTER result*/
    double result = myOLEW->FWa_0(Mz*Mz).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V1WIM, result, delta);    
}

void testclass::FWn_0_real() {
    double V2ZWW = -1.359807910427424; /* ZFITTER result*/
    double result = myOLEW->FWn_0(Mz*Mz).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V2ZWW, result, delta);    
}

void testclass::FWn_0_imag() {
    double ZFITTER = 0.000000000000000; /* ZFITTER result*/
    double result = myOLEW->FWn_0(Mz*Mz).imag();
    //double delta = fabs(epsilon_Li2*result);
    double delta = pow(10.0, -10.0);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);    
}

void testclass::FWa_t_real() {
    double WWv11 = -0.695682293850053; /* ZFITTER result*/
    double result = myOLEW->FWa_t(Mz*Mz).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv11, result, delta);    
}

void testclass::FbarWa_t_real() {
    double WWv12 = 3.496220092867632; /* ZFITTER result*/
    double result = myOLEW->FbarWa_t(Mz*Mz).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv12, result, delta);    
}

void testclass::FWn_t_real() {
    double WWv2 = -1.460637714980599; /* ZFITTER result*/
    double result = myOLEW->FWn_t(Mz*Mz).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv2, result, delta);      
}

void testclass::TEST_FWn_0_real() {
    double V2ZWW = -1.359807910427424; /* ZFITTER result*/
    double result = myOLEW->TEST_FWn(Mz*Mz, 0.0).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(V2ZWW, result, delta);    
}

void testclass::TEST_FWn_t_real() {
    double WWv2 = -1.460637714980599; /* ZFITTER result*/
    double result = myOLEW->TEST_FWn(Mz*Mz, mySM->getQuarks(mySM->TOP).getMass()).real()
                    - myOLEW->TEST_FWn(Mz*Mz, 0.0).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(WWv2, result, delta);    
}

void testclass::FW_diff_bb_dd_real() {
    double VTB = -2.649420523952066; /* ZFITTER result*/
    double result = myOLEW->FW(Mz*Mz, mySM->BOTTOM).real()
                    - myOLEW->FW(Mz*Mz, mySM->DOWN).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(VTB, result, delta);    
}

void testclass::f_convert() {
    double RENORM = 1.000000000000000; /* ZFITTER result*/    
    double result = myEWSMC->GetF_AlphaToGF();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(RENORM, result, delta);  
}

void testclass::Xt_GF() {
    double TOPX2 = 0.003198951884234; /* ZFITTER result*/
    double result = myEWSMC->GetXt_GF();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(TOPX2, result, delta);    
}

void testclass::DeltaRho_CORRHO() {
    double CORRHO = -0.012013822484982; /* ZFITTER result*/  
    double result = myOLEW->DeltaRho();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-CORRHO, result, delta);
} 

void testclass::DeltaRho_G() {
    double DROBLO = 0.012013822484982; /* ZFITTER result*/  
    double result = myOLEW->DeltaRho()*myEWSMC->GetF_AlphaToGF();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DROBLO, result, delta);
} 

void testclass::DeltaRbar_rem() {
    double DRREMD = 0.004147533855948; /* ZFITTER result*/   
    double result = myOLEW->DeltaRbar_rem();
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(DRREMD, result, delta);
}

void testclass::deltaRho_rem_NEUTRINO_1_real() {
    double ZFITTER = -0.002782064884822; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->NEUTRINO_1).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_NEUTRINO_2_real() {
    double ZFITTER = -0.002782064884822; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->NEUTRINO_2).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_NEUTRINO_3_real() {
    double ZFITTER = -0.002782064884822; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->NEUTRINO_3).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_ELECTRON_real() {
    double ZFITTER = -0.005521418517950; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->ELECTRON).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_MU_real() {
    double ZFITTER = -0.005521418517950; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->MU).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_TAU_real() {
    double ZFITTER = -0.005521418517950; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->TAU).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_UP_real() {
    double ZFITTER = -0.004833425293458; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->UP).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_DOWN_real() {
    double ZFITTER = -0.003920307415749; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->DOWN).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_CHARM_real() {
    double ZFITTER = -0.004833425293458; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->CHARM).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_STRANGE_real() {
    double ZFITTER = -0.003920307415749; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->STRANGE).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_NEUTRINO_1_imag() {
    double ZFITTER = -0.000306748406150; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->NEUTRINO_1).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_NEUTRINO_2_imag() {
    double ZFITTER = -0.000306748406150; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->NEUTRINO_2).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_NEUTRINO_3_imag() {
    double ZFITTER = -0.000306748406150; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->NEUTRINO_3).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_ELECTRON_imag() {
    double ZFITTER = -0.004905987392756; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->ELECTRON).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_MU_imag() {
    double ZFITTER = -0.004905987392756; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->MU).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_TAU_imag() {
    double ZFITTER = -0.004905987392756; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_l(mySM->TAU).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_UP_imag() {
    double ZFITTER = -0.003730010419090; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->UP).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_DOWN_imag() {
    double ZFITTER = -0.002196930756888; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->DOWN).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_CHARM_imag() {
    double ZFITTER = -0.003730010419090; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->CHARM).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_STRANGE_imag() {
    double ZFITTER = -0.002196930756888; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->STRANGE).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_NEUTRINO_1_real() {
    double ZFITTER = -0.000632677138537; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->NEUTRINO_1).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_NEUTRINO_2_real() {
    double ZFITTER = -0.000632677138537; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->NEUTRINO_2).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_NEUTRINO_3_real() {
    double ZFITTER = -0.000632677138537; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->NEUTRINO_3).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_ELECTRON_real() {
    double ZFITTER = 0.000905843167939; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->ELECTRON).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_MU_real() {
    double ZFITTER = 0.000905843167939; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->MU).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_TAU_real() {
    double ZFITTER = 0.000905843167939; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->TAU).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_UP_real() {
    double ZFITTER = 0.000468044616853; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->UP).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_DOWN_real() {
    double ZFITTER = -0.000044795485306; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->DOWN).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_CHARM_real() {
    double ZFITTER = 0.000468044616853; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->CHARM).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_STRANGE_real() {
    double ZFITTER = -0.000044795485306; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->STRANGE).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_NEUTRINO_1_imag() {
    double ZFITTER = 0.012600622680476; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->NEUTRINO_1).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_NEUTRINO_2_imag() {
    double ZFITTER = 0.012600622680476; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->NEUTRINO_2).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_NEUTRINO_3_imag() {
    double ZFITTER = 0.012600622680476; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->NEUTRINO_3).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_ELECTRON_imag() {
    double ZFITTER = 0.015168069190180; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->ELECTRON).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_MU_imag() {
    double ZFITTER = 0.015168069190180; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->MU).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_TAU_imag() {
    double ZFITTER = 0.015168069190180; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_l(mySM->TAU).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_UP_imag() {
    double ZFITTER = 0.014431287916457; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->UP).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_DOWN_imag() {
    double ZFITTER = 0.013575472413222; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->DOWN).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_CHARM_imag() {
    double ZFITTER = 0.014431287916457; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->CHARM).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_STRANGE_imag() {
    double ZFITTER = 0.013575472413222; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->STRANGE).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_BOTTOM_real() {
    double ZFITTER = -0.018424186901021; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->BOTTOM).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaRho_rem_BOTTOM_imag() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myOLEW->deltaRho_rem_q(mySM->BOTTOM).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_BOTTOM_real() {
    double ZFITTER = 0.007207144257329; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->BOTTOM).real();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::deltaKappa_rem_BOTTOM_imag() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myOLEW->deltaKappa_rem_q(mySM->BOTTOM).imag();
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);
}

void testclass::rho_GammaW_NuE_E() {
    double ZFITTER = 0.0; /* ZFITTER result*/
    double result = myOLEW->rho_GammaW(mySM->NEUTRINO_1, mySM->ELECTRON);
    double delta = fabs(epsilon_Li2*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}




