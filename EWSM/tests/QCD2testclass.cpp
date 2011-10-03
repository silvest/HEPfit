/*
 * File:   QCD2testclass.cpp
 * Author: mishima
 */

#include "QCD2testclass.h"


CPPUNIT_TEST_SUITE_REGISTRATION(QCD2testclass);

QCD2testclass::QCD2testclass() {
}

QCD2testclass::~QCD2testclass() {
}

void QCD2testclass::setUp() {
    mySM = new StandardModel();
    QCD2testclass::setSMparameters(*mySM);   
    
    myEWSMC = new EWSMcommon(*mySM);
    myEWSMC->SetConstants();
    myEWSMC->Compute(mySM->Mw_tree());
    
    myQCD2 = new TwoLoopQCD(*myEWSMC);

    Mw = myEWSMC->GetMw();
    Mw2 = Mw*Mw;
    Mz = myEWSMC->GetSM().getMz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = myEWSMC->GetSM().getQuarks(StandardModel::TOP).getMass();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 

}

void QCD2testclass::tearDown() {
    delete myQCD2;        
    delete myEWSMC;
    delete mySM;  
}

void QCD2testclass::setSMparameters(StandardModel& SM_i) {
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

void QCD2testclass::F1() {
    double ZFITTER = -1.027558816350545; /* ZFITTER result*/
    double result = myQCD2->F1(Mw*Mw/Mt/Mt);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void QCD2testclass::F1_0() {
    double ZFITTER = -1.188052388163498; /* ZFITTER result*/
    double result = myQCD2->F1(0.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void QCD2testclass::V1() {
    double ZFITTER = 0.289346402688953; /* ZFITTER result*/
    double result = myQCD2->V1(Mz*Mz/4.0/Mt/Mt);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void QCD2testclass::A1() {
    double ZFITTER = -6.747476785054744; /* ZFITTER result*/
    double result = myQCD2->A1(Mz*Mz/4.0/Mt/Mt);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}

void QCD2testclass::A1_0() {
    double ZFITTER = -6.897143619502220; /* ZFITTER result*/
    double result = myQCD2->A1(0.0);
    double delta = fabs(epsilon*result);
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ZFITTER, result, delta);   
}




