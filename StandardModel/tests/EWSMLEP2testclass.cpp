/*
 * File:   EWSMLEP2testclass.cpp
 * Author: mishima
 */

#include "EWSMLEP2testclass.h"


CPPUNIT_TEST_SUITE_REGISTRATION(EWSMLEP2testclass);

EWSMLEP2testclass::EWSMLEP2testclass() {
}

EWSMLEP2testclass::~EWSMLEP2testclass() {
}

void EWSMLEP2testclass::setUp() {
    const bool bDebug = true;
    mySM = new StandardModel(bDebug);
    setModelParameters(*mySM);
    mySM->InitializeModel();
    myLEP2 = new EWSMTwoFermionsLEP2(*mySM, bDebug);

    sqrt_s = 200.0;
    Mw = 80.360848365211552;
    // alpha(s)=0.007821563447202 
    
    
    Mw2 = Mw*Mw;
    Mz = mySM->getMz();
    Mz2 = Mz*Mz;
    cW2 = Mw2/Mz2;
    sW2 = 1.0 - cW2;
    Mt = mySM->getMtpole();
    
    /* accuracy for CPPUNIT_ASSERT_DOUBLES_EQUAL */
    epsilon = 1.0e-7; 
}

void EWSMLEP2testclass::tearDown() {
    delete myLEP2; 
    delete mySM; 
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
    
    Model_i.Init(Parameters);
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


//XCHI1 =         0.4714952963       -0.0162861206