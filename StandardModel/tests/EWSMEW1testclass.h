/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#ifndef EWSMEW1TESTCLASS_H
#define	EWSMEW1TESTCLASS_H

#include <cppunit/extensions/HelperMacros.h>

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <map>
#include "EWSMOneLoopEW.h"
using namespace std;


class EWSMEW1testclass : public CPPUNIT_NS::TestFixture {
    CPPUNIT_TEST_SUITE(EWSMEW1testclass);
    CPPUNIT_TEST(DeltaAlpha_l);
    CPPUNIT_TEST(SigmaWW_bos_Mw_0);
    CPPUNIT_TEST(SigmaWW_fer_Mw_0);
    CPPUNIT_TEST(SigmaWW_bos_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaWW_bos_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaWW_fer_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaWW_fer_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaZZ_bos_Mw_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_bos_Mw_Mz2_imag);
    CPPUNIT_TEST(SigmaZZ_fer_Mw_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_fer_Mw_Mz2_imag);
    CPPUNIT_TEST(PiZgamma_bos_Mw_Mz2_real);
    CPPUNIT_TEST(PiZgamma_bos_Mw_Mz2_imag);
    CPPUNIT_TEST(PiZgamma_fer_Mw_Mz2_real);
    CPPUNIT_TEST(PiZgamma_fer_Mw_Mz2_imag);
    CPPUNIT_TEST(PiGammaGamma_fer_Mw_0);    
    CPPUNIT_TEST(B0_diff_Mw2_Mz2_Mw2);
    CPPUNIT_TEST(B1_diff_Mw2_Mz2_Mw2);
    CPPUNIT_TEST(B11_diff_Mw2_Mz2_Mw2);
    CPPUNIT_TEST(Bf_diff_Mw2_Mz2_Mw2);
    CPPUNIT_TEST(B0p_diff_Mw2_0_Mw2);
    CPPUNIT_TEST(SigmaWW_bos_diff_0_real);
    CPPUNIT_TEST(SigmaWW_fer_diff_0_real);    
    CPPUNIT_TEST(SigmaWW_bos_diff_Mw2_real);
    CPPUNIT_TEST(SigmaWW_fer_diff_Mw2_real);
    CPPUNIT_TEST(SigmaZZ_bos_diff_Mz2_real);
    CPPUNIT_TEST(SigmaZZ_fer_diff_Mz2_real);
    CPPUNIT_TEST(PiGammaGamma_fer_diff_0_real);
    CPPUNIT_TEST(PiGammaGamma_bos_diff_Mz2_real);
    CPPUNIT_TEST(PiGammaGamma_fer_diff_Mz2_real);
    CPPUNIT_TEST(PiZgamma_bos_diff_Mz2_real);
    CPPUNIT_TEST(PiZgamma_fer_diff_Mz2_real);    
    CPPUNIT_TEST(DeltaRhobar_bos_Mw);
    CPPUNIT_TEST(DeltaRhobar_fer_Mw);
    CPPUNIT_TEST(DeltaRhobarW_bos_Mw);
    CPPUNIT_TEST(DeltaRhobarW_fer_Mw);
    CPPUNIT_TEST(DeltaRhobar_bos_Mz);
    CPPUNIT_TEST(DeltaRhobar_fer_Mz);
    CPPUNIT_TEST(TEST_DeltaRhobar_bos_Mw);
    CPPUNIT_TEST(TEST_DeltaRhobarW_bos_Mw);
    CPPUNIT_TEST(SCALE);
    CPPUNIT_TEST(DeltaRho);
    CPPUNIT_TEST(DeltaR_rem);
    CPPUNIT_TEST(SigmaPrimeWW_bos_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaPrimeWW_bos_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaPrimeWW_fer_Mw_Mw2_real);
    CPPUNIT_TEST(SigmaPrimeWW_fer_Mw_Mw2_imag);
    CPPUNIT_TEST(SigmaPrimeZZ_bos_Mw_Mz2_real); 
    CPPUNIT_TEST(SigmaPrimeZZ_bos_Mw_Mz2_imag); 
    CPPUNIT_TEST(SigmaPrimeZZ_fer_Mw_Mz2_real); 
    CPPUNIT_TEST(SigmaPrimeZZ_fer_Mw_Mz2_imag);
    CPPUNIT_TEST(C0_Mz2_Mt2_Mw2_Mt2_real);
    CPPUNIT_TEST(C0_Mz2_Mt2_Mw2_Mt2_imag);
    CPPUNIT_TEST(C0_Mz2_0_Mw2_0_real);
    CPPUNIT_TEST(C0_Mz2_0_Mw2_0_imag);
    CPPUNIT_TEST(C0_Mz2_Mw2_Mt2_Mw2_real);
    CPPUNIT_TEST(C0_Mz2_Mw2_Mt2_Mw2_imag);
    CPPUNIT_TEST(C0_Mz2_Mw2_0_Mw2_real);
    CPPUNIT_TEST(C0_Mz2_Mw2_0_Mw2_imag);
    CPPUNIT_TEST(FZa_0_real);
    CPPUNIT_TEST(FZa_0_imag);
    CPPUNIT_TEST(FWa_0_real);
    CPPUNIT_TEST(FWa_0_imag);
    CPPUNIT_TEST(FWn_0_real);
    CPPUNIT_TEST(FWn_0_imag);
    CPPUNIT_TEST(FWa_t_real);
    CPPUNIT_TEST(FbarWa_t_real);
    CPPUNIT_TEST(FWn_t_real);
    CPPUNIT_TEST(TEST_FWn_0_real);    
    CPPUNIT_TEST(TEST_FWn_t_real);    
    CPPUNIT_TEST(FW_diff_bb_dd_real);
    CPPUNIT_TEST(f_convert);
    CPPUNIT_TEST(Xt_GF);
    CPPUNIT_TEST(DeltaRho_CORRHO);
    CPPUNIT_TEST(DeltaRho_G);
    CPPUNIT_TEST(DeltaRbar_rem);
    CPPUNIT_TEST(deltaRho_rem_NEUTRINO_1_real);;    
    CPPUNIT_TEST(deltaRho_rem_NEUTRINO_2_real);
    CPPUNIT_TEST(deltaRho_rem_NEUTRINO_3_real);    
    CPPUNIT_TEST(deltaRho_rem_ELECTRON_real);
    CPPUNIT_TEST(deltaRho_rem_MU_real);    
    CPPUNIT_TEST(deltaRho_rem_TAU_real);    
    CPPUNIT_TEST(deltaRho_rem_UP_real);
    CPPUNIT_TEST(deltaRho_rem_DOWN_real);    
    CPPUNIT_TEST(deltaRho_rem_CHARM_real);   
    CPPUNIT_TEST(deltaRho_rem_STRANGE_real);    
    CPPUNIT_TEST(deltaRho_rem_NEUTRINO_1_imag);    
    CPPUNIT_TEST(deltaRho_rem_NEUTRINO_2_imag);
    CPPUNIT_TEST(deltaRho_rem_NEUTRINO_3_imag);    
    CPPUNIT_TEST(deltaRho_rem_ELECTRON_imag);
    CPPUNIT_TEST(deltaRho_rem_MU_imag);    
    CPPUNIT_TEST(deltaRho_rem_TAU_imag);    
    CPPUNIT_TEST(deltaRho_rem_UP_imag);
    CPPUNIT_TEST(deltaRho_rem_DOWN_imag);    
    CPPUNIT_TEST(deltaRho_rem_CHARM_imag);   
    CPPUNIT_TEST(deltaRho_rem_STRANGE_imag);    
    CPPUNIT_TEST(deltaKappa_rem_NEUTRINO_1_real);;    
    CPPUNIT_TEST(deltaKappa_rem_NEUTRINO_2_real);
    CPPUNIT_TEST(deltaKappa_rem_NEUTRINO_3_real);    
    CPPUNIT_TEST(deltaKappa_rem_ELECTRON_real);
    CPPUNIT_TEST(deltaKappa_rem_MU_real);    
    CPPUNIT_TEST(deltaKappa_rem_TAU_real);    
    CPPUNIT_TEST(deltaKappa_rem_UP_real);
    CPPUNIT_TEST(deltaKappa_rem_DOWN_real);    
    CPPUNIT_TEST(deltaKappa_rem_CHARM_real);   
    CPPUNIT_TEST(deltaKappa_rem_STRANGE_real);    
    CPPUNIT_TEST(deltaKappa_rem_NEUTRINO_1_imag);    
    CPPUNIT_TEST(deltaKappa_rem_NEUTRINO_2_imag);
    CPPUNIT_TEST(deltaKappa_rem_NEUTRINO_3_imag);    
    CPPUNIT_TEST(deltaKappa_rem_ELECTRON_imag);
    CPPUNIT_TEST(deltaKappa_rem_MU_imag);    
    CPPUNIT_TEST(deltaKappa_rem_TAU_imag);    
    CPPUNIT_TEST(deltaKappa_rem_UP_imag);
    CPPUNIT_TEST(deltaKappa_rem_DOWN_imag);    
    CPPUNIT_TEST(deltaKappa_rem_CHARM_imag);   
    CPPUNIT_TEST(deltaKappa_rem_STRANGE_imag);
    CPPUNIT_TEST(deltaRho_rem_BOTTOM_real);   
    CPPUNIT_TEST(deltaRho_rem_BOTTOM_imag);       
    CPPUNIT_TEST(deltaKappa_rem_BOTTOM_real);
    CPPUNIT_TEST(deltaKappa_rem_BOTTOM_imag);
    CPPUNIT_TEST(rho_GammaW_leptons);
    CPPUNIT_TEST(rho_GammaW_quarks);
    
    
    CPPUNIT_TEST_SUITE_END();

public:
    EWSMEW1testclass();
    virtual ~EWSMEW1testclass();
    void setUp();
    void tearDown();

private:
    StandardModel* mySM;
    EWSMcache* myCache;
    EWSMOneLoopEW* myEW1;
    
    double epsilon, epsilon_Li2;
    double Mw, Mw2, Mz, Mz2, cW2, sW2, Mt, Mt2;
    
    void DeltaAlpha_l();
    
    void SigmaWW_bos_Mw_0();
    void SigmaWW_fer_Mw_0();    
    void SigmaWW_bos_Mw_Mw2_real();
    void SigmaWW_bos_Mw_Mw2_imag();
    void SigmaWW_fer_Mw_Mw2_real();     
    void SigmaWW_fer_Mw_Mw2_imag();     
    void SigmaZZ_bos_Mw_Mz2_real();
    void SigmaZZ_bos_Mw_Mz2_imag();
    void SigmaZZ_fer_Mw_Mz2_real();
    void SigmaZZ_fer_Mw_Mz2_imag();
    void PiZgamma_bos_Mw_Mz2_real();
    void PiZgamma_bos_Mw_Mz2_imag();
    void PiZgamma_fer_Mw_Mz2_real();
    void PiZgamma_fer_Mw_Mz2_imag();
    void PiGammaGamma_fer_Mw_0();
    
    void B0_diff_Mw2_Mz2_Mw2();
    void B1_diff_Mw2_Mz2_Mw2();
    void B11_diff_Mw2_Mz2_Mw2();
    void Bf_diff_Mw2_Mz2_Mw2();
    void B0p_diff_Mw2_0_Mw2();
    
    void SigmaWW_bos_diff_0_real();    
    void SigmaWW_fer_diff_0_real();
    void SigmaWW_bos_diff_Mw2_real();
    void SigmaWW_fer_diff_Mw2_real();
    void SigmaZZ_bos_diff_Mz2_real();
    void SigmaZZ_fer_diff_Mz2_real();
    void PiGammaGamma_fer_diff_0_real();
    void PiGammaGamma_bos_diff_Mz2_real();
    void PiGammaGamma_fer_diff_Mz2_real();
    void PiZgamma_bos_diff_Mz2_real();
    void PiZgamma_fer_diff_Mz2_real();
    
    void DeltaRhobar_bos_Mw();
    void DeltaRhobar_fer_Mw();
    void DeltaRhobarW_bos_Mw();
    void DeltaRhobarW_fer_Mw();

    void DeltaRhobar_bos_Mz();
    void DeltaRhobar_fer_Mz();    
    
    void TEST_DeltaRhobar_bos_Mw();
    void TEST_DeltaRhobarW_bos_Mw();    
    
    void SCALE();
    
    void DeltaRho();
    void DeltaR_rem();
    
    void SigmaPrimeWW_bos_Mw_Mw2_real();
    void SigmaPrimeWW_bos_Mw_Mw2_imag();
    void SigmaPrimeWW_fer_Mw_Mw2_real();
    void SigmaPrimeWW_fer_Mw_Mw2_imag();
    void SigmaPrimeZZ_bos_Mw_Mz2_real(); 
    void SigmaPrimeZZ_bos_Mw_Mz2_imag(); 
    void SigmaPrimeZZ_fer_Mw_Mz2_real(); 
    void SigmaPrimeZZ_fer_Mw_Mz2_imag();     
    
    void C0_Mz2_Mt2_Mw2_Mt2_real();
    void C0_Mz2_Mt2_Mw2_Mt2_imag();
    void C0_Mz2_0_Mw2_0_real();
    void C0_Mz2_0_Mw2_0_imag();
    void C0_Mz2_Mw2_Mt2_Mw2_real();
    void C0_Mz2_Mw2_Mt2_Mw2_imag();
    void C0_Mz2_Mw2_0_Mw2_real();
    void C0_Mz2_Mw2_0_Mw2_imag();

    void FZa_0_real();
    void FZa_0_imag();
    void FWa_0_real();
    void FWa_0_imag();
    void FWn_0_real();
    void FWn_0_imag();
    void FWa_t_real();
    void FbarWa_t_real();
    void FWn_t_real();
    void TEST_FWn_0_real();
    void TEST_FWn_t_real();
    void FW_diff_bb_dd_real();
    
    void f_convert();
    void Xt_GF();
    void DeltaRho_CORRHO();
    void DeltaRho_G();
    void DeltaRbar_rem();

    void deltaRho_rem_NEUTRINO_1_real();
    void deltaRho_rem_NEUTRINO_2_real();
    void deltaRho_rem_NEUTRINO_3_real();    
    void deltaRho_rem_ELECTRON_real();
    void deltaRho_rem_MU_real();    
    void deltaRho_rem_TAU_real();    
    void deltaRho_rem_UP_real();
    void deltaRho_rem_DOWN_real();    
    void deltaRho_rem_CHARM_real();   
    void deltaRho_rem_STRANGE_real();    
    
    void deltaRho_rem_NEUTRINO_1_imag();        
    void deltaRho_rem_NEUTRINO_2_imag();
    void deltaRho_rem_NEUTRINO_3_imag();    
    void deltaRho_rem_ELECTRON_imag();
    void deltaRho_rem_MU_imag();    
    void deltaRho_rem_TAU_imag();    
    void deltaRho_rem_UP_imag();
    void deltaRho_rem_DOWN_imag();    
    void deltaRho_rem_CHARM_imag();   
    void deltaRho_rem_STRANGE_imag(); 

    void deltaKappa_rem_NEUTRINO_1_real();
    void deltaKappa_rem_NEUTRINO_2_real();
    void deltaKappa_rem_NEUTRINO_3_real();    
    void deltaKappa_rem_ELECTRON_real();
    void deltaKappa_rem_MU_real();    
    void deltaKappa_rem_TAU_real();    
    void deltaKappa_rem_UP_real();
    void deltaKappa_rem_DOWN_real();    
    void deltaKappa_rem_CHARM_real();   
    void deltaKappa_rem_STRANGE_real();    
    
    void deltaKappa_rem_NEUTRINO_1_imag();        
    void deltaKappa_rem_NEUTRINO_2_imag();
    void deltaKappa_rem_NEUTRINO_3_imag();    
    void deltaKappa_rem_ELECTRON_imag();
    void deltaKappa_rem_MU_imag();    
    void deltaKappa_rem_TAU_imag();    
    void deltaKappa_rem_UP_imag();
    void deltaKappa_rem_DOWN_imag();    
    void deltaKappa_rem_CHARM_imag();   
    void deltaKappa_rem_STRANGE_imag();    

    void deltaRho_rem_BOTTOM_real();    
    void deltaRho_rem_BOTTOM_imag();
    void deltaKappa_rem_BOTTOM_real();    
    void deltaKappa_rem_BOTTOM_imag();    
    
    void rho_GammaW_leptons();
    void rho_GammaW_quarks();    
    
};

#endif	/* EWSMEW1TESTCLASS_H */

