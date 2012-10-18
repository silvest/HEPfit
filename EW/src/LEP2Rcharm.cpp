/* 
 * File:   LEP2Rcharm.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2Rcharm.h"


double LEP2Rcharm::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        double sigma_c, sigma_had;
        if (!bRCs[LEP2TwoFermions::ISR]) {
            sigma_c = myTwoFermions.sigma_q(StandardModel::CHARM, s, Mw, GammaZ, bRCs);
            sigma_had = myTwoFermions.sigma_q(StandardModel::UP, s, Mw, GammaZ, bRCs)
                      + myTwoFermions.sigma_q(StandardModel::DOWN, s, Mw, GammaZ, bRCs)
                      + sigma_c
                      + myTwoFermions.sigma_q(StandardModel::STRANGE, s, Mw, GammaZ, bRCs)
                      + myTwoFermions.sigma_q(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
        } else {
            ROOT::Math::Functor1D wf(this, &LEP2Rcharm::IntegrandISR_sigma_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-6); // desired relative error
            sigma_c = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << sigma_c << std::endl;            
            
            myLEP2sigmaHadron.setFlag("Weak", bRCs[LEP2TwoFermions::Weak]);
            myLEP2sigmaHadron.setFlag("WeakBox", bRCs[LEP2TwoFermions::WeakBox]);
            myLEP2sigmaHadron.setFlag("ISR", bRCs[LEP2TwoFermions::ISR]);
            myLEP2sigmaHadron.setFlag("QEDFSR", bRCs[LEP2TwoFermions::QEDFSR]);
            myLEP2sigmaHadron.setFlag("QCDFSR", bRCs[LEP2TwoFermions::QCDFSR]);
            sigma_had = myLEP2sigmaHadron.getThValue()/GeVminus2_to_nb/1000.0;            
        }
        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s);
    }
    
    return R_c;
}


