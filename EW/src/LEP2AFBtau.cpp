/* 
 * File:   LEP2AFBtau.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2AFBtau.h"


double LEP2AFBtau::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        if (!bRCs[LEP2TwoFermions::ISR])
            SMresult_cache = myTwoFermions.AFB_l(StandardModel::TAU, s, Mw, GammaZ, bRCs);
        else {
            ROOT::Math::Functor1D wf(this, &LEP2AFBtau::IntegrandISR_AFB_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-11); // desired absolute error
            ig.SetRelTolerance(1.E-6); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << " "
            //          << myTwoFermions.G_3prime_l(l_flavor, s, Mw, GammaZ, bRCs)/s
            //          << std::endl;
            
            // cross section
            ROOT::Math::Functor1D wf2(this, &LEP2AFBtau::IntegrandISR_sigma_l);
            ROOT::Math::Integrator ig2(wf2, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig2.SetAbsTolerance(1.E-15); // desired absolute error
            ig2.SetRelTolerance(1.E-6); // desired relative error
            double sigma = ig2.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << sigma << " " 
            //          << myTwoFermions.sigma_l(l_flavor, s, Mw, GammaZ, bRCs)
            //          << std::endl;
            
            SMresult_cache *= M_PI*SM.getAle()*SM.getAle()/sigma;
        }    
    }
    double AFB_tau = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        AFB_tau += myLEP2oblique.AFB_l_LEP2_NP(StandardModel::TAU, s);
    }
    
    return AFB_tau;
}
        