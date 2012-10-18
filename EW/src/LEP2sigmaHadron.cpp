/* 
 * File:   LEP2sigmaHadron.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaHadron.h"


double LEP2sigmaHadron::getThValue() { 
    double s = sqrt_s*sqrt_s;
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        if (!bRCs[LEP2TwoFermions::ISR]) {
            SMresult_cache = myTwoFermions.sigma_q(StandardModel::UP, s, Mw, GammaZ, bRCs)
                           + myTwoFermions.sigma_q(StandardModel::DOWN, s, Mw, GammaZ, bRCs)
                           + myTwoFermions.sigma_q(StandardModel::CHARM, s, Mw, GammaZ, bRCs)
                           + myTwoFermions.sigma_q(StandardModel::STRANGE, s, Mw, GammaZ, bRCs)
                           + myTwoFermions.sigma_q(StandardModel::BOTTOM, s, Mw, GammaZ, bRCs);
        } else {
            q_flavor = StandardModel::UP;
            ROOT::Math::Functor1D wf_UP(this, &LEP2sigmaHadron::IntegrandISR_sigma_q);
            ROOT::Math::Integrator ig_UP(wf_UP, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_UP.SetAbsTolerance(1.E-14); // desired absolute error
            ig_UP.SetRelTolerance(1.E-6); // desired relative error
            double sigma_UP = ig_UP.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "UP " << sigma_UP << std::endl;
            //
            q_flavor = StandardModel::DOWN;
            ROOT::Math::Functor1D wf_DOWN(this, &LEP2sigmaHadron::IntegrandISR_sigma_q);
            ROOT::Math::Integrator ig_DOWN(wf_DOWN, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_DOWN.SetAbsTolerance(1.E-14); // desired absolute error
            ig_DOWN.SetRelTolerance(1.E-6); // desired relative error
            double sigma_DOWN = ig_DOWN.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "DOWN " << sigma_DOWN << std::endl;
            //
            q_flavor = StandardModel::CHARM;
            ROOT::Math::Functor1D wf_CHARM(this, &LEP2sigmaHadron::IntegrandISR_sigma_q);
            ROOT::Math::Integrator ig_CHARM(wf_CHARM, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_CHARM.SetAbsTolerance(1.E-14); // desired absolute error
            ig_CHARM.SetRelTolerance(1.E-6); // desired relative error
            double sigma_CHARM = ig_CHARM.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "CHARM " << sigma_CHARM << std::endl;
            //
            q_flavor = StandardModel::STRANGE;
            ROOT::Math::Functor1D wf_STRANGE(this, &LEP2sigmaHadron::IntegrandISR_sigma_q);
            ROOT::Math::Integrator ig_STRANGE(wf_STRANGE, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_STRANGE.SetAbsTolerance(1.E-14); // desired absolute error
            ig_STRANGE.SetRelTolerance(1.E-6); // desired relative error
            double sigma_STRANGE = ig_STRANGE.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "STRANGE " << sigma_STRANGE << std::endl;
            //
            q_flavor = StandardModel::BOTTOM;
            ROOT::Math::Functor1D wf_BOTTOM(this, &LEP2sigmaHadron::IntegrandISR_sigma_q);
            ROOT::Math::Integrator ig_BOTTOM(wf_BOTTOM, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_BOTTOM.SetAbsTolerance(1.E-14); // desired absolute error
            ig_BOTTOM.SetRelTolerance(1.E-6); // desired relative error
            double sigma_BOTTOM = ig_BOTTOM.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "BOTTOM " << sigma_BOTTOM << std::endl;
            //
            SMresult_cache = sigma_UP + sigma_DOWN + sigma_CHARM 
                             + sigma_STRANGE + sigma_BOTTOM;
        }
    }
    double sigmaH = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        sigmaH += myLEP2oblique.sigma_q_LEP2_NP(StandardModel::UP, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::DOWN, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::CHARM, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::STRANGE, s)
                + myLEP2oblique.sigma_q_LEP2_NP(StandardModel::BOTTOM, s);
    }
    
    return ( sigmaH*GeVminus2_to_nb*1000.0 );
}
        

