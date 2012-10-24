/* 
 * File:   LEP2sigmaCharm.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaCharm.h"


double LEP2sigmaCharm::getThValue() { 
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_q();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaCharm::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-13); // desired absolute error
            ig.SetRelTolerance(1.E-5); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "sigma = " << SMresult_cache << std::endl;
        }
        
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaCharm::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-14); // desired absolute error
            ig_box.SetRelTolerance(1.E-5); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            //std::cout << "sigma_box = " << sigma_box << std::endl;
            SMresult_cache += sigma_box;
        }
    }
    double sigma_charm = SMresult_cache;
    
    if ( myEW.checkModelForSTU() && !bSigmaForAFB)
        sigma_charm += myLEP2oblique.sigma_q_LEP2_NP(StandardModel::CHARM, s);
    
    return ( sigma_charm*GeVminus2_to_nb*1000.0 );
}
        
