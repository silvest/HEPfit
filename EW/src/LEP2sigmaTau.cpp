/* 
 * File:   LEP2sigmaTau.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaTau.h"


double LEP2sigmaTau::getThValue() { 
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_l();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaTau::Integrand_sigmaWithISR_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-15); // desired absolute error
            ig.SetRelTolerance(1.E-6); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << std::endl;
        }
                
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaTau::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-14); // desired absolute error
            ig_box.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            //std::cout << sigma_box << std::endl;
            SMresult_cache += sigma_box;
        }
    }
    double sigma_tau = SMresult_cache;
    
    if ( myEW.checkModelForSTU() )
        sigma_tau += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::TAU, s);
    
    return ( sigma_tau*GeVminus2_to_nb*1000.0 );
}
        

