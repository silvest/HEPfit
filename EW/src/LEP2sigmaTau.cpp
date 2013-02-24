/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaTau.h"


double LEP2sigmaTau::getThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        ml_cache = m_l(SM.TAU);        

        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_l();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaTau::Integrand_sigmaWithISR_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << " +- " << ig.Error() << std::endl;
            // results: 6.8e-9 -- 2.2e-8
        }
                
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaTau::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-17); // desired absolute error
            ig_box.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            SMresult_cache += sigma_box;
            //std::cout << sigma_box << " +- " << ig_box.Error() << std::endl;
            // results: 3.6e-12 -- 2.9e-10
        }

        if ( myEW.checkModelForSTU() && !bSigmaForAFB && SM.IsFlagFixedAllSMparams() ) {
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.sigma_l_LEP2_NP(StandardModel::TAU, s, ml_cache, ObParam);
            }
        }
    }
    double sigma_tau = SMresult_cache;
    
    #ifdef LEP2TEST
    sigma_tau = myTEST.sigmaTauTEST(sqrt_s)/GeVminus2_to_nb/1000.0;
    #endif 
    
    if ( myEW.checkModelForSTU() && !bSigmaForAFB ) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            sigma_tau += Coeff_cache[myLEP2oblique.Shat]*myEW.Shat()
                       + Coeff_cache[myLEP2oblique.That]*myEW.That()
                       + Coeff_cache[myLEP2oblique.Uhat]*myEW.Uhat()
                       + Coeff_cache[myLEP2oblique.V]*myEW.V()
                       + Coeff_cache[myLEP2oblique.W]*myEW.W()
                       + Coeff_cache[myLEP2oblique.X]*myEW.X()
                       + Coeff_cache[myLEP2oblique.Y]*myEW.Y();
        } else {
            double ObParam[7] = {myEW.Shat(), myEW.That(), myEW.Uhat(),
                                 myEW.V(), myEW.W(), myEW.X(), myEW.Y()};
            sigma_tau += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::TAU, s, ml_cache, ObParam);
        }
    }
    
    return ( sigma_tau*GeVminus2_to_nb*1000.0 );
}
        

