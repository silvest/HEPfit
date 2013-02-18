/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaBottom.h"


double LEP2sigmaBottom::getThValue() { 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mq_cache = m_q(SM.BOTTOM, sqrt_s);
                
        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_q();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaBottom::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << " +- " << ig.Error() << std::endl;
        }
        
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaBottom::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            SMresult_cache += sigma_box;
            //std::cout << sigma_box << " +- " << ig_box.Error() << std::endl;
        }

        if ( myEW.checkModelForSTU() && !bSigmaForAFB && SM.IsFlagFixedAllSMparams()) {
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.sigma_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
            }
        }
    }
    double sigma_bottom = SMresult_cache;
    
    if ( myEW.checkModelForSTU() && !bSigmaForAFB) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            sigma_bottom += Coeff_cache[myLEP2oblique.Shat]*myEW.Shat()
                          + Coeff_cache[myLEP2oblique.That]*myEW.That()
                          + Coeff_cache[myLEP2oblique.Uhat]*myEW.Uhat()
                          + Coeff_cache[myLEP2oblique.V]*myEW.V()
                          + Coeff_cache[myLEP2oblique.W]*myEW.W()
                          + Coeff_cache[myLEP2oblique.X]*myEW.X()
                          + Coeff_cache[myLEP2oblique.Y]*myEW.Y();
        } else {
            double ObParam[7] = {myEW.Shat(), myEW.That(), myEW.Uhat(),
                                 myEW.V(), myEW.W(), myEW.X(), myEW.Y()};
            sigma_bottom += myLEP2oblique.sigma_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
        }
    }
    
    return ( sigma_bottom*GeVminus2_to_nb*1000.0 );
}
        
