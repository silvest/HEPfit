/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaCharm.h"


double LEP2sigmaCharm::computeThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mq_cache = m_q(SM.CHARM, sqrt_s);
        
        if (!flag[ISR])
            SMresult_cache = sigma_NoISR_q();
        else {
            ROOT::Math::Functor1D wf(this, &LEP2sigmaCharm::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-13); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            SMresult_cache = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << SMresult_cache << " +- " << ig.Error() << std::endl;
        }
        
        if (flag[WeakBox]) {
            ROOT::Math::Functor1D wf_box(this, &LEP2sigmaCharm::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box(wf_box, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box = ig_box.Integral(-1.0, 1.0); // interval
            SMresult_cache += sigma_box;
            //std::cout << sigma_box << " +- " << ig_box.Error() << std::endl;
        }
        
        if ( myEW.checkSTUVWXY() && !bSigmaForAFB && SM.IsFlagFixedAllSMparams()) {
            double ObParam[7];
             for (int i=0; i<7; i++) {
                 SetObParam((LEP2oblique::Oblique)i, ObParam);
                 Coeff_cache[i] 
                         = myLEP2oblique.sigma_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
             }
        }
    }
    double sigma_charm = SMresult_cache;
    
    if ( myEW.checkSTUVWXY() && !bSigmaForAFB) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            sigma_charm += Coeff_cache[myLEP2oblique.Shat]*SM.obliqueShat()
                         + Coeff_cache[myLEP2oblique.That]*SM.obliqueThat()
                         + Coeff_cache[myLEP2oblique.Uhat]*SM.obliqueUhat()
                         + Coeff_cache[myLEP2oblique.V]*SM.obliqueV()
                         + Coeff_cache[myLEP2oblique.W]*SM.obliqueW()
                         + Coeff_cache[myLEP2oblique.X]*SM.obliqueX()
                         + Coeff_cache[myLEP2oblique.Y]*SM.obliqueY();
        } else {
            double ObParam[7] = {SM.obliqueShat(), SM.obliqueThat(), SM.obliqueUhat(),
                                 SM.obliqueV(), SM.obliqueW(), SM.obliqueX(), SM.obliqueY()};
            sigma_charm += myLEP2oblique.sigma_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
        }
    }
    
    return ( sigma_charm*GeVminus2_to_nb*1000.0 );
}
        
