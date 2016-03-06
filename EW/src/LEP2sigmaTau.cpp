/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaTau.h"


double LEP2sigmaTau::computeThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = SM.Gamma_Z();

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
    }
    double sigma_tau = SMresult_cache;
    
    #ifdef LEP2TEST
    sigma_tau = myTEST.sigmaTauTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif 
    
    if ( checkLEP2NP() && !bSigmaForAFB ) {
        double obliqueShat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueShat();
        double obliqueThat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueThat();
        double obliqueUhat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueUhat();
        double obliqueV = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueV();
        double obliqueW = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueW();
        double obliqueX = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueX();
        double obliqueY = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueY();
        double ObParam[7] = {obliqueShat, obliqueThat, obliqueUhat,
                             obliqueV, obliqueW, obliqueX, obliqueY};
        sigma_tau += myLEP2oblique.sigma_l_LEP2_NP(StandardModel::TAU, s, ml_cache, ObParam);
    }
    
    if ( checkLEP2GIMR() && !bSigmaForAFB ) {
        double myCLL = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLL_tau();
        double myCLR = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLR_tau();
        double myCRL = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRL_tau();
        double myCRR = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRR_tau();
        
        double deltaGammaZ = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGamma_Z();
        double deltaGR_l = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getLeptons(StandardModel::TAU));
        double deltaGL_l = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getLeptons(StandardModel::TAU));
        double deltaGR_e = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getLeptons(StandardModel::ELECTRON));
        double deltaGL_e = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getLeptons(StandardModel::ELECTRON));
        double deltaMz2=(dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaMz2();
        
        double GIMRParam[10] ={myCLL,myCLR,myCRL,myCRR,deltaGammaZ,deltaGL_l,deltaGR_l,deltaGL_e,deltaGR_e,deltaMz2};
        
        sigma_tau += myLEP2GIMR.sigma_l_LEP2_GIMR(StandardModel::TAU, s, GIMRParam);
    }
    
    return ( sigma_tau*SM.GeVminus2_to_nb*1000.0 );
}
        

