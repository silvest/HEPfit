/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2AFBtau.h"


double LEP2AFBtau::computeThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = SM.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        ml_cache = m_l(SM.TAU);

        double AFB_noBox, sigma = 0.0;
        if (!flag[ISR])
            AFB_noBox = AFB_NoISR_l();
        else {
            // numerator
            ROOT::Math::Functor1D wf(this, &LEP2AFBtau::Integrand_AFBnumeratorWithISR_l);
            ROOT::Math::Integrator ig(wf, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig.SetAbsTolerance(1.E-14); // desired absolute error
            ig.SetRelTolerance(1.E-4); // desired relative error
            double numerator = ig.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << numerator << " +- " << ig.Error() << std::endl;
            // results: 3.7e-9 -- 1.5e-8
            
            // denominator
            myLEP2sigmaTau.setFlags(flag);
            sigma = myLEP2sigmaTau.computeThValue()/SM.GeVminus2_to_nb/1000.0;
            
            AFB_noBox = numerator/sigma;
        }    
        SMresult_cache = AFB_noBox;
        
        if (flag[WeakBox]) {
            // numerator
            ROOT::Math::Functor1D wf_F(this, &LEP2AFBtau::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_F(wf_F, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_F.SetAbsTolerance(1.E-17); // desired absolute error
            ig_F.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_F = ig_F.Integral(0.0, 1.0); // interval
            //std::cout << sigma_box_F << " +- " << ig_F.Error()  << std::endl;
            // results: 6.0e-12 -- 2.0e-10
            //
            ROOT::Math::Functor1D wf_B(this, &LEP2AFBtau::Integrand_dsigmaBox_l);
            ROOT::Math::Integrator ig_B(wf_B, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_B.SetAbsTolerance(1.E-18); // desired absolute error
            ig_B.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_B = ig_B.Integral(-1.0, 0.0); // interval
            //std::cout << sigma_box_B << " +- " << ig_B.Error()  << std::endl;
            // results: 3.4e-13 -- 8.7e-11
            
            // denominator
            if (!flag[ISR]) {
                myLEP2sigmaTau.setFlags(flag);
                sigma = myLEP2sigmaTau.computeThValue()/SM.GeVminus2_to_nb/1000.0;
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    }
    double AFB_tau = SMresult_cache;
    
    #ifdef LEP2TEST
    AFB_tau = myTEST.AFBtauTEST(sqrt_s);
    #endif
            
    if ( checkLEP2NP() ) {
        double obliqueShat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueShat();
        double obliqueThat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueThat();
        double obliqueUhat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueUhat();
        double obliqueV = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueV();
        double obliqueW = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueW();
        double obliqueX = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueX();
        double obliqueY = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueY();
        double ObParam[7] = {obliqueShat, obliqueThat, obliqueUhat,
                             obliqueV, obliqueW, obliqueX, obliqueY};
        AFB_tau += myLEP2oblique.AFB_l_LEP2_NP(SM.TAU, s, ml_cache, ObParam);
    }
    
    if ( checkLEP2GIMR() ) {
        double deltaAFB = 0., sigma = 0.0,dFminusB =0.0;
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
        
        myLEP2sigmaTau.setFlags(flag);
        sigma = myLEP2sigmaTau.computeThValue()/SM.GeVminus2_to_nb/1000.0;
        dFminusB = (myLEP2GIMR.sigmaF_l_LEP2_GIMR(StandardModel::TAU, s, GIMRParam)
                    -myLEP2GIMR.sigmaB_l_LEP2_GIMR(StandardModel::TAU, s, GIMRParam));
        
        deltaAFB = AFB_tau*(1-myLEP2GIMR.sigma_l_LEP2_GIMR(StandardModel::TAU, s, GIMRParam)/sigma)+dFminusB/sigma;
        
        AFB_tau += deltaAFB;
    }

    return AFB_tau;
}
        
