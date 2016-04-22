/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2sigmaHadron.h"


double LEP2sigmaHadron::computeThValue()
{ 
    Mw = SM.Mw(); 
    GammaZ = SM.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mqForHad_cache[SM.UP] = m_q(SM.UP, sqrt_s);
        mqForHad_cache[SM.DOWN] = m_q(SM.DOWN, sqrt_s);
        mqForHad_cache[SM.CHARM] = m_q(SM.CHARM, sqrt_s);
        mqForHad_cache[SM.STRANGE] = m_q(SM.STRANGE, sqrt_s);
        mqForHad_cache[SM.BOTTOM] = m_q(SM.BOTTOM, sqrt_s);

        if (!flag[ISR]) {
            q_flavor = QCD::UP;
            mq_cache = mqForHad_cache[SM.UP];
            SMresult_cache = sigma_NoISR_q();
            //
            q_flavor = QCD::DOWN;
            mq_cache = mqForHad_cache[SM.DOWN];
            SMresult_cache += sigma_NoISR_q();
            //
            q_flavor = QCD::CHARM;
            mq_cache = mqForHad_cache[SM.CHARM];
            SMresult_cache += sigma_NoISR_q();
            //
            q_flavor = QCD::STRANGE;
            mq_cache = mqForHad_cache[SM.STRANGE];
            SMresult_cache += sigma_NoISR_q();
            //
            q_flavor = QCD::BOTTOM;
            mq_cache = mqForHad_cache[SM.BOTTOM];
            SMresult_cache += sigma_NoISR_q();
        } else {
            q_flavor = QCD::UP;
            mq_cache = mqForHad_cache[SM.UP];
            ROOT::Math::Functor1D wf_UP(this, &LEP2sigmaHadron::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig_UP(wf_UP, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_UP.SetAbsTolerance(1.E-13); // desired absolute error
            ig_UP.SetRelTolerance(1.E-4); // desired relative error
            double sigma_UP = ig_UP.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "UP " << sigma_UP << " +- " << ig_UP.Error() << std::endl;
            // results: 1.2e-8 -- 4.7e-8
            //
            q_flavor = QCD::DOWN;
            mq_cache = mqForHad_cache[SM.DOWN];
            ROOT::Math::Functor1D wf_DOWN(this, &LEP2sigmaHadron::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig_DOWN(wf_DOWN, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_DOWN.SetAbsTolerance(1.E-14); // desired absolute error
            ig_DOWN.SetRelTolerance(1.E-4); // desired relative error
            double sigma_DOWN = ig_DOWN.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "DOWN " << sigma_DOWN << " +- " << ig_DOWN.Error() << std::endl;
            // results: 7.1e-9 -- 4.0e-8
            //
            q_flavor = QCD::CHARM;
            mq_cache = mqForHad_cache[SM.CHARM];
            ROOT::Math::Functor1D wf_CHARM(this, &LEP2sigmaHadron::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig_CHARM(wf_CHARM, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_CHARM.SetAbsTolerance(1.E-13); // desired absolute error
            ig_CHARM.SetRelTolerance(1.E-4); // desired relative error
            double sigma_CHARM = ig_CHARM.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "CHARM " << sigma_CHARM << " +- " << ig_CHARM.Error() << std::endl;
            // results: 1.2e-8 -- 4.7e-8
            //
            q_flavor = QCD::STRANGE;
            mq_cache = mqForHad_cache[SM.STRANGE];
            ROOT::Math::Functor1D wf_STRANGE(this, &LEP2sigmaHadron::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig_STRANGE(wf_STRANGE, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_STRANGE.SetAbsTolerance(1.E-14); // desired absolute error
            ig_STRANGE.SetRelTolerance(1.E-4); // desired relative error
            double sigma_STRANGE = ig_STRANGE.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "STRANGE " << sigma_STRANGE << " +- " << ig_STRANGE.Error() << std::endl;
            // results: 7.1e-9 -- 4.0e-8
            //
            q_flavor = QCD::BOTTOM;
            mq_cache = mqForHad_cache[SM.BOTTOM];
            ROOT::Math::Functor1D wf_BOTTOM(this, &LEP2sigmaHadron::Integrand_sigmaWithISR_q);
            ROOT::Math::Integrator ig_BOTTOM(wf_BOTTOM, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_BOTTOM.SetAbsTolerance(1.E-14); // desired absolute error
            ig_BOTTOM.SetRelTolerance(1.E-4); // desired relative error
            double sigma_BOTTOM = ig_BOTTOM.Integral(0.0, 1.0-0.85*0.85); // interval
            //std::cout << "BOTTOM " << sigma_BOTTOM << " +- " << ig_BOTTOM.Error() << std::endl;
            // results: 7.0e-9 -- 4.0e-8
            //
            SMresult_cache = sigma_UP + sigma_DOWN + sigma_CHARM 
                             + sigma_STRANGE + sigma_BOTTOM;
        }
        
        if (flag[WeakBox]) {
            q_flavor = QCD::UP;
            mq_cache = mqForHad_cache[SM.UP];
            ROOT::Math::Functor1D wf_box_UP(this, &LEP2sigmaHadron::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box_UP(wf_box_UP, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box_UP.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box_UP.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_UP = ig_box_UP.Integral(-1.0, 1.0); // interval
            //std::cout << "UP_box " << sigma_box_UP << " +- " << ig_box_UP.Error() << std::endl;
            // results: 6.2e-11 -- 4.2e-10
            //
            q_flavor = QCD::DOWN;
            mq_cache = mqForHad_cache[SM.DOWN];
            ROOT::Math::Functor1D wf_box_DOWN(this, &LEP2sigmaHadron::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box_DOWN(wf_box_DOWN, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box_DOWN.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box_DOWN.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_DOWN = ig_box_DOWN.Integral(-1.0, 1.0); // interval
            //std::cout << "DOWN_box " << sigma_box_DOWN << " +- " << ig_box_DOWN.Error() << std::endl;
            // results: 1.4e-11 -- 6.8e-10
            //
            q_flavor = QCD::CHARM;
            mq_cache = mqForHad_cache[SM.CHARM];
            ROOT::Math::Functor1D wf_box_CHARM(this, &LEP2sigmaHadron::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box_CHARM(wf_box_CHARM, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box_CHARM.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box_CHARM.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_CHARM = ig_box_CHARM.Integral(-1.0, 1.0); // interval
            //std::cout << "CHARM_box " << sigma_box_CHARM << " +- " << ig_box_CHARM.Error() << std::endl;
            // results: 6.2e-11 -- 4.2e-10
            //
            q_flavor = QCD::STRANGE;
            mq_cache = mqForHad_cache[SM.STRANGE];
            ROOT::Math::Functor1D wf_box_STRANGE(this, &LEP2sigmaHadron::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box_STRANGE(wf_box_STRANGE, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box_STRANGE.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box_STRANGE.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_STRANGE = ig_box_STRANGE.Integral(-1.0, 1.0); // interval
            //std::cout << "STRANGE_box " << sigma_box_STRANGE << " +- " << ig_box_STRANGE.Error() << std::endl;
            // results: 1.4e-11 -- 6.8e-10
            //
            q_flavor = QCD::BOTTOM;
            mq_cache = mqForHad_cache[SM.BOTTOM];
            ROOT::Math::Functor1D wf_box_BOTTOM(this, &LEP2sigmaHadron::Integrand_dsigmaBox_q);
            ROOT::Math::Integrator ig_box_BOTTOM(wf_box_BOTTOM, ROOT::Math::IntegrationOneDim::kADAPTIVESINGULAR);
            ig_box_BOTTOM.SetAbsTolerance(1.E-16); // desired absolute error
            ig_box_BOTTOM.SetRelTolerance(1.E-4); // desired relative error
            double sigma_box_BOTTOM = ig_box_BOTTOM.Integral(-1.0, 1.0); // interval
            //std::cout << "BOTTOM_box " << sigma_box_BOTTOM << " +- " << ig_box_BOTTOM.Error() << std::endl;
            // results: 1.4e-11 -- 2.5e-10
            //
            SMresult_cache += sigma_box_UP + sigma_box_DOWN + sigma_box_CHARM 
                              + sigma_box_STRANGE + sigma_box_BOTTOM;
        }        
    }
    double sigmaH = SMresult_cache;
    
    #ifdef LEP2TEST
    sigmaH = myTEST.sigmaHadronTEST(sqrt_s)/SM.GeVminus2_to_nb/1000.0;
    #endif
    
    if ( checkLEP2NP() && !bSigmaForR) {
        double obliqueShat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueShat();
        double obliqueThat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueThat();
        double obliqueUhat = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueUhat();
        double obliqueV = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueV();
        double obliqueW = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueW();
        double obliqueX = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueX();
        double obliqueY = (dynamic_cast<const NPSTUVWXY*> (&SM))->obliqueY();
        double ObParam[7] = {obliqueShat, obliqueThat, obliqueUhat,
                             obliqueV, obliqueW, obliqueX, obliqueY};
        sigmaH += myLEP2oblique.sigma_q_LEP2_NP(QCD::UP, s, mqForHad_cache[SM.UP], ObParam)
                  + myLEP2oblique.sigma_q_LEP2_NP(QCD::DOWN, s, mqForHad_cache[SM.DOWN], ObParam)
                  + myLEP2oblique.sigma_q_LEP2_NP(QCD::CHARM, s, mqForHad_cache[SM.CHARM], ObParam) 
                  + myLEP2oblique.sigma_q_LEP2_NP(QCD::STRANGE, s, mqForHad_cache[SM.STRANGE], ObParam) 
                  + myLEP2oblique.sigma_q_LEP2_NP(QCD::BOTTOM, s, mqForHad_cache[SM.BOTTOM], ObParam);
    }
    
    if ( checkLEP2GIMR() && !bSigmaForAFB ) {
        
        
        double deltaGammaZ = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGamma_Z();
        
        double deltaGR_u = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getQuarks(SM.UP));
        double deltaGL_u = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getQuarks(SM.UP));
        double myCLL_u = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLL_up();
        double myCLR_u = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLR_up();
        double myCRL_u = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRL_up();
        double myCRR_u = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRR_up();
        
        double deltaGR_d = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getQuarks(SM.DOWN));
        double deltaGL_d = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getQuarks(SM.DOWN));
        double myCLL_d = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLL_down();
        double myCLR_d = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLR_down();
        double myCRL_d = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRL_down();
        double myCRR_d = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRR_down();
        
        double deltaGR_s = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getQuarks(SM.STRANGE));
        double deltaGL_s = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getQuarks(SM.STRANGE));
        double myCLL_s = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLL_strange();
        double myCLR_s = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLR_strange();
        double myCRL_s = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRL_strange();
        double myCRR_s = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRR_strange();
        
        double deltaGR_c = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getQuarks(SM.CHARM));
        double deltaGL_c = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getQuarks(SM.CHARM));
        double myCLL_c = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLL_charm();
        double myCLR_c = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLR_charm();
        double myCRL_c = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRL_charm();
        double myCRR_c = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRR_charm();
        
        double deltaGR_b = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getQuarks(SM.BOTTOM));
        double deltaGL_b = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getQuarks(SM.BOTTOM));
        double myCLL_b = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLL_bottom();
        double myCLR_b = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CLR_bottom();
        double myCRL_b = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRL_bottom();
        double myCRR_b = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->CRR_bottom();
        
        double deltaGR_e = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGR_f(SM.getLeptons(StandardModel::ELECTRON));
        double deltaGL_e = (dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaGL_f(SM.getLeptons(StandardModel::ELECTRON));
        double deltaMz2=(dynamic_cast<const NPEffectiveGIMR*> (&SM))->deltaMz2();
        
        double GIMRParamUP[10] ={myCLL_u,myCLR_u,myCRL_u,myCRR_u,deltaGammaZ,deltaGL_u,deltaGR_u,deltaGL_e,deltaGR_e,deltaMz2};
        double GIMRParamDOWN[10] ={myCLL_d,myCLR_d,myCRL_d,myCRR_d,deltaGammaZ,deltaGL_d,deltaGR_d,deltaGL_e,deltaGR_e,deltaMz2};
        double GIMRParamCHARM[10] ={myCLL_c,myCLR_c,myCRL_c,myCRR_c,deltaGammaZ,deltaGL_c,deltaGR_c,deltaGL_e,deltaGR_e,deltaMz2};
        double GIMRParamSTRANGE[10] ={myCLL_s,myCLR_s,myCRL_s,myCRR_s,deltaGammaZ,deltaGL_s,deltaGR_s,deltaGL_e,deltaGR_e,deltaMz2};
        double GIMRParamBOTTOM[10] ={myCLL_b,myCLR_b,myCRL_b,myCRR_b,deltaGammaZ,deltaGL_b,deltaGR_b,deltaGL_e,deltaGR_e,deltaMz2};
        
        sigmaH += myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::UP, s, GIMRParamUP)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::DOWN, s, GIMRParamDOWN)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::CHARM, s, GIMRParamCHARM)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::STRANGE, s, GIMRParamSTRANGE)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamBOTTOM);
    }

    return ( sigmaH*SM.GeVminus2_to_nb*1000.0 );
}
        

