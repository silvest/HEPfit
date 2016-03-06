/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2Rbottom.h"


double LEP2Rbottom::computeThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = SM.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mq_cache = m_q(SM.BOTTOM, sqrt_s);
        
        myLEP2sigmaBottom.setFlags(flag);
        double sigma_b = myLEP2sigmaBottom.computeThValue();
        myLEP2sigmaHadron.setFlags(flag);
        double sigma_had = myLEP2sigmaHadron.computeThValue();

        SMresult_cache = sigma_b/sigma_had;
    }
    double R_b = SMresult_cache;
    
    #ifdef LEP2TEST
    R_b = myTEST.RbottomTEST(sqrt_s);
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
        R_b += myLEP2oblique.R_q_LEP2_NP(QCD::BOTTOM, s, mq_cache, ObParam);
    }
    
    if ( checkLEP2GIMR() ) {
        double deltaR_b = 0., dsigma_b = 0.0, dsigma_had =0.0;
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
        
        myLEP2sigmaHadron.setFlags(flag);
        double sigma_had = myLEP2sigmaHadron.computeThValue();
        
        dsigma_had = myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::UP, s, GIMRParamUP)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::DOWN, s, GIMRParamDOWN)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::CHARM, s, GIMRParamCHARM)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::STRANGE, s, GIMRParamSTRANGE)+
        myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamBOTTOM);
        
        dsigma_b = myLEP2GIMR.sigma_q_LEP2_GIMR(QCD::BOTTOM, s, GIMRParamBOTTOM);
        
        deltaR_b = R_b*(1-dsigma_had/sigma_had)+dsigma_b/sigma_had;
        
        R_b += deltaR_b;
    }

    return R_b;
}





