/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2Rcharm.h"


double LEP2Rcharm::computeThValue() 
{ 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mq_cache = m_q(SM.CHARM, sqrt_s);
        
        myLEP2sigmaCharm.setFlags(flag);
        double sigma_c = myLEP2sigmaCharm.computeThValue();
        myLEP2sigmaHadron.setFlags(flag);        
        double sigma_had = myLEP2sigmaHadron.computeThValue();

        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    #ifdef LEP2TEST
    R_c = myTEST.RcharmTEST(sqrt_s);
    #endif
    
    if ( checkLEP2NP() ) {
        double obliqueShat = (static_cast<const NPbase*> (&SM))->obliqueShat();
        double obliqueThat = (static_cast<const NPbase*> (&SM))->obliqueThat();
        double obliqueUhat = (static_cast<const NPbase*> (&SM))->obliqueUhat();
        double obliqueV = (static_cast<const NPbase*> (&SM))->obliqueV();
        double obliqueW = (static_cast<const NPbase*> (&SM))->obliqueW();
        double obliqueX = (static_cast<const NPbase*> (&SM))->obliqueX();
        double obliqueY = (static_cast<const NPbase*> (&SM))->obliqueY();
        double ObParam[7] = {obliqueShat, obliqueThat, obliqueUhat,
                             obliqueV, obliqueW, obliqueX, obliqueY};
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
    }
    
    return R_c;
}


