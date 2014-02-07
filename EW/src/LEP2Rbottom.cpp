/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
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
    GammaZ = myEW.Gamma_Z();

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
        double obliqueShat = (static_cast<const NPbase*> (&SM))->obliqueShat();
        double obliqueThat = (static_cast<const NPbase*> (&SM))->obliqueThat();
        double obliqueUhat = (static_cast<const NPbase*> (&SM))->obliqueUhat();
        double obliqueV = (static_cast<const NPbase*> (&SM))->obliqueV();
        double obliqueW = (static_cast<const NPbase*> (&SM))->obliqueW();
        double obliqueX = (static_cast<const NPbase*> (&SM))->obliqueX();
        double obliqueY = (static_cast<const NPbase*> (&SM))->obliqueY();
        double ObParam[7] = {obliqueShat, obliqueThat, obliqueUhat,
                             obliqueV, obliqueW, obliqueX, obliqueY};
        R_b += myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
    }
    
    return R_b;
}
        




