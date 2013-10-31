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

        if ( myEW.checkSTUVWXY() && SM.IsFlagFixedAllSMparams() ) {
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
            }
        }
    }
    double R_b = SMresult_cache;
    
    #ifdef LEP2TEST
    R_b = myTEST.RbottomTEST(sqrt_s);
    #endif
            
    if ( myEW.checkSTUVWXY() ) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            R_b += Coeff_cache[myLEP2oblique.Shat]*SM.obliqueShat()
                 + Coeff_cache[myLEP2oblique.That]*SM.obliqueThat()
                 + Coeff_cache[myLEP2oblique.Uhat]*SM.obliqueUhat()
                 + Coeff_cache[myLEP2oblique.V]*SM.obliqueV()
                 + Coeff_cache[myLEP2oblique.W]*SM.obliqueW()
                 + Coeff_cache[myLEP2oblique.X]*SM.obliqueX()
                 + Coeff_cache[myLEP2oblique.Y]*SM.obliqueY();
        } else {
            double ObParam[7] = {SM.obliqueShat(), SM.obliqueThat(), SM.obliqueUhat(),
                                 SM.obliqueV(), SM.obliqueW(), SM.obliqueX(), SM.obliqueY()};
            R_b += myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s, mq_cache, ObParam);
        }
    }
    
    return R_b;
}
        




