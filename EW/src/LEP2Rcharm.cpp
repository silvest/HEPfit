/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2Rcharm.h"


double LEP2Rcharm::getThValue() { 
    Mw = SM.Mw(); 
    GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        mq_cache = m_q(SM.CHARM, sqrt_s);
        
        myLEP2sigmaCharm.setFlags(flag);
        double sigma_c = myLEP2sigmaCharm.getThValue();
        myLEP2sigmaHadron.setFlags(flag);        
        double sigma_had = myLEP2sigmaHadron.getThValue();

        SMresult_cache = sigma_c/sigma_had;

        if ( myEW.checkModelForSTU() && SM.IsFlagFixedAllSMparams() ) {
            double ObParam[7];
            for (int i=0; i<7; i++) {
                SetObParam((LEP2oblique::Oblique)i, ObParam);
                Coeff_cache[i] 
                    = myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
            }
        }
    }
    double R_c = SMresult_cache;
    
    #ifdef LEP2TEST
    R_c = myTEST.RcharmTEST(sqrt_s);
    #endif
    
    if ( myEW.checkModelForSTU() ) {
        if ( SM.IsFlagFixedAllSMparams() ) {
            R_c += Coeff_cache[myLEP2oblique.Shat]*myEW.Shat()
                 + Coeff_cache[myLEP2oblique.That]*myEW.That()
                 + Coeff_cache[myLEP2oblique.Uhat]*myEW.Uhat()
                 + Coeff_cache[myLEP2oblique.V]*myEW.V()
                 + Coeff_cache[myLEP2oblique.W]*myEW.W()
                 + Coeff_cache[myLEP2oblique.X]*myEW.X()
                 + Coeff_cache[myLEP2oblique.Y]*myEW.Y();                 
        } else {
            double ObParam[7] = {myEW.Shat(), myEW.That(), myEW.Uhat(),
                                 myEW.V(), myEW.W(), myEW.X(), myEW.Y()};
            R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s, mq_cache, ObParam);
        }
    }
    
    return R_c;
}


