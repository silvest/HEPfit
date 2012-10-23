/* 
 * File:   LEP2Rcharm.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2Rcharm.h"


double LEP2Rcharm::getThValue() { 
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        myLEP2sigmaCharm.setFlags(flag);
        double sigma_c = myLEP2sigmaCharm.getThValue();
        myLEP2sigmaHadron.setFlags(flag);        
        double sigma_had = myLEP2sigmaHadron.getThValue();

        SMresult_cache = sigma_c/sigma_had;
    }
    double R_c = SMresult_cache;
    
    if ( myEW.checkModelForSTU() )
        R_c += myLEP2oblique.R_q_LEP2_NP(StandardModel::CHARM, s);
    
    return R_c;
}


