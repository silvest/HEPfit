/* 
 * File:   LEP2Rbottom.cpp
 * Author: mishima
 */

#include <Math/Functor.h>
#include <Math/Integrator.h>
#include <Math/AllIntegrationTypes.h>
#include "LEP2Rbottom.h"


double LEP2Rbottom::getThValue() { 
    double Mw = SM.Mw(); 
    double GammaZ = myEW.Gamma_Z();

    if (!checkSMparams(s, Mw, GammaZ)) {
        myLEP2sigmaBottom.setFlags(flag);
        double sigma_b = myLEP2sigmaBottom.getThValue();
        myLEP2sigmaHadron.setFlags(flag);
        double sigma_had = myLEP2sigmaHadron.getThValue();

        SMresult_cache = sigma_b/sigma_had;
    }
    double R_b = SMresult_cache;
    
    if ( myEW.checkModelForSTU() ) {
        R_b += myLEP2oblique.R_q_LEP2_NP(StandardModel::BOTTOM, s);
    }
    
    return R_b;
}
        




