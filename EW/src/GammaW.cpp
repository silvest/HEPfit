/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaW.h"


double GammaW::getThValue() 
{  
    double Gamma_W;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        Gamma_W = myEW.getMyEW_CHMN().GammaW();
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        throw std::runtime_error("GammaW::getThValue() is not implemented for EW::EWABC");  
    else {
        Gamma_W = SM.GammaW();
        
        double Wbar = (SM.obliqueV() - SM.obliqueW())/SM.alphaMz();

        if(myEWTYPE==EW::EWBURGESS) {
            Gamma_W *= 1.0 - 0.00723*SM.obliqueS() + 0.0111*SM.obliqueT()
                       + 0.00849*SM.obliqueU() + 0.00781*Wbar;
            return Gamma_W;
        }

        if (!SM.IsFlagNotLinearizedNP() ) {
            double alpha = myEW.alpha();
            double c2 = myEW.cW2_SM();
            double s2 = myEW.sW2_SM();
                
            //Gamma_W *= 1.0 - alpha/2.0/(c2-s2) /* corrected on May 27, 2013 */
            Gamma_W *= 1.0 - 3.0*alpha/4.0/(c2-s2)
                       *( SM.obliqueS() - 2.0*c2*SM.obliqueT()
                          - (c2-s2)*SM.obliqueU()/2.0/s2 - 2.0*(c2 - s2)*Wbar )
                        //- s2/(c2-s2)*SM.DeltaGF(); /* corrected on May 27, 2013 */
                        - 3.0*s2/2.0/(c2-s2)*SM.DeltaGF();
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("GammaW::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //Gamma_W -= SM.GammaW();
    }
 
    return Gamma_W;
}
