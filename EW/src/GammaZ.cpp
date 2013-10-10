/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GammaZ.h"


double GammaZ::getThValue() 
{ 
    double Gamma_Z;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        Gamma_Z = myEW.getMyEW_CHMN().GammaZ();
    else if (myEWTYPE==EW::EWABC) 
        Gamma_Z = myEW.getMyEW_ABC().GammaZ(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double Gamma_T0 = 2.48946*(1.0 + 0.73*delta_als - 0.35*delta_alpha);
        Gamma_Z = Gamma_T0*(1.0 + 1.35*SM.epsilon1() - 0.46*SM.epsilon3() + 0.35*SM.epsilonb());
    } else {
        Gamma_Z = myEW.Gamma_Z();

        /* Theoretical uncertainty */
        Gamma_Z += SM.getDelGammaZ();

        if(myEWTYPE==EW::EWBURGESS) {
            Gamma_Z += - 0.00961*SM.obliqueS() + 0.0263*SM.obliqueT();
            return Gamma_Z;
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            bool nonZeroNP = false;

            double delGVl[6], delGAl[6], delGVq[6], delGAq[6];
            for (int p=0; p<6; ++p) {
                delGVl[p] = SM.deltaGVl((StandardModel::lepton)p);
                delGAl[p] = SM.deltaGAl((StandardModel::lepton)p);
                delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
                delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
                if (delGVl[p]!=0.0 || delGAl[p]!=0.0
                        || delGVq[p]!=0.0 || delGAq[p]!=0.0)
                    nonZeroNP = true;
            }

            if (nonZeroNP) {
                double gVf, gAf;
                double deltaGl[6], deltaGq[6];
                double delGammaZ = 0.0;
                for (int p=0; p<6; ++p) {
                    gVf = SM.StandardModel::gVl((StandardModel::lepton)p).real();
                    gAf = SM.StandardModel::gAl((StandardModel::lepton)p).real();
                    deltaGl[p] = 2.0*(gVf*delGVl[p] + gAf*delGAl[p]);

                    gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
                    gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
                    deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

                    delGammaZ += deltaGl[p] + 3.0*deltaGq[p];
                }

                Gamma_Z += SM.alphaMz()*SM.getMz()/12.0/myEW.sW2_SM()/myEW.cW2_SM()
                           * delGammaZ;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("GammaZ::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //Gamma_Z -= myEW.Gamma_Z();
    }
      
    return Gamma_Z;
}
        
