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

        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //Gamma_Z = 2.4954; 

                Gamma_Z += - 0.00961*myEW.S() + 0.0263*myEW.T();
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                
                Gamma_Z += alpha*alpha*SM.getMz()/72.0/c2/s2/(c2-s2)
                           *( -10.0*(3.0-8.0*s2)*myEW.S() 
                              + (63.0-126.0*s2-40.0*s4)*myEW.T() );
            }
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
        }
        
        /* TEST */
        //Gamma_Z -=  myEW.Gamma_Z();
    }
      
    return Gamma_Z;
}
        
