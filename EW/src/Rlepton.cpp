/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rlepton.h"


double Rlepton::getThValue() 
{
    double R0_l;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        R0_l = myEW.getMyEW_CHMN().R_l(SM.ELECTRON);
    else if (myEWTYPE==EW::EWABC) 
        R0_l = myEW.getMyEW_ABC().R_l(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double R_0 = 20.8228*(1.0 + 1.05*delta_als - 0.28*delta_alpha);
        R0_l = R_0*(1.0 + 0.28*SM.epsilon1() - 0.36*SM.epsilon3() + 0.50*SM.epsilonb());
    } else {       
        R0_l = myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
        
        /* Oblique NP */
        if ( myEW.checkSTU() && !SM.IsFlagNotLinearizedNP() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //R0_l = 20.740;

                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                double delta_l = - 0.000192*myEW.S() + 0.000790*myEW.T();
                R0_l *= 1.0 + delta_had/myEW.Gamma_had()
                        - delta_l/myEW.Gamma_l(SM.ELECTRON);
            } else {
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                
                R0_l += 8.0*alpha*(3.0-2.0*s2)*(1.0-5.0*s2)
                        /3.0/pow(1.0-4.0*s2+8.0*s4, 2.0)/(c2-s2)
                        *( myEW.S() - 4.0*c2*s2*myEW.T() );                
            }
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            bool nonZeroNP = false;

            double delGVe = SM.deltaGVl(SM.ELECTRON);
            double delGAe = SM.deltaGAl(SM.ELECTRON);
            if (delGVe!=0.0 || delGAe!=0.0) nonZeroNP = true;

            double delGVq[6], delGAq[6];
            for (int p=0; p<6; ++p) {
                delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
                delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
                if (delGVq[p]!=0.0 || delGAq[p]!=0.0) nonZeroNP = true;
            }

            if (nonZeroNP) {
                double gVe = SM.StandardModel::gVl(SM.ELECTRON).real();
                double gAe = SM.StandardModel::gAl(SM.ELECTRON).real();
                double Ge = gVe*gVe + gAe*gAe;
                double deltaGe = 2.0*(gVe*delGVe + gAe*delGAe);

                double Gq[6], deltaGq[6];
                double gVq, gAq;
                double Gq_sum = 0.0, delGq_sum = 0.0;
                for (int p=0; p<6; ++p) {
                    gVq = SM.StandardModel::gVq((StandardModel::quark)p).real();
                    gAq = SM.StandardModel::gAq((StandardModel::quark)p).real();
                    Gq[p] = gVq*gVq + gAq*gAq;
                    deltaGq[p] = 2.0*(gVq*delGVq[p] + gAq*delGAq[p]);

                    Gq_sum += 3.0*Gq[p];
                    delGq_sum += 3.0*deltaGq[p];
                }

                R0_l += delGq_sum/Ge - Gq_sum*deltaGe/Ge/Ge;
            }
        }
        
        /* TEST */
        //R0_l -= myEW.Gamma_had()/myEW.Gamma_l(SM.ELECTRON);
    }
 
    return R0_l;
}
        

