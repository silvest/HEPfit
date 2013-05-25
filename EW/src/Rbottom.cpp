/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rbottom.h"
#include "sigmaHadron.h"
#include <EWSM.h>


double Rbottom::getThValue() 
{ 
    double R0_b;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        R0_b = myEW.getMyEW_CHMN().R_b();
    else if (myEWTYPE==EW::EWABC) 
        R0_b = myEW.getMyEW_ABC().R_b(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else if (myEWTYPE==EW::EWABC2) {
        double R_b0 = 0.2182355;
        R0_b = R_b0*(1.0 - 0.06*SM.epsilon1() + 0.07*SM.epsilon3() + 1.79*SM.epsilonb());
    } else {    
        if (SM.IsFlagApproximateGqOverGb() 
                && !SM.IsFlagRhoZbFromGuOverGb()
                && !SM.IsFlagRhoZbFromGdOverGb()
                && !SM.IsFlagTestSubleadingTwoLoopEW()) {
            /* We use this part in the case where rhoZb is not derived from 
             * the approximate formula of either Gu/Gb or Gd/Gb, or where 
             * it is not calculated from the input delRhoZb. */
            double Gu_over_Gb = SM.getEWSM()->Gu_over_Gb_SM();
            double Gd_over_Gb = SM.getEWSM()->Gd_over_Gb_SM();
            R0_b = 1.0/(1.0 + 2.0*(Gd_over_Gb + Gu_over_Gb));
        } else
            R0_b = myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();

        if(myEWTYPE==EW::EWBURGESS) {
            double delta_b = - 0.00171*SM.obliqueS() + 0.00416*SM.obliqueT();
            double delta_had = - 0.00901*SM.obliqueS() + 0.0200*SM.obliqueT();
            R0_b *= 1.0 + delta_b/myEW.Gamma_q(SM.BOTTOM)
                    - delta_had/myEW.Gamma_had();
            return R0_b;
        }

        /* NP contribution to the Zff vertex */
        if ( !SM.IsFlagNotLinearizedNP() ) {
            bool nonZeroNP = false;
            double delGVq[6], delGAq[6];
            for (int p=0; p<6; ++p) {
                delGVq[p] = SM.deltaGVq((StandardModel::quark)p);
                delGAq[p] = SM.deltaGAq((StandardModel::quark)p);
                if (delGVq[p]!=0.0 || delGAq[p]!=0.0) nonZeroNP = true;
            }

            if (nonZeroNP) {
                double gVf, gAf;
                double Gq[6], deltaGq[6];
                double Gq_sum = 0.0, delGq_sum = 0.0;
                for (int p=0; p<6; ++p) {
                    gVf = SM.StandardModel::gVq((StandardModel::quark)p).real();
                    gAf = SM.StandardModel::gAq((StandardModel::quark)p).real();
                    Gq[p] = gVf*gVf + gAf*gAf;
                    deltaGq[p] = 2.0*(gVf*delGVq[p] + gAf*delGAq[p]);

                    Gq_sum += Gq[p]; /* without the color factor */
                    delGq_sum += deltaGq[p]; /* without the color factor */
                }

                R0_b += deltaGq[(int)SM.BOTTOM]/Gq_sum
                        - Gq[(int)SM.BOTTOM]*delGq_sum/Gq_sum/Gq_sum;
            }
        } else
            if (SM.obliqueS()!=0.0 || SM.obliqueT()!=0.0 || SM.obliqueU()!=0.0)
                throw std::runtime_error("Rbottom::getThValue(): The oblique corrections STU cannot be used with flag NotLinearizedNP=1");

        /* Debug: extract pure NP contribution */
        //R0_b -= myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();
    }
    
    return R0_b;
}
        

