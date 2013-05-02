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
        
        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //R0_b = 0.21474;
                
                double delta_b = - 0.00171*myEW.S() + 0.00416*myEW.T();
                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                R0_b *= 1.0 + delta_b/myEW.Gamma_q(SM.BOTTOM) 
                        - delta_had/myEW.Gamma_had();
            } else {
                double alpha = myEW.alpha0();  
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                R0_b += 6.0*alpha*(9.0-36.0*s2+16.0*s4)
                        /pow(45.0-84.0*s2+88.0*s4, 2.0)/(c2-s2)
                        *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }
        
        if (SM.IsFlagNPZbbbarLinearize() && (SM.deltaGVb()!=0.0 || SM.deltaGAb()!=0.0) ) {
            double gVb0 = SM.getQuarks(SM.BOTTOM).getIsospin() 
                          - 2.0*SM.getQuarks(SM.BOTTOM).getCharge()*myEW.s02();
            double gAb0 = SM.getQuarks(SM.BOTTOM).getIsospin();        
            double gVu0 = SM.getQuarks(SM.UP).getIsospin() 
                          - 2.0*SM.getQuarks(SM.UP).getCharge()*myEW.s02();
            double gAu0 = SM.getQuarks(SM.UP).getIsospin();        
            double Nc = 3.0;
            double sum = Nc*2.0*(gVu0*gVu0 + gAu0*gAu0)
                         + Nc*3.0*(gVb0*gVb0 + gAb0*gAb0);
            double R0b0 = Nc*(gVb0*gVb0 + gAb0*gAb0)/sum;
            double coeff = 2.0*Nc*(1.0 - R0b0)/sum;
            double coeffV = coeff*gVb0;
            double coeffA = coeff*gAb0;
            //std::cout << "cV: " << coeffV << std::endl;
            //std::cout << "cA: " << coeffA << std::endl;
            //std::cout << "cL: " << coeffV+coeffA << std::endl;
            //std::cout << "cR: " << coeffV-coeffA << std::endl;

            R0_b += coeffV*SM.deltaGVb() + coeffA*SM.deltaGAb();
        }        
        
        /* TEST */
        //R0_b -= myEW.Gamma_q(SM.BOTTOM)/myEW.Gamma_had();
    }
    
    return R0_b;
}
        

