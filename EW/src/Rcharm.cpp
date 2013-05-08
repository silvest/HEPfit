/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Rcharm.h"
#include <EWSM.h>


double Rcharm::getThValue() 
{   
    double R0_c;
    EW::EWTYPE myEWTYPE = myEW.getEWTYPE();

    if (myEWTYPE==EW::EWCHMN)  
        R0_c = myEW.getMyEW_CHMN().R_c();
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        R0_c = myEW.getMyEW_ABC().R_c(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {    
        if (SM.IsFlagApproximateGqOverGb() 
                && !SM.IsFlagRhoZbFromGuOverGb()
                && !SM.IsFlagRhoZbFromGdOverGb()
                && !SM.IsFlagTestSubleadingTwoLoopEW()) {
            double Gu_over_Gb = SM.getEWSM()->Gu_over_Gb_SM();
            double Gd_over_Gb = SM.getEWSM()->Gd_over_Gb_SM();
            R0_c = Gu_over_Gb/(1.0 + 2.0*(Gd_over_Gb + Gu_over_Gb));
        } else
            R0_c = myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();
        
        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //R0_c = 0.17223;
                    
                double delta_c_over_Gamma_c = - 0.00649*myEW.S() + 0.0124*myEW.T();
                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                R0_c *= 1.0 + delta_c_over_Gamma_c - delta_had/myEW.Gamma_had();
            } else {            
                double alpha = myEW.alpha();  
                double c2 = myEW.cW2_SM();
                double s2 = myEW.sW2_SM();
                double s4 = s2*s2;
                
                R0_c -= 9.0*alpha*(9.0-36.0*s2+16.0*s4)
                        /pow(45.0-84.0*s2+88.0*s4, 2.0)/(c2-s2)
                        *( myEW.S() - 4.0*c2*s2*myEW.T() );
            }
        }

        if (SM.IsFlagNPZbbbarLinearize() && (SM.deltaGVb()!=0.0 || SM.deltaGAb()!=0.0) ) {
            double gVb0 = SM.getQuarks(SM.BOTTOM).getIsospin() 
                          - 2.0*SM.getQuarks(SM.BOTTOM).getCharge()*myEW.sW2_SM();
            double gAb0 = SM.getQuarks(SM.BOTTOM).getIsospin();        
            double gVc0 = SM.getQuarks(SM.CHARM).getIsospin() 
                          - 2.0*SM.getQuarks(SM.CHARM).getCharge()*myEW.sW2_SM();
            double gAc0 = SM.getQuarks(SM.CHARM).getIsospin();        
            double Nc = 3.0;
            double sum = Nc*2.0*(gVc0*gVc0 + gAc0*gAc0)
                         + Nc*3.0*(gVb0*gVb0 + gAb0*gAb0);
            double R0c0 = Nc*(gVc0*gVc0 + gAc0*gAc0)/sum;
            double coeff = - 2.0*Nc*R0c0/sum;
            double coeffV = coeff*gVb0;
            double coeffA = coeff*gAb0;
            //std::cout << "cV: " << coeffV << std::endl;
            //std::cout << "cA: " << coeffA << std::endl;
            //std::cout << "cL: " << coeffV+coeffA << std::endl;
            //std::cout << "cR: " << coeffV-coeffA << std::endl;

            R0_c += coeffV*SM.deltaGVb() + coeffA*SM.deltaGAb();
        }   
    
        /* TEST */
        //R0_c -= myEW.Gamma_q(SM.CHARM)/myEW.Gamma_had();
    }

    return R0_c;
}
        

