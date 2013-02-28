/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AFBbottom.h"


double AFBbottom::getThValue() 
{   
    double AFB_b;
    if (myEWTYPE==EW::EWCHMN) 
        AFB_b = myEW.getMyEW_CHMN().AFB_q(SM.BOTTOM);
    else if (myEWTYPE==EW::EWABC || myEWTYPE==EW::EWABC2) 
        AFB_b = myEW.getMyEW_ABC().AFB_b(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else {
        AFB_b = 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
        
        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //AFB_b = 0.1032;
                
                AFB_b += - 0.0188*myEW.S() + 0.0131*myEW.T();
            } else {
                double alpha = myEW.alpha0();          
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                double s6 = s4*s2;        
                double s8 = s6*s2;
            
                AFB_b -= 18.0*alpha*s2*(15.0-76.0*s2+152.0*s4-160.0*s6+64.0*s8)
                         /pow(1.0-4.0*s2+8.0*s4, 2.0)
                         /pow(9.0-12.0*s2+8.0*s4, 2.0)/(c2-s2)
                         *( myEW.S() - 4.0*c2*s2*myEW.T() ); 
            }
        }
        
        if (SM.IsFlagNPZbbbarLinearize() && (SM.deltaGVb()!=0.0 || SM.deltaGAb()!=0.0) ) {
            double gVb0 = SM.getQuarks(SM.BOTTOM).getIsospin() 
                          - 2.0*SM.getQuarks(SM.BOTTOM).getCharge()*myEW.s02();
            double gAb0 = SM.getQuarks(SM.BOTTOM).getIsospin();        
            double gVe0 = SM.getLeptons(SM.ELECTRON).getIsospin() 
                          - 2.0*SM.getLeptons(SM.ELECTRON).getCharge()*myEW.s02();
            double gAe0 = SM.getLeptons(SM.ELECTRON).getIsospin();        
            double coeff = - 3.0*gVe0*gAe0*(gVb0*gVb0 - gAb0*gAb0)
                           /(gVe0*gVe0 + gAe0*gAe0)
                           /(gVb0*gVb0 + gAb0*gAb0)/(gVb0*gVb0 + gAb0*gAb0);
            double coeffV = coeff*gAb0;
            double coeffA = - coeff*gVb0;
            //std::cout << "cV: " << coeffV << std::endl;
            //std::cout << "cA: " << coeffA << std::endl;
            //std::cout << "cL: " << coeffV+coeffA << std::endl;
            //std::cout << "cR: " << coeffV-coeffA << std::endl;

            AFB_b += coeffV*SM.deltaGVb() + coeffA*SM.deltaGAb();
        }
        
        /* TEST */
        //AFB_b -= 3.0/4.0*myEW.A_l(SM.ELECTRON)*myEW.A_q(SM.BOTTOM);
    }

    //std::cout << "EWTYPE = " << myEW.getEWTYPE() << std::endl; // TEST

    return AFB_b;
}
