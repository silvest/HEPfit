/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "sigmaHadron.h"


double sigmaHadron::getThValue() 
{ 
    double sigma_had;
    if (myEWTYPE==EW::EWCHMN)  
        sigma_had = myEW.getMyEW_CHMN().sigma0_had();
    else if (myEWTYPE==EW::EWABC) 
        sigma_had = myEW.getMyEW_ABC().sigma0_had(SM.epsilon1(),SM.epsilon3(),SM.epsilonb());
    else if (myEWTYPE==EW::EWABC2) {
        double delta_als = (SM.Als(SM.getMz(),FULLNNLO) - 0.119)/M_PI;
        double delta_alpha = (SM.alphaMz() - 1.0/128.90)/SM.getAle();
        double sigma_h0 = 41.420*(1.0 - 0.41*delta_als + 0.03*delta_alpha)/GeVminus2_to_nb;
        sigma_had = sigma_h0*(1.0 - 0.03*SM.epsilon1() + 0.04*SM.epsilon3() - 0.20*SM.epsilonb());
    } else {   
        sigma_had = myEW.sigma0_had();
        
        if ( myEW.checkModelForSTU() ) {
            if(myEWTYPE==EW::EWBURGESS) {
                // TEST: the fit result by Gfitter in arXiv:1209.2716, 
                //       corresponding to MH=125.7 and Mt=173.52 
                //sigma_had = 41.479/GeVminus2_to_nb;                
                
                double delta_l = - 0.000192*myEW.S() + 0.000790*myEW.T();
                double delta_had = - 0.00901*myEW.S() + 0.0200*myEW.T();
                double delta_Z = - 0.00961*myEW.S() + 0.0263*myEW.T();
                sigma_had *= 1.0 + delta_l/myEW.Gamma_l(SM.ELECTRON)
                             + delta_had/myEW.Gamma_had()
                             - 2.0*delta_Z/myEW.Gamma_Z();
            } else {
                double alpha = myEW.alpha0();  
                double Mz = SM.getMz();
                double c2 = myEW.c02();
                double s2 = myEW.s02();
                double s4 = s2*s2;
                double s6 = s4*s2;        
                double s8 = s6*s2;
                
                sigma_had -= 72.0*M_PI*alpha
                             *(729.0-4788.0*s2+8352.0*s4-6176.0*s6+640.0*s8)
                             /Mz/Mz/pow(63.0-120.0*s2+160.0*s4, 3.0)/(c2-s2)
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
            double gVnu0 = SM.getLeptons(SM.NEUTRINO_1).getIsospin() 
                          - 2.0*SM.getLeptons(SM.NEUTRINO_1).getCharge()*myEW.s02();
            double gAnu0 = SM.getLeptons(SM.NEUTRINO_1).getIsospin();        
            double gVe0 = SM.getLeptons(SM.ELECTRON).getIsospin() 
                          - 2.0*SM.getLeptons(SM.ELECTRON).getCharge()*myEW.s02();
            double gAe0 = SM.getLeptons(SM.ELECTRON).getIsospin();        
            double Nc = 3.0;
            double sum_q = Nc*2.0*(gVu0*gVu0 + gAu0*gAu0)
                           + Nc*3.0*(gVb0*gVb0 + gAb0*gAb0);
            double sum_f = 3.0*(gVnu0*gVnu0 + gAnu0*gAnu0) 
                           + 3.0*(gVe0*gVe0 + gAe0*gAe0)
                           + sum_q;
            double sigmah0 = 12.0*M_PI/SM.getMz()/SM.getMz()
                             *(gVe0*gVe0 + gAe0*gAe0)*sum_q/sum_f/sum_f;
            double coeff = 2.0*Nc*sigmah0*(1.0/sum_q - 2.0/sum_f);
            double coeffV = coeff*gVb0;
            double coeffA = coeff*gAb0;
            //std::cout << "cV: " << coeffV << std::endl;
            //std::cout << "cA: " << coeffA << std::endl;
            //std::cout << "cL: " << coeffV+coeffA << std::endl;
            //std::cout << "cR: " << coeffV-coeffA << std::endl;

            sigma_had += coeffV*SM.deltaGVb() + coeffA*SM.deltaGAb();
        }      
        
        /* TEST */
        //sigma_had -= myEW.sigma0_had();
    }
    
    return ( sigma_had*GeVminus2_to_nb );
}
        


