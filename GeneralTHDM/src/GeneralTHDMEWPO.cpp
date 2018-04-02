/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


/*Rb as defined in arXiv:hep-ph/9909335,  arXiv:hep-ph/9801355 and arXiv:1006.0470 */



#include "GeneralTHDMEWPO.h"
#include "StandardModel.h"
#include <gsl/gsl_sf_dilog.h>
#include <math.h>
#include "GeneralTHDM.h"





Rb0GTHDM::Rb0GTHDM(const StandardModel& SM_i)
: ThObservable(SM_i),myGTHDM(static_cast<const GeneralTHDM&> (SM_i)) 
{}




double Rb0GTHDM::computeThValue()
{
    
    double GF = myGTHDM.getGF();
    double as = myGTHDM.getAlsMz();
    double a = myGTHDM.getAle();
//    double MW=myGTHDM.Mw();
    double MZ=myGTHDM.getMz();
    
    /*square mass of the charged Higgs*/
    
    double mHp2=myGTHDM.getmHp2();
    
    /*up and down couplings */
    gslpp::complex u = myGTHDM.getNu_11();
    gslpp::complex d = myGTHDM.getNd_22();
    
    /*The following inputs are needed but are not defined before*/
    /** From arXiv: hep-ex/0509008 and arXiv: 1002.1071 **/
    double gbLSM = -0.42112;
//    double upgbLSM = 0.00035;
//    double umgbLSM = 0.00018;
    double gbRSM = 0.07744;
//    double upgbRSM = 0.00006;
//    double umgbRSM = 0.00008;
    
    /** From: arXiv: 1002.1071v2 and "QCD corrections to e+e- cross section and Z decay rates: Concepts and results" **/
    
    double Sb = 1.3214;
    double Cb = 1.0086;
    double Qb = -1.0/3.0;
    
    /**Masses of the top and bottom quarks at the scale of the W mass (mu) **/
    
    double mtMZ = 171.1;
    double mbMZ = 2.89;
//    double mu = MW;

    double xH = (mtMZ*mtMZ)/(mHp2);
    double pi=M_PI;



    
    if (!myGTHDM.getATHDMflag())
    {
        std::cout << "Rb is only available in the ATHDM at the moment.";
        double Rb0GGTHDM = 0.0;
        return Rb0GGTHDM;
    }
    else
    {
        
        /** Loop functions **/
        
        double f1 =(xH*xH - xH - xH*log(xH))/((1.0 - xH)*(1.0 - xH));
        
        double f2 =  - (6.0*xH*(xH-2.0)*gsl_sf_dilog(1.0-1./xH))/((xH-1.0)*(xH-1.0))+ xH*(-27.0+11.0*xH)/((xH-1.0)*(xH-1.0))+ xH*(25.0-9.0*xH)*log(xH)/((xH-1.0)*(xH-1.0)*(xH-1.0))+ (6.0*xH*(3.0-xH)/((xH-1.0)*(xH-1.0)) - (12.0*xH*log(xH))/((xH-1.0)*(xH-1.0)*(xH-1.0)))*log((mtMZ*mtMZ)/(MZ*MZ)) - 3.0*f1;
        
        /*Couplings*/
        
        double gbL =  gbLSM + (sqrt(2)*GF*mtMZ*mtMZ*(u*u).abs())/(16*pi*pi)*(f1 + (as/(3*pi))*f2);
        
        double gbR = gbRSM - (sqrt(2)*GF*mbMZ*mbMZ*(d*d).abs())/(16*pi*pi)*(f1+ (as/(3*pi))*f2);
        double sb = ((gbL- gbR)*(gbL- gbR)+ (gbL+ gbR)*(gbL+ gbR))*(1.0 + 3.0*a*Qb*Qb/(4*pi));
        
        /*Rb*/
        
        double Rb0GGTHDM = 1.0/(1.0 + (Sb*Cb)/sb);
        return Rb0GGTHDM;
    }
       /* double DeltaRb0=1.0;
        double Rb0SM=myGTHDM.R0_f(SM.getQuarks(SM.BOTTOM));
        return Rb0SM+DeltaRb0;*/
    
}
