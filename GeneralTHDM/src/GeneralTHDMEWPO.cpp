/* 
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */


/*FALTA PER ARREGLAR:
 2. Cb continua com un imput ací
 3. Fer-ho general*/




/*Rb as defined in arXiv:hep-ph/9909335,  arXiv:hep-ph/9801355 and arXiv:1006.0470 */


#include "GeneralTHDMcache.h"
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
    
    if (!myGTHDM.getATHDMflag())
    {
        throw std::runtime_error("Rb only aviable in the A2HDM.");
    }
    
    double GF = myGTHDM.getGF();
    double as = myGTHDM.getAlsMz();
    double a = myGTHDM.getAle();

    double MW=myGTHDM.Mw();

    double MZ=myGTHDM.getMz();
    
    /*square mass of the charged Higgs*/
    
    double mHp2=myGTHDM.getmHp2();
    
    /*up and down couplings */
    gslpp::complex u = myGTHDM.getNu_11();
    gslpp::complex d = myGTHDM.getNd_11();
    
    /*The following inputs are needed but are not defined before*/
    
    /** From: arXiv: 1002.1071v2 and "QCD corrections to e+e- cross section and Z decay rates: Concepts and results" **/
    
    double Cb = 1.0086;
 
    /**Masses of the top and bottom quarks at the scale of the W mass (mu) **/
    
    double mu = MW;
    double mtMZ = myGTHDM.Mrun(mu, myGTHDM.getQuarks(QCD::TOP).getMass_scale(),
                              myGTHDM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    double mbMZ = myGTHDM.Mrun(mu, myGTHDM.getQuarks(QCD::BOTTOM).getMass_scale(),
                              myGTHDM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);


    double xH = (mtMZ*mtMZ)/(mHp2);
    double pi=M_PI;
    
    
    double Qb = myGTHDM.getQuarks(QCD::BOTTOM).getCharge();
    double Qt = myGTHDM.getQuarks(QCD::TOP).getCharge();
    

       
    double gAb = SM.getQuarks(QCD::BOTTOM).getIsospin();
    double gVb = SM.getQuarks(QCD::BOTTOM).getIsospin()-2.0*SM.getQuarks(QCD::BOTTOM).getCharge()*SM.sW2();
    double gbRSM = (gVb - gAb)/2.;
    double gbLSM = (gVb + gAb)/2.;

    double gVU = SM.getQuarks(QCD::UP).getIsospin();
    double gAU = SM.getQuarks(QCD::UP).getIsospin()-2.0*SM.getQuarks(QCD::UP).getCharge()*SM.sW2();
    double gbLSMU = (1.0/2.0)*(gVU+gAU);
    double gbRSMU = (1.0/2.0)*(gVU-gAU);
    
    double gVD = SM.getQuarks(QCD::DOWN).getIsospin();
    double gAD = SM.getQuarks(QCD::DOWN).getIsospin()-2.0*SM.getQuarks(QCD::DOWN).getCharge()*SM.sW2();
    double gbLSMD = (1.0/2.0)*(gVD+gAD);
    double gbRSMD = (1.0/2.0)*(gVD-gAD);
    
    double gVC = SM.getQuarks(QCD::CHARM).getIsospin();
    double gAC = SM.getQuarks(QCD::CHARM).getIsospin()-2.0*SM.getQuarks(QCD::CHARM).getCharge()*SM.sW2();
    double gbLSMC = (1.0/2.0)*(gVC+gAC);
    double gbRSMC = (1.0/2.0)*(gVC-gAC);
    
    double gVS = SM.getQuarks(QCD::STRANGE).getIsospin();
    double gAS = SM.getQuarks(QCD::STRANGE).getIsospin()-2.0*SM.getQuarks(QCD::STRANGE).getCharge()*SM.sW2();
    double gbLSMS = (1.0/2.0)*(gVS+gAS);
    double gbRSMS = (1.0/2.0)*(gVS-gAS);
    

    
    double su = ((gbLSMU - gbRSMU)*(gbLSMU - gbRSMU)+ (gbLSMU+ gbRSMU)*(gbLSMU+ gbRSMU))*(1.0 + 3.0*a*Qt*Qt/(4*pi));
    double sd = ((gbLSMD- gbRSMD)*(gbLSMD- gbRSMD)+ (gbLSMD+ gbRSMD)*(gbLSMD+ gbRSMD))*(1.0 + 3.0*a*Qb*Qb/(4*pi));
    double sc = ((gbLSMC- gbRSMC)*(gbLSMC- gbRSMC)+ (gbLSMC+ gbRSMC)* (gbLSMC+ gbRSMC))*(1.0 + 3.0*a*Qt*Qt/(4*pi));
    double ss = ((gbLSMS- gbRSMS)*(gbLSMS- gbRSMS)+ (gbLSMS+ gbRSMS)* (gbLSMS+ gbRSMS))*(1.0 + 3.0*a*Qb*Qb/(4*pi));
    
    double Sb = su+sd+sc+ss;
    
  
        
        /** Loop functions **/
        
        double f1 =(xH*xH - xH - xH*log(xH))/((1.0 - xH)*(1.0 - xH));
        
        double f2 =  - (6.0*xH*(xH-2.0)*gsl_sf_dilog(1.0-1./xH))/((xH-1.0)*(xH-1.0))+ xH*(-27.0+11.0*xH)/((xH-1.0)*(xH-1.0))+ xH*(25.0-9.0*xH)*log(xH)/((xH-1.0)*(xH-1.0)*(xH-1.0))+ (6.0*xH*(3.0-xH)/((xH-1.0)*(xH-1.0)) - (12.0*xH*log(xH))/((xH-1.0)*(xH-1.0)*(xH-1.0)))*log((mtMZ*mtMZ)/(MZ*MZ)) - 3.0*f1;
        
        /*Couplings*/
        
        double gbL= gbLSM + (sqrt(2)*GF*mtMZ*mtMZ*(u).abs2())/(16*pi*pi)*(f1 + (as/(3*pi))*f2);
        
        double gbR = gbRSM - (sqrt(2)*GF*mbMZ*mbMZ*(d).abs2())/(16*pi*pi)*(f1+ (as/(3*pi))*f2);
        double sb = ((gbL- gbR)*(gbL- gbR)+ (gbL+ gbR)*(gbL+ gbR))*(1.0 + 3.0*a*Qb*Qb/(4*pi));
        
        /*Rb*/
        
        double Rb0GGTHDM = 1.0/(1.0 + (Sb*Cb)/sb);
        return Rb0GGTHDM;

       
       /* double DeltaRb0=1.0;
        double Rb0SM=myGTHDM.R0_f(SM.getQuarks(SM.BOTTOM));
        return Rb0SM+DeltaRb0;*/
    
}
