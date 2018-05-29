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
    
    gslpp::complex gVB = myGTHDM.gV_f(myGTHDM.getQuarks(QCD::BOTTOM));
    gslpp::complex gAB = myGTHDM.gA_f(myGTHDM.getQuarks(QCD::BOTTOM));
    

    gslpp::complex gbLSM = (1.0/2.0)*(gVB+gAB);
    gslpp::complex gbRSM = (1.0/2.0)*(gVB-gAB);
    
    
    gslpp::complex gVU = myGTHDM.gV_f(myGTHDM.getQuarks(QCD::UP));
    gslpp::complex gAU = myGTHDM.gA_f(myGTHDM.getQuarks(QCD::UP));
    gslpp::complex gbLSMU = (1.0/2.0)*(gVU+gAU);
    gslpp::complex gbRSMU = (1.0/2.0)*(gVU-gAU);
    
    gslpp::complex gVD = myGTHDM.gV_f(myGTHDM.getQuarks(QCD::DOWN));
    gslpp::complex gAD = myGTHDM.gA_f(myGTHDM.getQuarks(QCD::DOWN));
    gslpp::complex gbLSMD = (1.0/2.0)*(gVD+gAD);
    gslpp::complex gbRSMD = (1.0/2.0)*(gVD-gAD);
    
    gslpp::complex gVC = myGTHDM.gV_f(myGTHDM.getQuarks(QCD::CHARM));
    gslpp::complex gAC = myGTHDM.gA_f(myGTHDM.getQuarks(QCD::CHARM));
    gslpp::complex gbLSMC = (1.0/2.0)*(gVC+gAC);
    gslpp::complex gbRSMC = (1.0/2.0)*(gVC-gAC);
    
    gslpp::complex gVS = myGTHDM.gV_f(myGTHDM.getQuarks(QCD::STRANGE));
    gslpp::complex gAS = myGTHDM.gA_f(myGTHDM.getQuarks(QCD::STRANGE));
    gslpp::complex gbLSMS = (1.0/2.0)*(gVS+gAS);
    gslpp::complex gbRSMS = (1.0/2.0)*(gVS-gAS);
    
    double su = ((gbLSMU - gbRSMU).abs2()+ (gbLSMU+ gbRSMU).abs2())*(1.0 + 3.0*a*Qt*Qt/(4*pi));
    double sd = ((gbLSMD- gbRSMD).abs2()+ (gbLSMD+ gbRSMD).abs2())*(1.0 + 3.0*a*Qb*Qb/(4*pi));
    double sc = ((gbLSMC- gbRSMC).abs2()+ (gbLSMC+ gbRSMC).abs2())*(1.0 + 3.0*a*Qt*Qt/(4*pi));
    double ss = ((gbLSMS- gbRSMS).abs2()+ (gbLSMS+ gbRSMS).abs2())*(1.0 + 3.0*a*Qb*Qb/(4*pi));
    
    double Sb = su+sd+sc+ss;

    
    if (!myGTHDM.getATHDMflag())
    {
        std::cout << "Rb is only available in the ATHDM at the moment.";
        double Rb0GGTHDM = 0.0;
        return Rb0GGTHDM;
    }
    else{
        
        /** Loop functions **/
        
        double f1 =(xH*xH - xH - xH*log(xH))/((1.0 - xH)*(1.0 - xH));
        
        double f2 =  - (6.0*xH*(xH-2.0)*gsl_sf_dilog(1.0-1./xH))/((xH-1.0)*(xH-1.0))+ xH*(-27.0+11.0*xH)/((xH-1.0)*(xH-1.0))+ xH*(25.0-9.0*xH)*log(xH)/((xH-1.0)*(xH-1.0)*(xH-1.0))+ (6.0*xH*(3.0-xH)/((xH-1.0)*(xH-1.0)) - (12.0*xH*log(xH))/((xH-1.0)*(xH-1.0)*(xH-1.0)))*log((mtMZ*mtMZ)/(MZ*MZ)) - 3.0*f1;
        
        /*Couplings*/
        
        double gbL =  gbLSM.real() + (sqrt(2)*GF*mtMZ*mtMZ*(u*u).abs())/(16*pi*pi)*(f1 + (as/(3*pi))*f2);
        
        double gbR = gbRSM.real() - (sqrt(2)*GF*mbMZ*mbMZ*(d*d).abs())/(16*pi*pi)*(f1+ (as/(3*pi))*f2);
        double sb = ((gbL- gbR)*(gbL- gbR)+ (gbL+ gbR)*(gbL+ gbR))*(1.0 + 3.0*a*Qb*Qb/(4*pi));
        
        /*Rb*/
        
        double Rb0GGTHDM = 1.0/(1.0 + (Sb*Cb)/sb);
        return Rb0GGTHDM;

        }
       /* double DeltaRb0=1.0;
        double Rb0SM=myGTHDM.R0_f(SM.getQuarks(SM.BOTTOM));
        return Rb0SM+DeltaRb0;*/
    
}
