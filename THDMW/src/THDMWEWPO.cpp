/*
 * Copyright (C) 2018 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */



#include "THDMWEWPO.h"
#include "THDMWcache.h"
#include "THDMW.h"
#include <gsl/gsl_sf_dilog.h>
#include <math.h>


Rb0THDMW::Rb0THDMW(const StandardModel& SM_i)
: ThObservable(SM_i),myTHDMW(static_cast<const THDMW&> (SM_i)) 
{}




double Rb0THDMW::computeThValue()
{
    //Review this at some point -> use HEP-FIT values, not the ones included by hand. For doing so a SM fit should be done
    // For the implementation we follow the reference arxiv:1002.1071v2
    double as = myTHDMW.getAlsMz();
    double a = myTHDMW.getAle();
    //double a=1/127.916;//fine-structure at the mass of the Z boson
    double MW=myTHDMW.Mw();
    double T3b=myTHDMW.getQuarks(QCD::BOTTOM).getIsospin();
    double Qb=myTHDMW.getQuarks(QCD::BOTTOM).getCharge();
    double CF=4./3.;
    double CR1=CF;
    //double CR1=1;
    double CR2=3.;
    //double CR2=0;
    double Cb=1.0086;
    double Rb=0.;
    //double rho_SM=1.00578;//up and charm value
    //double rho_SM=1.00675;//down and strange value
    //double rhob_SM=0.99426;
    double rhob_SM=myTHDMW.getTHDMW_rho_b();
    double Sb=myTHDMW.getTHDMW_S_b();
    //double sW2_SM=myTHDMW.sW2();
    double sW2b_SM= 0.23298;
    double mSpsq=myTHDMW.getMyTHDMWCache() -> mSpsq;
    
    double eta_U=myTHDMW.getTHDMW_etaU();
    double eta_D=myTHDMW.getTHDMW_etaD();
    double MZ = myTHDMW.getMz();
    double mu=MZ;
    double mtMZ = myTHDMW.Mrun(mu, myTHDMW.getQuarks(QCD::TOP).getMass_scale(),
                              myTHDMW.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    double mbMZ = myTHDMW.Mrun(mu, myTHDMW.getQuarks(QCD::BOTTOM).getMass_scale(),
                              myTHDMW.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double x=pow(mtMZ,2)/mSpsq;
    
    //std::cout<<"orders_EW_size="<< SM.orders_EW_size << std::endl;
    //std::cout<<"mu="<< mu << std::endl;
    //std::cout<<"log "<< log(10) << std::endl;
    //std::cout<<"fine-structure "<< a << std::endl;
    //std::cout<<"mSpsq= "<< mSpsq << std::endl;
    
    double f1=x/(x-1.)- x*log(x)/pow((x-1.),2);
    double f2=-6.*x*(x-2.)*gsl_sf_dilog(1.0-1./x)/pow((x-1.),2)+x*(-27.+11.*x)/pow((x-1.),2)+x*(25.-9.*x)*log(x)/pow(x-1.,3)+(6.*x*(3.-x)/pow(x-1.,2)-12.*x*log(x)/pow(x-1.,3))*log(pow(mtMZ,2)/pow(mu,2))-3.*f1;
    double f3=3.*x*gsl_sf_dilog(1.0-1./x)/(x-1.)+3.*x*(1.-2.*x+pow(x,2)+pow(log(x),2))/pow(x-1.,3)-6.*x*log(x)/pow(x-1.,2);
    //std::cout<<"x= "<< x << std::endl;
    //std::cout<<"(1.0-1./x)= "<< (1.0-1./x) << std::endl;
    //std::cout<<"gsl_sf_dilog(1.0-1./x)= "<< gsl_sf_dilog(1.0-1./x) << std::endl;
    
    
    double delrhob=(1/T3b)*(a/(4*M_PI*sW2b_SM))*CR1*(pow(sqrt(pow(eta_U,2))*mtMZ/(sqrt(2)*MW),2)+pow(sqrt(pow(eta_D,2))*mbMZ/(sqrt(2)*MW),2))*(f1+(as/(4*M_PI))*(CF*f2+CR2*f3));
    double delsW2b=-0.5*delrhob*sW2b_SM+(1/(2*Qb))*(a/(4*M_PI*sW2b_SM))*CR1*pow(sqrt(pow(eta_D,2))*mbMZ/(sqrt(2)*MW),2)*(f1+(as/(4*M_PI))*(CF*f2+CR2*f3));
    
    //double rhob=rho_SM+delrhob;
    //double sW2b=sW2_SM+delsW2b;

    
    
    
    
    
    //double gLu=sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::UP).getIsospin()-myTHDMW.getQuarks(QCD::UP).getCharge()*sW2_SM);
    //double gRu = -sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::UP).getCharge()*sW2_SM);

    //double gLd=sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::DOWN).getIsospin()-myTHDMW.getQuarks(QCD::DOWN).getCharge()*sW2_SM);
    //double gRd = -sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::DOWN).getCharge()*sW2_SM);
    
    //double gLc=sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::CHARM).getIsospin()-myTHDMW.getQuarks(QCD::CHARM).getCharge()*sW2_SM);
    //double gRc = -sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::CHARM).getCharge()*sW2_SM);

    //double gLs=sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::STRANGE).getIsospin()-myTHDMW.getQuarks(QCD::STRANGE).getCharge()*sW2_SM);
    //double gRs = -sqrt(rho_SM)*(myTHDMW.getQuarks(QCD::STRANGE).getCharge()*sW2_SM);
    
    double gLb=sqrt(rhob_SM)*(T3b-Qb*(sW2b_SM+delsW2b))+sqrt(rhob_SM)*(delrhob/(2.*rhob_SM))*(T3b-Qb*sW2b_SM);
    double gRb = -sqrt(rhob_SM)*Qb*(sW2b_SM+delsW2b+sW2b_SM*delrhob/(2.*rhob_SM));

    
    
    
    
    //Values taken from https://arxiv.org/pdf/hep-ex/0509008.pdf
//    double gLu=0.34674;
//    double gRu = -0.15470;
//    
//    double gLd=-0.42434;
//    double gRd =0.077345;
//    
//    double gLc=0.34674;
//    double gRc = -0.15470;
//
//    double gLs=-0.42434;
//    double gRs =0.077345;
    

    
    
    
    
    
//    double su= (pow(gLu,2)+pow(gRu,2))*(1.+pow(myTHDMW.getQuarks(QCD::UP).getCharge(),2)*3.*a/(4.*M_PI));
//    double sd= (pow(gLd,2)+pow(gRd,2))*(1.+pow(myTHDMW.getQuarks(QCD::DOWN).getCharge(),2)*3.*a/(4.*M_PI));
//    double sc= (pow(gLc,2)+pow(gRc,2))*(1.+pow(myTHDMW.getQuarks(QCD::CHARM).getCharge(),2)*3.*a/(4.*M_PI));
//    double ss= (pow(gLs,2)+pow(gRs,2))*(1.+pow(myTHDMW.getQuarks(QCD::STRANGE).getCharge(),2)*3.*a/(4.*M_PI));
    //double Sb=su+sd+sc+ss;
    //std::cout<<"glu= "<< gLu << std::endl;
    //std::cout<<"gru= "<< gRu << std::endl;
    //std::cout<<"gld= "<< gLd << std::endl;
    //std::cout<<"grd= "<< gRd << std::endl;
    //std::cout<<"glc= "<< gLc << std::endl;
    //std::cout<<"grc= "<< gRc << std::endl;
    //std::cout<<"gls= "<< gLs << std::endl;
    //std::cout<<"grs= "<< gRs << std::endl;
    //std::cout<<"glb= "<< gLb << std::endl;
    //std::cout<<"grb= "<< gRb << std::endl;
    //std::cout<<"sW2= "<< sW2_SM << std::endl;
    //std::cout<<"Sb= "<< Sb << std::endl;
    double sb= (pow(gLb,2)+pow(gRb,2))*(1.+pow(Qb,2)*3.*a/(4.*M_PI));
    
    Rb=1/(1+Sb*Cb/sb);
    
    
    
    return Rb;
}