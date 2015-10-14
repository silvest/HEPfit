/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SUSYMassInsertionMatching.h"
#include "SUSYMassInsertion.h" 

SUSYMassInsertionMatching::SUSYMassInsertionMatching(const SUSYMassInsertion & SUSYMassInsertion_i):
        
    StandardModelMatching(SUSYMassInsertion_i), 
    SusyMI(SUSYMassInsertion_i), 
        
    drNDRLRI(5, 5, 0.), 
    mcd2(5, NDR, NLO), 
    mcd1(10, NDR, NLO), 
    mcbd(5, NDR, NLO), mcbs(5, NDR, NLO), mck2(5, NDR, NLO) {
    
    Nf = 6;
    
    double Nc = SusyMI.getNc();
    
    drNDRLRI(0,0) = -(((-1. + Nc)*(-7. + log(4096.))) / Nc);                   
    drNDRLRI(1,1) = (-2.*(-1. + 6.*Nc*Nc - 8.*log(2.) + Nc*(-13. + log(1024.))))/(3.*Nc);
    drNDRLRI(1,2) = (-2.*(13. - 10.*log(2.) + Nc*(-5. + log(256.))))/(3.*Nc);   
    drNDRLRI(2,1) = (-8. + 6.*Nc*Nc + 20.*log(2.) - 8.*Nc*(1. + log(4.)))/(3.*Nc);
    drNDRLRI(2,2) = (2.*(4. + Nc - 10.*Nc*log(2.) + log(256.)))/(3.*Nc);
    drNDRLRI(3,3) = (2. - 4.*Nc*Nc + log(4.))/Nc;
    drNDRLRI(3,4) = 2. - log(4.);
    drNDRLRI(4,3) = -2.*(1. + log(2.));
    drNDRLRI(4,4) = (2. + log(4.))/Nc;

}

// delta F = 1 loop functions

double SUSYMassInsertionMatching::B1(double x) const {
    return ((1. + 4.*x - 5.*x*x + 4.*x*log(x) + 2.*x*x*log(x))/(8.*pow(x-1.,4.)));
}

double SUSYMassInsertionMatching::B2(double x) const {
    return (x*(5. - 4.*x - x*x + 2.*log(x) + 4.*x*log(x))/(2.*pow(x-1.,4.)));
}

double SUSYMassInsertionMatching::P1(double x) const {
    return ((1. - 6.*x + 18.*x*x - 10.*x*x*x - 3.*x*x*x*x + 
            12.*x*x*x*log(x))/(18.*pow(x-1.,5.)));
}

double SUSYMassInsertionMatching::P2(double x) const {
    return ((7. - 18.*x + 9.*x*x + 2.*x*x*x +
             3.*log(x) - 9.*x*x*log(x))/(9.*pow(x-1.,5.)));
}

double SUSYMassInsertionMatching::M1(double x) const{
  return((1. + 4. * x - 5. * x*x + 4.*x*log(x) + 2.*x*x*log(x))/(2.*pow(x-1.,4.)));
}

double SUSYMassInsertionMatching::M2(double x) const{
  return((-5. + 4.*x + x*x - 2.*log(x) - 4.*x*log(x))/(2.*pow(x-1.,4.)));
}

double SUSYMassInsertionMatching::M3(double x) const{
  return((-1. + 9.*x + 9.*x*x - 17.*x*x*x + 18.*x*x*log(x) + 6.*x*x*x*log(x))/(12.*pow(x-1.,5.)));
}

double SUSYMassInsertionMatching::M4(double x) const {
  return((-1.-9.*x + 9.*x*x + x*x*x - 6.*x*log(x) - 6.*x*x* log(x))/(6.*pow(x-1.,5.)));
}

// delta F = 2 loop functions

double SUSYMassInsertionMatching::f6(double x){
    return(x*(17. - 9.*x - 9.*x*x + x*x*x + 18.*x*log(x) +
           6.*log(x))/(6.*pow(x-1.,5.)));
}

double SUSYMassInsertionMatching::f6t(double x){
    return((1. + 9.*x - 9.*x*x - x*x*x + 6.*x*x*log(x) +
	   6.*x*log(x))/(3.*pow(x-1.,5.)));
}

//the LO terms have to be multiplyed by als^2(mumatch)/Msq^2
/**
 * Wilson Coefficients in RI-MOM scheme
 *  
 */

double SUSYMassInsertionMatching::C0LO(double x){
    //return(16./pow(1.-x,5.)*(11./1728. + 133.*x/1728. - 13.*x*x/192. - 29.*x*x*x/1728. +
    //       x*x*x*x/864 + 13.*x*log(x)/288. + 17.*x*x*log(x)/288.));
 double pow5 = pow(1.-x,5.); 
 return( (1./pow5) *
         (11./108.+133./108.*x-13./12.*x*x-29./108.*x*x*x+1./54.*x*x*x*x+(13./18.*x+17./18.*x*x)*log(x)
	  ) );
}


double SUSYMassInsertionMatching::C1LO(double x){
    return(16.*(289.*x/1728. - 17.*x*x/192. - 17.*x*x*x/192. + 17.*x*x*x*x/1728. + 
            17.*x*log(x)/288. + 17.*x*x*log(x)/96. )/pow(1.-x,5.));
}

double SUSYMassInsertionMatching::C2LO(double x){
    return(16.*(-17.*x/576. + x*x/64. + x*x*x/64. - x*x*x*x/576. - x*log(x)/96. - 
           x*x*log(x)/32.)/(pow(1.-x,5)));
}

double SUSYMassInsertionMatching::C3LOA(double x){
    return(16.*(-11./864. - 11.*x/96. + 11.*x*x/96. + 11.*x*x*x/864. - 11.*x*log(x)/144. - 
           11.*x*x*log(x)/144.)/pow(1.-x,5.) );
}

double SUSYMassInsertionMatching::C3LOB(double x){
    return(16.*(-1./144. + 101.*x/288. - 5.*x*x/32. - 61.*x*x*x/288. + 7.*x*x*x*x/288. 
           + 5.*x*log(x)/48. + 19.*x*x*log(x)/48.)/(pow(1.-x,5.)));
}

double SUSYMassInsertionMatching::C4LOA(double x){    
    return (16.*(-5./288. - 5.*x/32. + 5.*x*x/32. + 5.*x*x*x/288.
            - 5.*x*log(x)/48. - 5.*x*x*log(x)/48.)/(pow(1.-x,5.))           
            );
}

double SUSYMassInsertionMatching::C4LOB(double x){
    return(16.*(5./432. + 107.*x/864. - 11.*x*x/96. - 19.*x*x*x/864. + x*x*x*x/864.
           + 11.*x*log(x)/144. + 13.*x*x*log(x)/144.)/(pow(1.-x,5.)));
}  

//to multiplyed by als^3(mumatch)/M_PI/Msq^2
/**
 * Wilson Coefficients in RI-MOM scheme
 *  
 */

double SUSYMassInsertionMatching::C0NLO(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    double pow4 = pow(1.-x,4.);
    
    return  ( 64.*
            ((-11./(20736.*pow7) + 
            4459./(165888.*pow6) + (209.*Nf)/(27648.*
                  pow6) - 1./(1728.*pow5) + 
            151./(41472.*pow4) - (89.*x)/(20736.*
                  pow7) + (13235.*x)/(62208.*pow6) - (83.*
                  Nf*x)/(6912.*pow6) + (3853.*x)/(15552.*
                  pow5) - (173.*x)/(27648.*pow4) + (77.*
                  x*x)/(2592.*pow7) + (33361.*
                  x*x)/(248832.*pow6) - (33.*Nf*
                  x*x)/(512.*pow6) - (9355.*
                  x*x)/(41472.*pow5) + (217.*
                  x*x)/(82944.*pow4) - (419.*
                  x*x*x)/(10368.*pow7) - (7.*
                  x*x*x)/(18.*pow6) + (481.*Nf*
                  x*x*x)/(6912.*pow6) - (503.*
                  x*x*x)/(20736.*pow5) + (233.*
                  x*x*x*x)/(20736.*pow7) + (10085.*
                  x*x*x*x)/(497664.*pow6) - (19.*Nf*
                  x*x*x*x)/(27648.*pow6) + (331.*
                  x*x*x*x)/(124416.*pow5) + (95.*
                  x*x*x*x*x)/(20736.*pow7) - (91.*
                  x*x*x*x*x)/(15552.*pow6) - 
            x*x*x*x*x*x/(3456.*pow7) + 
            x*x*x*x*x*x/(1296.*pow6) - (203.*Li2)/(13824.*
                  pow6) - (11.*Nf*Li2)/(2304.*
                  pow6) + (47.*Li2)/(4608.*
                  pow5) + (53.*x*Li2)/(864.*
                  pow6) - (Nf*x*Li2)/(192.*
                  pow6) + (367.*x*Li2)/(6912.*
                  pow5) - (85.*x*Li2)/(13824.*
                  pow4) - (127.*x*x*Li2)/(13824.*
                  pow6) + (19.*Nf*x*x*Li2)/(768.*
                  pow6) - (35.*x*x*Li2)/(13824.*
                  pow5) - (259.*x*x*x*Li2)/(6912.*
                  pow6) - (17.*Nf*x*x*x*Li2)/(1152.*
                  pow6) + (11.*a)/(3456.*
                  pow5) + (133.*x*a)/(3456.*
                  pow5) - (13.*x*x*a)/(384.*
                  pow5) - (29.*x*x*x*a)/(3456.*
                  pow5) + (x*x*x*x*a)/(1728.*
                  pow5) - (13.*x*log(x))/(3456.*
                  pow7) + (7.*x*log(x))/(256.*
                  pow6) - (Nf*x*log(x))/(1152.*
                  pow6) + (6431.*x*log(x))/(82944.*
                  pow5) - (425.*x*log(x))/(82944.*
                  pow4) + (29.*x*x*log(x))/(2592.*
                  pow7) + (10535.*x*x*log(x))/(41472.*
                  pow6) - (Nf*x*x*log(x))/(48.*
                  pow6) + (1589.*x*x*log(x))/(6912.*
                  pow5) + (55.*x*x*x*log(x))/(2592.*
                  pow7) + (3751.*x*x*x*log(x))/(10368.*
                  pow6) - (25.*Nf*x*x*x*log(x))/(576.*
                  pow6) + (1681.*x*x*x*log(x))/(82944.*
                  pow5) - (5.*x*x*x*x*log(x))/(192.*
                  pow7) - (91.*x*x*x*x*log(x))/(4608.*
                  pow6) - (29.*x*x*x*x*log(x))/(20736.*
                  pow5) - (29.*x*x*x*x*x*log(x))/(10368.*
                  pow7) + (139.*x*x*x*x*x*log(x))/(20736.*
                  pow6) + (x*x*x*x*x*x*log(x))/(5184.*
                  pow7) - (x*x*x*x*x*x*log(x))/(1296.*
                  pow6) + (13.*x*a*log(x))/(576.*
                  pow5) + (17.*x*x*a*log(x))/(576.*
                  pow5) - (73.*x*pow(log(x), 2.))/(3072.*
                  pow6) - (29.*x*pow(log(x), 2.))/(1152.*
                  pow5) - (43.*x*pow(log(x), 2.))/(9216.*
                  pow4) - (1075.*x*x*
                  pow(log(x), 2.))/(4608.*pow6) - (205.*
                  x*x*pow(log(x), 2.))/(2304.*pow5) + (13.*
                  x*x*x*pow(log(x), 2.))/(1728.*
                  pow7) - (3053.*x*x*x*
                  pow(log(x), 2.))/(27648.*pow6) + (17.*
                  x*x*x*x*pow(log(x), 2.))/(1728.*pow7) + 
            log(Ms2/mumatch2)*((11.*Nf)/(13824.*pow6) - 
                  55./(6912.*pow5) - (1153.*x)/(6912.*
                        pow6) + (97.*Nf*x)/(3456.*
                        pow6) - (665.*x)/(6912.*
                        pow5) - (1819.*x*x)/(20736.*
                        pow6) + (Nf*x*x)/(192.*
                        pow6) + (65.*x*x)/(768.*
                        pow5) + (2789.*x*x*x)/(10368.*
                        pow6) - (113.*Nf*x*x*x)/(3456.*
                        pow6) + (145.*x*x*x)/(6912.*
                        pow5) - (143.*x*x*x*x)/(6912.*
                        pow6) - (19.*Nf*x*x*x*x)/(13824.*
                        pow6) - (5.*x*x*x*x)/(3456.*
                        pow5) + (145.*x*x*x*x*x)/(20736.*
                        pow6) - 
                  x*x*x*x*x*x/(1296.*pow6) - (13.*x*log(x))/(256.*
                        pow6) + (13.*Nf*x*log(x))/(1152.*
                        pow6) - (65.*x*log(x))/(1152.*
                        pow5) - (953.*x*x*log(x))/(3456.*
                        pow6) + (5.*Nf*x*x*log(x))/(128.*
                        pow6) - (85.*x*x*log(x))/(1152.*
                        pow5) - (593.*x*x*x*log(x))/(6912.*
                        pow6) + (17.*Nf*x*x*x*log(x))/(1152.*
                        pow6))))
            
            );
}


double SUSYMassInsertionMatching::C1NLO(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    
    return(64.*(((-289.*x)/(20736.*pow7) + (134447.*x)/(248832.*
                  pow6) + (323.*Nf*x)/(6912.*
                  pow6) + (43957.*x)/(82944.*
                  pow5) + (1309.*x*x)/(20736.*
                  pow7) + (35261.*x*x)/(124416.*
                  pow6) - (187.*Nf*x*x)/(768.*
                  pow6) - (30427.*x*x)/(82944.*
                  pow5) - (221.*x*x*x)/(3456.*
                  pow7) - (8107.*x*x*x)/(10368.*
                  pow6) + (51.*Nf*x*x*x)/(256.*
                  pow6) - (15493.*x*x*x)/(82944.*
                  pow5) - (85.*x*x*x*x)/(10368.*
                  pow7) - (1813.*x*x*x*x)/(124416.*
                  pow6) - (17.*Nf*x*x*x*x)/(6912.*
                  pow6) + (1963.*x*x*x*x)/(82944.*
                  pow5) + (527.*x*x*x*x*x)/(20736.*
                  pow7) - (8407.*x*x*x*x*x)/(248832.*
                  pow6) - (17.*x*x*x*x*x*x)/(6912.*
                  pow7) + (17.*x*x*x*x*x*x)/(2592.*
                  pow6) + (89.*x*Li2)/(3456.*
                  pow6) - (17.*Nf*x*Li2)/(384.*
                  pow6) + (2.*x*Li2)/(27.*
                  pow5) + (61.*x*x*Li2)/(576.*
                  pow6) + (17.*Nf*x*x*Li2)/(192.*
                  pow6) - (127.*x*x*Li2)/(864.*
                  pow5) - (455.*x*x*x*Li2)/(3456.*
                  pow6) - (17.*Nf*x*x*x*Li2)/(384.*
                  pow6) + (1411.*x*a)/(31104.*
                  pow5) - (83.*x*x*a)/(3456.*
                  pow5) - (83.*x*x*x*a)/(3456.*
                  pow5) + (83.*x*x*x*x*a)/(31104.*
                  pow5) - (17.*x*log(x))/(3456.*
                  pow7) + (377.*x*log(x))/(3456.*
                  pow6) - (17.*Nf*x*log(x))/(2304.*
                  pow6) + (6443.*x*log(x))/(41472.*
                  pow5) + (17.*x*x*log(x))/(3456.*
                  pow7) + (35.*x*x*log(x))/(81.*
                  pow6) - (17.*Nf*x*x*log(x))/(1152.*
                  pow6) + (16055.*x*x*log(x))/(41472.*
                  pow5) + (187.*x*x*x*log(x))/(2592.*
                  pow7) + (251.*x*x*x*log(x))/(288.*
                  pow6) - (289.*Nf*x*x*x*log(x))/(2304.*
                  pow6) + (35.*x*x*x*log(x))/(576.*
                  pow5) - (17.*x*x*x*x*log(x))/(288.*
                  pow7) - (31.*x*x*x*x*log(x))/(2304.*
                  pow6) - (35.*x*x*x*x*log(x))/(5184.*
                  pow5) - (17.*x*x*x*x*x*log(x))/(1152.*
                  pow7) + (701.*x*x*x*x*x*log(x))/(20736.*
                  pow6) + (17.*x*x*x*x*x*x*log(x))/(10368.*
                  pow7) - (17.*x*x*x*x*x*x*log(x))/(2592.*
                  pow6) + (83.*x*a*log(x))/(5184.*
                  pow5) + (83.*x*x*a*log(x))/(1728.*
                  pow5) - (7.*x*pow(log(x), 2.))/(1536.*
                  pow6) - (17.*x*pow(log(x), 2.))/(1728.*
                  pow5) - (2599.*x*x*
                  pow(log(x), 2.))/(6912.*pow6) - (167.*
                  x*x*pow(log(x), 2.))/(864.*pow5) + (17.*
                  x*x*x*pow(log(x), 2.))/(1728.*
                  pow7) - (4673.*x*x*x*
                  pow(log(x), 2.))/(13824.*pow6) + (17.*
                  x*x*x*x*pow(log(x), 2.))/(576.*pow7) + 
            log(Ms2/mumatch2)*((-12869.*x)/(41472.*pow6) + (85.*Nf*
                        x)/(1728.*pow6) - (697.*x)/(5184.*
                        pow5) - (7429.*x*x)/(20736.*
                        pow6) + (17.*Nf*x*x)/(384.*
                        pow6) + (41.*x*x)/(576.*
                        pow5) + (3859.*x*x*x)/(5184.*
                        pow6) - (17.*Nf*x*x*x)/(192.*
                        pow6) + (41.*x*x*x)/(576.*
                        pow5) - (2363.*x*x*x*x)/(20736.*
                        pow6) - (17.*Nf*x*x*x*x)/(3456.*
                        pow6) - (41.*x*x*x*x)/(5184.*
                        pow5) + (1853.*x*x*x*x*x)/(41472.*
                        pow6) - (17.*x*x*x*x*x*x)/(2592.*
                        pow6) - (595.*x*log(x))/(6912.*
                        pow6) + (17.*Nf*x*log(x))/(1152.*
                        pow6) - (41.*x*log(x))/(864.*
                        pow5) - (2159.*x*x*log(x))/(3456.*
                        pow6) + (17.*Nf*x*x*log(x))/(192.*
                        pow6) - (41.*x*x*log(x))/(288.*
                        pow5) - (1547.*x*x*x*log(x))/(6912.*
                        pow6) + (17.*Nf*x*x*x*log(x))/(384.*
                        pow6))))   
            );
}

double SUSYMassInsertionMatching::C2NLO(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    
    return (64.*(((17.*x)/(6912.*pow7) - (10861.*x)/(248832.*
                  pow6) - (19.*Nf*x)/(2304.*
                  pow6) - (18253.*x)/(248832.*
                  pow5) - (77.*x*x)/(6912.*
                  pow7) - (19159.*x*x)/(124416.*
                  pow6) + (11.*Nf*x*x)/(256.*
                  pow6) + (4289.*x*x)/(82944.*
                  pow5) + (13.*x*x*x)/(1152.*
                  pow7) + (659.*x*x*x)/(3456.*
                  pow6) - (9.*Nf*x*x*x)/(256.*
                  pow6) + (77.*x*x*x)/(3072.*
                  pow5) + (5.*x*x*x*x)/(3456.*
                  pow7) + (287.*x*x*x*x)/(124416.*
                  pow6) + (Nf*x*x*x*x)/(2304.*
                  pow6) - (851.*x*x*x*x)/(248832.*
                  pow5) - (31.*x*x*x*x*x)/(6912.*
                  pow7) + (1445.*x*x*x*x*x)/(248832.*
                  pow6) + x*x*x*x*x*x/(2304.*pow7) - 
            x*x*x*x*x*x/(864.*pow6) - (35.*x*Li2)/(3456.*
                  pow6) + (Nf*x*Li2)/(128.*
                  pow6) - (7.*x*Li2)/(432.*
                  pow5) - (13.*x*x*Li2)/(1728.*
                  pow6) - (Nf*x*x*Li2)/(64.*
                  pow6) - (5.*x*x*Li2)/(864.*
                  pow5) + (61.*x*x*x*Li2)/(3456.*
                  pow6) + (Nf*x*x*x*Li2)/(128.*
                  pow6) + (731.*x*a)/(31104.*
                  pow5) - (43.*x*x*a)/(3456.*
                  pow5) - (43.*x*x*x*a)/(3456.*
                  pow5) + (43.*x*x*x*x*a)/(31104.*
                  pow5) + (x*log(x))/(1152.*pow7) + (35.*
                  x*log(x))/(1152.*pow6) + (Nf*x*log(x))/(768.*
                  pow6) - (401.*x*log(x))/(41472.*
                  pow5) - (x*x*log(x))/(1152.*
                  pow7) - (53.*x*x*log(x))/(648.*
                  pow6) + (Nf*x*x*log(x))/(384.*
                  pow6) - (3605.*x*x*log(x))/(41472.*
                  pow5) - (11.*x*x*x*log(x))/(864.*
                  pow7) - (533.*x*x*x*log(x))/(2592.*
                  pow6) + (17.*Nf*x*x*x*log(x))/(768.*
                  pow6) - (11.*x*x*x*log(x))/(576.*
                  pow5) + (x*x*x*x*log(x))/(96.*
                  pow7) + (229.*x*x*x*x*log(x))/(20736.*
                  pow6) + (11.*x*x*x*x*log(x))/(5184.*
                  pow5) + (x*x*x*x*x*log(x))/(384.*
                  pow7) - (143.*x*x*x*x*x*log(x))/(20736.*
                  pow6) - (x*x*x*x*x*x*log(x))/(3456.*
                  pow7) + (x*x*x*x*x*x*log(x))/(864.*
                  pow6) + (43.*x*a*log(x))/(5184.*
                  pow5) + (43.*x*x*a*log(x))/(1728.*
                  pow5) + (7.*x*pow(log(x), 2.))/(512.*
                  pow6) + (11.*x*pow(log(x), 2.))/(1728.*
                  pow5) + (215.*x*x*
                  pow(log(x), 2.))/(2304.*pow6) + (35.*x*x*
                  pow(log(x), 2.))/(864.*pow5) - (x*x*x*
                  pow(log(x), 2.))/(576.*pow7) + (257.*x*x*x*
                  pow(log(x), 2.))/(4608.*pow6) - (x*x*x*x*
                  pow(log(x), 2.))/(192.*pow7) + 
            log(Ms2/mumatch2)*((757.*x)/(13824.*pow6) - (5.*Nf*x)/(576.*
                        pow6) + (187.*x)/(5184.*
                        pow5) + (437.*x*x)/(6912.*
                        pow6) - (Nf*x*x)/(128.*
                        pow6) - (11.*x*x)/(576.*
                        pow5) - (227.*x*x*x)/(1728.*
                        pow6) + (Nf*x*x*x)/(64.*
                        pow6) - (11.*x*x*x)/(576.*
                        pow5) + (139.*x*x*x*x)/(6912.*
                        pow6) + (Nf*x*x*x*x)/(1152.*
                        pow6) + (11.*x*x*x*x)/(5184.*
                        pow5) - (109.*x*x*x*x*x)/(13824.*
                        pow6) + 
                  x*x*x*x*x*x/(864.*pow6) + (35.*x*log(x))/(2304.*
                        pow6) - (Nf*x*log(x))/(384.*
                        pow6) + (11.*x*log(x))/(864.*
                        pow5) + (127.*x*x*log(x))/(1152.*
                        pow6) - (Nf*x*x*log(x))/(64.*
                        pow6) + (11.*x*x*log(x))/(288.*
                        pow5) + (91.*x*x*x*log(x))/(2304.*
                        pow6) - (Nf*x*x*x*log(x))/(128.*
                        pow6))))

            
            );
}

double SUSYMassInsertionMatching::C3NLOA(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    double pow4 = pow(1.-x,4.);
    
    return (64.*(11./(10368.*pow7) - 
                  8047./(248832.*pow6) - (209.*Nf)/(13824.*
                        pow6) - 15881./(248832.*pow5) - 
                  3./(256.*pow4) + (55.*x)/(10368.*
                        pow7) - (16187.*x)/(62208.*
                        pow6) + (121.*Nf*x)/(3456.*
                        pow6) - (38227.*x)/(82944.*
                        pow5) + (35.*x)/(1536.*
                        pow4) - (77.*x*x)/(1728.*
                        pow7) - (6361.*x*x)/(13824.*
                        pow6) + (55.*Nf*x*x)/(768.*
                        pow6) + (43003.*x*x)/(82944.*
                        pow5) - (17.*x*x)/(1536.*
                        pow4) + (341.*x*x*x)/(5184.*
                        pow7) + (56701.*x*x*x)/(62208.*
                        pow6) - (319.*Nf*x*x*x)/(3456.*
                        pow6) + (1553.*x*x*x)/(248832.*
                        pow5) - (253.*x*x*x*x)/(10368.*
                        pow7) - (40567.*x*x*x*x)/(248832.*
                        pow6) + (11.*Nf*x*x*x*x)/(13824.*
                        pow6) - (11.*x*x*x*x*x)/(3456.*
                        pow7) + (11.*x*x*x*x*x)/(2592.*
                        pow6) - (11.*Li2)/(6912.*
                        pow6) + (11.*Nf*Li2)/(1152.*
                        pow6) - (5.*x*Li2)/(108.*
                        pow6) + (61.*x*Li2)/(1728.*
                        pow5) + (x*Li2)/(256.*
                        pow4) - (31.*x*x*Li2)/(6912.*
                        pow6) - (11.*Nf*x*x*
                        Li2)/(384.*pow6) + (125.*x*x*
                        Li2)/(864.*pow5) + (181.*x*x*x*
                        Li2)/(3456.*pow6) + (11.*Nf*
                        x*x*x*Li2)/(576.*pow6) - (17.*
                        a)/(10368.*pow5) - (17.*x*
                        a)/(1152.*pow5) + (17.*x*x*
                        a)/(1152.*pow5) + (17.*x*x*x*
                        a)/(10368.*pow5) + (11.*x*
                        log(x))/(1728.*pow7) + (83.*x*
                        log(x))/(1152.*pow6) - (47.*x*log(x))/(432.*
                        pow5) + (5.*x*log(x))/(1536.*
                        pow4) - (55.*x*x*log(x))/(2592.*
                        pow7) - (1037.*x*x*log(x))/(1728.*
                        pow6) + (11.*Nf*x*x*log(x))/(288.*
                        pow6) - (455.*x*x*log(x))/(1536.*
                        pow5) - (11.*x*x*x*log(x))/(432.*
                        pow7) - (293.*x*x*x*log(x))/(648.*
                        pow6) + (11.*Nf*x*x*x*log(x))/(192.*
                        pow6) - (5.*x*x*x*log(x))/(512.*
                        pow5) + (11.*x*x*x*x*log(x))/(288.*
                        pow7) + (637.*x*x*x*x*log(x))/(10368.*
                        pow6) + (11.*x*x*x*x*x*log(x))/(5184.*
                        pow7) - (11.*x*x*x*x*x*log(x))/(2592.*
                        pow6) - (17.*x*a*log(x))/(1728.*
                        pow5) - (17.*x*x*a*
                        log(x))/(1728.*pow5) + (19.*x*
                        pow(log(x), 2.))/(192.*pow6) + (23.*x*
                        pow(log(x), 2.))/(256.*pow5) + (19.*x*
                        pow(log(x), 2.))/(1536.*pow4) + (337.*
                        x*x*pow(log(x), 2.))/(1152.*
                        pow6) + (45.*x*x*
                        pow(log(x), 2.))/(256.*pow5) - (11.*
                        x*x*x*pow(log(x), 2.))/(864.*
                        pow7) + (589.*x*x*x*
                        pow(log(x), 2.))/(3456.*pow6) - (11.*
                        x*x*x*x*pow(log(x), 2.))/(864.*
                        pow7) + 
                  log(Ms2/mumatch2)*(11./(2592.*pow6) - (11.*Nf)/(6912.*
                              pow6) + 
                        65./(13824.*pow5) + (1595.*x)/(5184.*
                              pow6) - (77.*Nf*x)/(1728.*
                              pow6) + (65.*x)/(1536.*
                              pow5) - (55.*x*x)/(10368.*
                              pow6) - (65.*x*x)/(1536.*
                              pow5) - (1705.*x*x*x)/(5184.*
                              pow6) + (77.*Nf*x*x*x)/(1728.*
                              pow6) - (65.*x*x*x)/(13824.*
                              pow5) + (275.*x*x*x*x)/(10368.*
                              pow6) + (11.*Nf*x*x*x*x)/(6912.*
                              pow6) - (11.*x*x*x*x*x)/(2592.*
                              pow6) + (385.*x*log(x))/(3456.*
                              pow6) - (11.*Nf*x*log(x))/(576.*
                              pow6) + (65.*x*log(x))/(2304.*
                              pow5) + (715.*x*x*
                              log(x))/(1728.*pow6) - (11.*Nf*
                              x*x*log(x))/(192.*
                              pow6) + (65.*x*x*
                              log(x))/(2304.*pow5) + (275.*
                              x*x*x*log(x))/(3456.*
                              pow6) - (11.*Nf*x*x*x*
                              log(x))/(576.*pow6)))
            );
}

double SUSYMassInsertionMatching::C3NLOB(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    double pow4 = pow(1.-x,4.);
    
    return (64.*
           
            (((5./(3456.*pow7) - 
                  1465./(82944.*pow6) - (95.*Nf)/(4608.*
                        pow6) - 4151./(82944.*pow5) - 
                  5./(768.*pow4) + (25.*x)/(3456.*
                        pow7) - (11465.*x)/(20736.*
                        pow6) + (55.*Nf*x)/(1152.*
                        pow6) - (12493.*x)/(27648.*
                        pow5) + (5.*x)/(512.*pow4) - (35.*
                        pow(x, 2.))/(576.*pow7) - (535.*
                        pow(x, 2.))/(4608.*pow6) + (25.*Nf*
                        pow(x, 2.))/(256.*pow6) + (13093.*
                        pow(x, 2.))/(27648.*pow5) - (5.*
                        pow(x, 2.))/(1536.*pow4) + (155.*
                        pow(x, 3.))/(1728.*pow7) + (14695.*
                        pow(x, 3.))/(20736.*pow6) - (145.*Nf*
                        pow(x, 3.))/(1152.*pow6) + (2351.*
                        pow(x, 3.))/(82944.*pow5) - (115.*
                        pow(x, 4.))/(3456.*pow7) - (2305.*
                        pow(x, 4.))/(82944.*pow6) + (5.*Nf*
                        pow(x, 4.))/(4608.*pow6) - (5.*
                        pow(x, 5.))/(1152.*pow7) + (5.*
                        pow(x, 5.))/(864.*pow6) - (5.*
                        Li2)/(2304.*pow6) + (5.*Nf*
                        Li2)/(384.*pow6) - (5.*x*
                        Li2)/(144.*pow6) - (95.*x*
                        Li2)/(576.*pow5) + (5.*x*
                        Li2)/(256.*pow4) - (145.*pow(x, 2.)*
                        Li2)/(2304.*pow6) - (5.*Nf*
                        pow(x, 2.)*Li2)/(128.*pow6) + (5.*
                        pow(x, 2.)*Li2)/(72.*pow5) + (115.*
                        pow(x, 3.)*Li2)/(1152.*pow6) + (5.*
                        Nf*pow(x, 3.)*Li2)/(192.*pow6) - 
                  a/(1152.*pow5) - (x*a)/(128.*
                        pow5) + (pow(x, 2.)*a)/(128.*
                        pow5) + (pow(x, 3.)*a)/(1152.*
                        pow5) + (5.*x*log(x))/(576.*
                        pow7) + (5.*x*log(x))/(384.*
                        pow6) - (139.*x*log(x))/(576.*
                        pow5) + (25.*x*log(x))/(1536.*
                        pow4) - (25.*pow(x, 2.)*log(x))/(864.*
                        pow7) - (95.*pow(x, 2.)*log(x))/(144.*
                        pow6) + (5.*Nf*pow(x, 2.)*log(x))/(96.*
                        pow6) - (1897.*pow(x, 2.)*log(x))/(4608.*
                        pow5) - (5.*pow(x, 3.)*log(x))/(144.*
                        pow7) - (545.*pow(x, 3.)*log(x))/(864.*
                        pow6) + (5.*Nf*pow(x, 3.)*log(x))/(64.*
                        pow6) - (35.*pow(x, 3.)*log(x))/(1536.*
                        pow5) + (5.*pow(x, 4.)*log(x))/(96.*
                        pow7) + (85.*pow(x, 4.)*log(x))/(3456.*
                        pow6) + (5.*pow(x, 5.)*log(x))/(1728.*
                        pow7) - (5.*pow(x, 5.)*log(x))/(864.*
                        pow6) - (x*a*log(x))/(192.*
                        pow5) - (pow(x, 2.)*a*log(x))/(192.*
                        pow5) + (5.*x*pow(log(x), 2.))/(64.*
                        pow6) + (15.*x*pow(log(x), 2.))/(256.*
                        pow5) + (5.*x*pow(log(x), 2.))/(512.*
                        pow4) + (175.*pow(x, 2.)*
                        pow(log(x), 2.))/(384.*pow6) + (45.*
                        pow(x, 2.)*pow(log(x), 2.))/(256.*
                        pow5) - (5.*pow(x, 3.)*
                        pow(log(x), 2.))/(288.*pow7) + (235.*
                        pow(x, 3.)*pow(log(x), 2.))/(1152.*
                        pow6) - (5.*pow(x, 4.)*
                        pow(log(x), 2.))/(288.*pow7) + 
                  log(Ms2/mumatch2)*(5./(864.*pow6) - (5.*Nf)/(2304.*
                              pow6) + 
                        95./(4608.*pow5) + (725.*x)/(1728.*
                              pow6) - (35.*Nf*x)/(576.*
                              pow6) + (95.*x)/(512.*
                              pow5) - (25.*pow(x, 2.))/(3456.*
                              pow6) - (95.*pow(x, 2.))/(512.*
                              pow5) - (775.*pow(x, 3.))/(1728.*
                              pow6) + (35.*Nf*pow(x, 3.))/(576.*
                              pow6) - (95.*pow(x, 3.))/(4608.*
                              pow5) + (125.*pow(x, 4.))/(3456.*
                              pow6) + (5.*Nf*pow(x, 4.))/(2304.*
                              pow6) - (5.*pow(x, 5.))/(864.*
                              pow6) + (175.*x*log(x))/(1152.*
                              pow6) - (5.*Nf*x*log(x))/(192.*
                              pow6) + (95.*x*log(x))/(768.*
                              pow5) + (325.*pow(x, 2.)*
                              log(x))/(576.*pow6) - (5.*Nf*
                              pow(x, 2.)*log(x))/(64.*pow6) + (95.*
                              pow(x, 2.)*log(x))/(768.*
                              pow5) + (125.*pow(x, 3.)*
                              log(x))/(1152.*pow6) - (5.*Nf*
                              pow(x, 3.)*log(x))/(192.*pow6)))))
            
            );    
}

double SUSYMassInsertionMatching::C4NLOA(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1.-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    double pow4 = pow(1.-x,4.);
    
    return (
             (5./(3456.*pow7) - 
                  1465./(82944.*pow6) - (95.*Nf)/(4608.*
                        pow6) - 4151./(82944.*pow5) - 
                  5./(768.*pow4) + (25.*x)/(3456.*
                        pow7) - (11465.*x)/(20736.*
                        pow6) + (55.*Nf*x)/(1152.*
                        pow6) - (12493.*x)/(27648.*
                        pow5) + (5.*x)/(512.*pow4) - (35.*
                        x*x)/(576.*pow7) - (535.*
                        x*x)/(4608.*pow6) + (25.*Nf*
                        x*x)/(256.*pow6) + (13093.*
                        x*x)/(27648.*pow5) - (5.*
                        x*x)/(1536.*pow4) + (155.*
                        x*x*x)/(1728.*pow7) + (14695.*
                        x*x*x)/(20736.*pow6) - (145.*Nf*
                        x*x*x)/(1152.*pow6) + (2351.*
                        x*x*x)/(82944.*pow5) - (115.*
                        x*x*x*x)/(3456.*pow7) - (2305.*
                        x*x*x*x)/(82944.*pow6) + (5.*Nf*
                        x*x*x*x)/(4608.*pow6) - (5.*
                        x*x*x*x*x)/(1152.*pow7) + (5.*
                        x*x*x*x*x)/(864.*pow6) - (5.*
                        Li2)/(2304.*pow6) + (5.*Nf*
                        Li2)/(384.*pow6) - (5.*x*
                        Li2)/(144.*pow6) - (95.*x*
                        Li2)/(576.*pow5) + (5.*x*
                        Li2)/(256.*pow4) - (145.*x*x*
                        Li2)/(2304.*pow6) - (5.*Nf*
                        x*x*Li2)/(128.*pow6) + (5.*
                        x*x*Li2)/(72.*pow5) + (115.*
                        x*x*x*Li2)/(1152.*pow6) + (5.*
                        Nf*x*x*x*Li2)/(192.*pow6) - 
                  a/(1152.*pow5) - (x*a)/(128.*
                        pow5) + (x*x*a)/(128.*
                        pow5) + (x*x*x*a)/(1152.*
                        pow5) + (5.*x*log(x))/(576.*
                        pow7) + (5.*x*log(x))/(384.*
                        pow6) - (139.*x*log(x))/(576.*
                        pow5) + (25.*x*log(x))/(1536.*
                        pow4) - (25.*x*x*log(x))/(864.*
                        pow7) - (95.*x*x*log(x))/(144.*
                        pow6) + (5.*Nf*x*x*log(x))/(96.*
                        pow6) - (1897.*x*x*log(x))/(4608.*
                        pow5) - (5.*x*x*x*log(x))/(144.*
                        pow7) - (545.*x*x*x*log(x))/(864.*
                        pow6) + (5.*Nf*x*x*x*log(x))/(64.*
                        pow6) - (35.*x*x*x*log(x))/(1536.*
                        pow5) + (5.*x*x*x*x*log(x))/(96.*
                        pow7) + (85.*x*x*x*x*log(x))/(3456.*
                        pow6) + (5.*x*x*x*x*x*log(x))/(1728.*
                        pow7) - (5.*x*x*x*x*x*log(x))/(864.*
                        pow6) - (x*a*log(x))/(192.*
                        pow5) - (x*x*a*log(x))/(192.*
                        pow5) + (5.*x*pow(log(x), 2.))/(64.*
                        pow6) + (15.*x*pow(log(x), 2.))/(256.*
                        pow5) + (5.*x*pow(log(x), 2.))/(512.*
                        pow4) + (175.*x*x*
                        pow(log(x), 2.))/(384.*pow6) + (45.*
                        x*x*pow(log(x), 2.))/(256.*
                        pow5) - (5.*x*x*x*
                        pow(log(x), 2.))/(288.*pow7) + (235.*
                        x*x*x*pow(log(x), 2.))/(1152.*
                        pow6) - (5.*x*x*x*x*
                        pow(log(x), 2.))/(288.*pow7) + 
                  log(Ms2/mumatch2)*(5./(864.*pow6) - (5.*Nf)/(2304.*
                              pow6) + 
                        95./(4608.*pow5) + (725.*x)/(1728.*
                              pow6) - (35.*Nf*x)/(576.*
                              pow6) + (95.*x)/(512.*
                              pow5) - (25.*x*x)/(3456.*
                              pow6) - (95.*x*x)/(512.*
                              pow5) - (775.*x*x*x)/(1728.*
                              pow6) + (35.*Nf*x*x*x)/(576.*
                              pow6) - (95.*x*x*x)/(4608.*
                              pow5) + (125.*x*x*x*x)/(3456.*
                              pow6) + (5.*Nf*x*x*x*x)/(2304.*
                              pow6) - (5.*x*x*x*x*x)/(864.*
                              pow6) + (175.*x*log(x))/(1152.*
                              pow6) - (5.*Nf*x*log(x))/(192.*
                              pow6) + (95.*x*log(x))/(768.*
                              pow5) + (325.*x*x*
                              log(x))/(576.*pow6) - (5.*Nf*
                              x*x*log(x))/(64.*pow6) + (95.*
                              x*x*log(x))/(768.*
                              pow5) + (125.*x*x*x*
                              log(x))/(1152.*pow6) - (5.*Nf*
                              x*x*x*log(x))/(192.*pow6)))
            );   
}

double SUSYMassInsertionMatching::C4NLOB(double x, double mumatch2, double Ms2){
    double Li2 = gsl_sf_dilog(1.-x);
    double a = log(2.);
    double pow7 = pow(1.-x,7.);
    double pow6 = pow(1.-x,6.);
    double pow5 = pow(1.-x,5.);
    double pow4 = pow(1.-x,4.);
    
    return (64.*(-5./(5184.*pow7) + 
                  149./(1944.*pow6) + (95.*Nf)/(6912.*
                        pow6) - 3653./(124416.*pow5) + 
                  29./(10368.*pow4) - (67.*x)/(10368.*
                        pow7) + (3971.*x)/(15552.*
                        pow6) - (91.*Nf*x)/(3456.*
                        pow6) + (42647.*x)/(82944.*
                        pow5) - (17.*x)/(3456.*
                        pow4) + (497.*x*x)/(10368.*
                        pow7) + (2059.*x*x)/(6912.*
                        pow6) - (3.*Nf*x*x)/(32.*
                        pow6) - (1397.*x*x)/(3072.*
                        pow5) + (11.*x*x)/(5184.*
                        pow4) - (349.*x*x*x)/(5184.*
                        pow7) - (41527.*x*x*x)/(62208.*
                        pow6) + (371.*Nf*x*x*x)/(3456.*
                        pow6) - (8525.*x*x*x)/(248832.*
                        pow5) + (55.*x*x*x*x)/(2592.*
                        pow7) + (2773.*x*x*x*x)/(62208.*
                        pow6) - (7.*Nf*x*x*x*x)/(6912.*
                        pow6) + (349.*x*x*x*x)/(82944.*
                        pow5) + (61.*x*x*x*x*x)/(10368.*
                        pow7) - (53.*x*x*x*x*x)/(6912.*
                        pow6) - 
                  x*x*x*x*x*x/(3456.*pow7) + 
                  x*x*x*x*x*x/(1296.*pow6) - (19.*Li2)/(432.*
                        pow6) - (5.*Nf*Li2)/(576.*
                        pow6) + (41.*Li2)/(1152.*
                        pow5) + (599.*x*Li2)/(3456.*
                        pow6) - (Nf*x*Li2)/(192.*
                        pow6) + (13.*x*Li2)/(576.*
                        pow5) - (7.*x*Li2)/(1728.*
                        pow4) - (179.*x*x*
                        Li2)/(1728.*pow6) + (7.*Nf*
                        x*x*Li2)/(192.*pow6) + (7.*
                        x*x*Li2)/(384.*pow5) - (89.*
                        x*x*x*Li2)/(3456.*pow6) - (13.*
                        Nf*x*x*x*Li2)/(576.*
                        pow6) - (7.*a)/(5184.*
                        pow5) + (401.*x*a)/(10368.*
                        pow5) - (17.*x*x*a)/(1152.*
                        pow5) - (265.*x*x*x*a)/(10368.*
                        pow5) + (31.*x*x*x*x*a)/(10368.*
                        pow5) - (11.*x*log(x))/(1728.*
                        pow7) + (41.*x*log(x))/(768.*
                        pow6) - (Nf*x*log(x))/(1152.*
                        pow6) + (5075.*x*log(x))/(41472.*
                        pow5) - (35.*x*log(x))/(10368.*
                        pow4) + (103.*x*x*log(x))/(5184.*
                        pow7) + (6985.*x*x*log(x))/(20736.*
                        pow6) - (7.*Nf*x*x*log(x))/(192.*
                        pow6) + (5741.*x*x*log(x))/(13824.*
                        pow5) + (41.*x*x*x*log(x))/(1296.*
                        pow7) + (4411.*x*x*x*log(x))/(6912.*
                        pow6) - (77.*Nf*x*x*x*log(x))/(1152.*
                        pow6) + (443.*x*x*x*log(x))/(20736.*
                        pow5) - (x*x*x*x*log(x))/(24.*
                        pow7) - (67.*x*x*x*x*log(x))/(2304.*
                        pow6) - (7.*x*x*x*x*log(x))/(5184.*
                        pow5) - (19.*x*x*x*x*x*log(x))/(5184.*
                        pow7) + (29.*x*x*x*x*x*log(x))/(3456.*
                        pow6) + (x*x*x*x*x*x*log(x))/(5184.*
                        pow7) - (x*x*x*x*x*x*log(x))/(1296.*
                        pow6) + (17.*x*a*log(x))/(1728.*
                        pow5) + (79.*x*x*a*
                        log(x))/(1728.*pow5) - (71.*x*
                        pow(log(x), 2.))/(1536.*pow6) - (89.*x*
                        pow(log(x), 2.))/(2304.*pow5) - (x*
                        pow(log(x), 2.))/(288.*pow4) - (1007.*
                        x*x*pow(log(x), 2.))/(2304.*
                        pow6) - (281.*x*x*
                        pow(log(x), 2.))/(2304.*pow5) + (11.*
                        x*x*x*pow(log(x), 2.))/(864.*
                        pow7) - (2551.*x*x*x*
                        pow(log(x), 2.))/(13824.*pow6) + (13.*
                        x*x*x*x*pow(log(x), 2.))/(864.*
                        pow7) + 
                  log(Ms2/mumatch2)*((5.*Nf)/(3456.*pow6) - 
                        95./(6912.*pow5) - (1927.*x)/(6912.*
                              pow6) + (5.*Nf*x)/(108.*
                              pow6) - (2033.*x)/(13824.*
                              pow5) - (1211.*x*x)/(10368.*
                              pow6) + (Nf*x*x)/(192.*
                              pow6) + (209.*x*x)/(1536.*
                              pow5) + (541.*x*x*x)/(1296.*
                              pow6) - (11.*Nf*x*x*x)/(216.*
                              pow6) + (361.*x*x*x)/(13824.*
                              pow5) - (103.*x*x*x*x)/(3456.*
                              pow6) - (7.*Nf*x*x*x*x)/(3456.*
                              pow6) - (19.*x*x*x*x)/(13824.*
                              pow5) + (181.*x*x*x*x*x)/(20736.*
                              pow6) - 
                        x*x*x*x*x*x/(1296.*pow6) - (11.*x*
                              log(x))/(128.*pow6) + (11.*Nf*x*
                              log(x))/(576.*pow6) - (209.*x*
                              log(x))/(2304.*pow5) - (769.*
                              x*x*log(x))/(1728.*
                              pow6) + (Nf*x*x*log(x))/(16.*
                              pow6) - (247.*x*x*
                              log(x))/(2304.*pow5) - (445.*
                              x*x*x*log(x))/(3456.*
                              pow6) + (13.*Nf*x*x*x*
                              log(x))/(576.*pow6)))
            
            );
}

 std::vector<WilsonCoefficient>& SUSYMassInsertionMatching::CMd1() {
    
    vmcd1.clear();
    vmcd1 = StandardModelMatching::CMd1();
    
    double x = pow(SusyMI.getM3() / SusyMI.getMsq(), 2.);
    double constLO = SusyMI.Als(SusyMI.getMuM()) / SusyMI.getMsq();
    double mcharm = SusyMI.getMuc();
    double mgluino = SusyMI.getM3();
    
    DLL = SusyMI.getDu_LL()(0,1);
    DLR = SusyMI.getDu_LR()(0,1);
    DRL = SusyMI.getDu_RL()(0,1);
    DRR = SusyMI.getDu_RR()(0,1);
    /*
    std::cout << " DLL --> " << DLL << std::endl;
    std::cout << " DLR --> " << DLR << std::endl;
    std::cout << " DRL --> " << DRL << std::endl;
    std::cout << " DRR --> " << DRR << std::endl;
    */
    
    switch (mcd1.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcd1.getScheme();
            throw std::runtime_error("SUSYMassInsertionMatching::CMD1(): scheme " + out.str() + "not implemented"); 
    }

    mcd1.setMu(SusyMI.getMuM());
    
    gslpp::vector<gslpp::complex> C0_B(10,0.), C0_M(10,0.), C0_M_eff(10,0.);
    
    C0_B.assign(2, constLO * constLO / 4. * (-1./9.*B1(x) - 5./9.*B2(x) 
        - 1./18.*P1(x) - 1./2.*P2(x))* DLL);
    C0_B.assign(3, constLO * constLO / 4. * (-7./3.*B1(x) + 1./3.*B2(x) 
        + 1./6.*P1(x) + 3./2.*P2(x))* DLL);
    C0_B.assign(4, constLO * constLO / 4. * (10./9.*B1(x) + 1./18.*B2(x) 
        - 1./18.*P1(x) - 1./2.*P2(x))* DLL);
    C0_B.assign(5, constLO * constLO / 4. * (-2./3.*B1(x) + 7./6.*B2(x) 
        + 1./6.*P1(x) + 3./2.*P2(x)) * DLL);
    C0_B.assign(6, constLO/SusyMI.getMsq() / 4. * M_PI * ((8./3.*M3(x)) * DLL
                + mgluino/mcharm * (8./3.*M1(x)) * DLR));
    C0_B.assign(7, constLO/SusyMI.getMsq() / 4. * M_PI * ((-1./3.*M3(x)-3.*M4(x)) * DLL 
            + mgluino/mcharm * (-1./3.*M1(x)-3.*M2(x)) * DLR));
    
    C0_M = RtoMisiak().transpose() * C0_B;
    C0_M_eff = EffectiveBase() * C0_M ;
    
  /*  for(int i= 0; i< 10; i ++){
        std::cout << "i" << "  "<< "C0_B(i)" << "  " <<  "C0_M" << "  " << "C0_M_eff" << std::endl;
        std::cout << i << "  "<< C0_B(i) << "  " <<  C0_M(i) << "  " << C0_M_eff(i) << std::endl;
    } */
    
     switch (mcd1.getOrder()) {
        case NLO: 
            for(int k = 0; k<10; k++){
            mcd1.setCoeff(k, 0., NLO);
            }
        case LO:
            for(int k = 0; k<10; k++){
            mcd1.setCoeff(k, C0_M_eff(k), LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcd1.getOrder();
            throw std::runtime_error("SUSYMassInsertionMatching::CMd1(): order " + out.str() + "not implemented"); 
    }

    vmcd1.push_back(mcd1);
    
    C0_B.assign(2, constLO * constLO / 4. * (-1./9.*B1(x) - 5./9.*B2(x) 
        - 1./18.*P1(x) - 1./2.*P2(x))* DRR);
    C0_B.assign(3, constLO * constLO / 4. * (-7./3.*B1(x) + 1./3.*B2(x) 
        + 1./6.*P1(x) + 3./2.*P2(x))* DRR);
    C0_B.assign(4, constLO * constLO / 4. * (10./9.*B1(x) + 1./18.*B2(x) 
        - 1./18.*P1(x) - 1./2.*P2(x))* DRR);
    C0_B.assign(5, constLO * constLO / 4. * (-2./3.*B1(x) + 7./6.*B2(x) 
        + 1./6.*P1(x) + 3./2.*P2(x)) * DRR);
    C0_B.assign(6, constLO/SusyMI.getMsq() / 4. * M_PI * ((8./3.*M3(x)) * DRR 
            + mgluino/mcharm * (8./3.*M1(x) * DRL)));
    C0_B.assign(7, (constLO/SusyMI.getMsq() / 4. * M_PI * ((-1./3.*M3(x)-3.*M4(x)) * DRR
             + mgluino/mcharm * (-1./3.*M1(x)-3.*M2(x))* DRL)));
    C0_M = RtoMisiak().transpose() * C0_B;
    C0_M_eff = EffectiveBase() * C0_M ;
    
  /*  for(int i= 0; i< 10; i ++){
        std::cout << "i" << "  "<< "C0_B(i)" << "  " <<  "C0_M" << "  " << "C0_M_eff" << std::endl;
        std::cout << i << "  "<< C0_B(i) << "  " <<  C0_M(i) << "  " << C0_M_eff(i) << std::endl;
    } */
    
    switch (mcd1.getOrder()) {
        case NLO:
            for(int k = 0; k<10; k++){
            mcd1.setCoeff(k, 0., NLO);
            }
        case LO:
            for(int k=0; k<10; k++){
             mcd1.setCoeff(k, C0_M_eff(k), LO);   
            }
            break;
        default:
            std::stringstream out;
            out << mcd1.getOrder();
            throw std::runtime_error("SUSYMassInsertionMatching::CMd1(): order " + out.str() + "not implemented"); 
    }

    vmcd1.push_back(mcd1);
    
    return(vmcd1);
}

 std::vector<WilsonCoefficient>& SUSYMassInsertionMatching::CMdd2 () {
    
    vmcd2.clear();
    vmcd2 = StandardModelMatching::CMdd2();
    
    double als = SusyMI.Als(SusyMI.getMuM())*(1.+SusyMI.Als(SusyMI.getM3())/4./M_PI);
    double x = pow(SusyMI.getM3() / SusyMI.getMsq(), 2.);
    double coLO = pow(als/ SusyMI.getMsq(), 2.);
    double coNLO = coLO* als / M_PI;
    MuM2 = SusyMI.getMuM()*SusyMI.getMuM();
    Ms2 = SusyMI.getMsq()*SusyMI.getMsq();
    DLL = SusyMI.getDu_LL()(0,1);
    DLR = SusyMI.getDu_LR()(0,1);
    DRL = SusyMI.getDu_RL()(0,1);
    DRR = SusyMI.getDu_RR()(0,1);
    
    switch (mcd2.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcd2.getScheme();
            throw std::runtime_error("StandardModel::CMbd2(): scheme " + out.str() + "not implemented"); 
    }
    mcd2.setMu(SusyMI.getMuM());
     switch (mcd2.getOrder()) {    
        case NLO:
            mcd2.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DLL * DLL, NLO);
            mcd2.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mcd2.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mcd2.setCoeff(3, coNLO * (C3NLOA(x,MuM2,Ms2) * DLR * DRL +
                          C3NLOB(x,MuM2,Ms2) * DLL * DRR), NLO);
            mcd2.setCoeff(4, coNLO * (C4NLOA(x,MuM2,Ms2) * DLR * DRL + 
                          C4NLOB(x,MuM2,Ms2) * DLL * DRR), NLO);
        case LO:
            mcd2.setCoeff(0, coLO * C0LO(x) * DLL * DLL, LO);
            mcd2.setCoeff(1, coLO * C1LO(x) * DRL * DRL, LO);
            mcd2.setCoeff(2, coLO * C2LO(x) * DRL * DRL, LO);
            mcd2.setCoeff(3, coLO * (C3LOA(x) * DLR * DRL + C3LOB(x) * DLL * DRR), LO);
            mcd2.setCoeff(4, coLO * (C4LOA(x) * DLR * DRL + C4LOB(x) * DLL * DRR), LO);
            break;
        default:
            std::stringstream out;
            out << mcd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMd2(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(1);

    vmcd2.push_back(mcd2);
    
    switch (mcd2.getOrder()) {
        case NLO:
            mcd2.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DRR * DRR, NLO);
            mcd2.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mcd2.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mcd2.setCoeff(3, 0., NLO);
            mcd2.setCoeff(4, 0., NLO);
        case LO:
            mcd2.setCoeff(0, coLO * C0LO(x) * DRR * DRR, LO);
            mcd2.setCoeff(1, coLO * C1LO(x) * DLR * DLR, LO);
            mcd2.setCoeff(2, coLO * C2LO(x) * DLR * DLR, LO);
            mcd2.setCoeff(3, 0., LO);
            mcd2.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMd2(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(1);

    vmcd2.push_back(mcd2);
    
    return(vmcd2);
}

 std::vector<WilsonCoefficient>& SUSYMassInsertionMatching::CMdbd2() {
    
    vmcdb.clear();
    
    vmcdb = StandardModelMatching::CMdbd2();    
     
    double als = SusyMI.Als(SusyMI.getMuM())*(1.+SusyMI.Als(SusyMI.getM3())/4./M_PI);
    double x = pow(SusyMI.getM3() / SusyMI.getMsq(), 2.);
    double coLO = pow(als/ SusyMI.getMsq(), 2.);
    double coNLO = coLO * als / M_PI;
    MuM2 = SusyMI.getMuM()*SusyMI.getMuM();
    Ms2 = SusyMI.getMsq()*SusyMI.getMsq(); 
    
    DLL = SusyMI.getDd_LL()(0,2);
    DLR = SusyMI.getDd_LR()(0,2);
    DRL = SusyMI.getDd_RL()(0,2);
    DRR = SusyMI.getDd_RR()(0,2);
    
    switch (mcbd.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcbd.getScheme();
            throw std::runtime_error("StandardModel::CMbbd(): scheme " + out.str() + "not implemented"); 
    }

    mcbd.setMu(SusyMI.getMuM());
    
    switch (mcbd.getOrder()) {
        case NLO:
            mcbd.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DLL * DLL, NLO);
            mcbd.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mcbd.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mcbd.setCoeff(3, coNLO * (C3NLOA(x,MuM2,Ms2) * DLR * DRL + 
                          C3NLOB(x,MuM2,Ms2) * DLL * DRR), NLO);
            mcbd.setCoeff(4, coNLO * (C4NLOA(x,MuM2,Ms2) * DLR * DRL + 
                          C4NLOB(x,MuM2,Ms2) * DLL * DRR), NLO);
        case LO:
            mcbd.setCoeff(0, coLO * C0LO(x) * DLL * DLL, LO);
            mcbd.setCoeff(1, coLO * C1LO(x) * DRL * DRL, LO);
            mcbd.setCoeff(2, coLO * C2LO(x) * DRL * DRL, LO);
            mcbd.setCoeff(3, coLO * (C3LOA(x) * DLR * DRL + C3LOB(x) * DLL * DRR), LO);
            mcbd.setCoeff(4, coLO * (C4LOA(x) * DLR * DRL + C4LOB(x) * DLL * DRR), LO);
            break;
        default:
            std::stringstream out;
            out << mcbd.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbd(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(2);

    vmcdb.push_back(mcbd);
    
    switch (mcbd.getOrder()) {
        case NLO:
            mcbd.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DRR * DRR, NLO);
            mcbd.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mcbd.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mcbd.setCoeff(3, 0., NLO);
            mcbd.setCoeff(4, 0., NLO);
        case LO:
            mcbd.setCoeff(0, coLO * C0LO(x) * DRR * DRR, LO);
            mcbd.setCoeff(1, coLO * C1LO(x) * DLR * DLR, LO);
            mcbd.setCoeff(2, coLO * C2LO(x) * DLR * DLR, LO);
            mcbd.setCoeff(3, 0., LO);
            mcbd.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcbd.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbd(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(2);

    vmcdb.push_back(mcbd);
    
    return(vmcdb);
}

 std::vector<WilsonCoefficient>& SUSYMassInsertionMatching::CMdbs2() {
    
    
    vmcds.clear();
    
    vmcds = StandardModelMatching::CMdbs2();
    
    double als = SusyMI.Als(SusyMI.getMuM())*(1.+SusyMI.Als(SusyMI.getM3())/4./M_PI);
    double x = pow(SusyMI.getM3() / SusyMI.getMsq(), 2.);
    double coLO = pow(als / SusyMI.getMsq(), 2.);
    double coNLO = coLO * als / M_PI;
    MuM2 = SusyMI.getMuM()*SusyMI.getMuM();
    Ms2 = SusyMI.getMsq()*SusyMI.getMsq(); 
    
    DLL = SusyMI.getDd_LL()(1,2);
    DLR = SusyMI.getDd_LR()(1,2);
    DRL = SusyMI.getDd_RL()(1,2);
    DRR = SusyMI.getDd_RR()(1,2);
    
    switch (mcbs.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcbs.getScheme();
            throw std::runtime_error("StandardModel::CMbbs(): scheme " + out.str() + "not implemented"); 
    }

    mcbs.setMu(SusyMI.getMuM());
    
    switch (mcbs.getOrder()) {
        case NLO:
            mcbs.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DLL * DLL, NLO);
            mcbs.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mcbs.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mcbs.setCoeff(3, coNLO * (C3NLOA(x,MuM2,Ms2) * DLR * DRL 
                          + C3NLOB(x,MuM2,Ms2) * DLL * DRR), NLO);
            mcbs.setCoeff(4, coNLO * (C4NLOA(x,MuM2,Ms2) * DLR * DRL + 
                          C4NLOB(x,MuM2,Ms2) * DLL * DRR), NLO);
        case LO:
            mcbs.setCoeff(0, coLO * C0LO(x) * DLL * DLL, LO);
            mcbs.setCoeff(1, coLO * C1LO(x) * DRL * DRL, LO);
            mcbs.setCoeff(2, coLO * C2LO(x) * DRL * DRL, LO);
            mcbs.setCoeff(3, coLO * (C3LOA(x) * DLR * DRL + C3LOB(x) * DLL * DRR), LO);
            mcbs.setCoeff(4, coLO * (C4LOA(x) * DLR * DRL + C4LOB(x) * DLL * DRR), LO);
            break;
        default:
            std::stringstream out;
            out << mcbs.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbs(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(3);

    vmcds.push_back(mcbs);
    
    switch (mcbs.getOrder()) {
        case NLO:
            mcbs.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DRR * DRR, NLO);
            mcbs.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mcbs.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DLR * DRL, NLO);
            mcbs.setCoeff(3, 0., NLO);
            mcbs.setCoeff(4, 0., NLO);
        case LO:
            mcbs.setCoeff(0, coLO * C0LO(x) * DRR * DRR, LO);
            mcbs.setCoeff(1, coLO * C1LO(x) * DLR * DLR, LO);
            mcbs.setCoeff(2, coLO * C2LO(x) * DLR * DLR, LO);
            mcbs.setCoeff(3, 0., LO);
            mcbs.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcbs.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbs(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(3);

    vmcds.push_back(mcbs);
    
    return(vmcds);
}

 std::vector<WilsonCoefficient>& SUSYMassInsertionMatching::CMdk2() {
    
    vmck2.clear();
    vmck2 = StandardModelMatching::CMdk2();
    
    double als = SusyMI.Als(SusyMI.getMuM())*(1.+SusyMI.Als(SusyMI.getM3())/4./M_PI);
    double x = pow(SusyMI.getM3() / SusyMI.getMsq(), 2.);
    double coLO = pow(als/ SusyMI.getMsq(), 2.);
    double coNLO = coLO * als / M_PI;
    MuM2 = SusyMI.getMuM()*SusyMI.getMuM();
    Ms2 = SusyMI.getMsq()*SusyMI.getMsq();   
    
    DLL = SusyMI.getDd_LL()(0,1);
    DLR = SusyMI.getDd_LR()(0,1);
    DRL = SusyMI.getDd_RL()(0,1);
    DRR = SusyMI.getDd_RR()(0,1);
    
    switch (mck2.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mck2.getScheme();
            throw std::runtime_error("StandardModel::CMbk2(): scheme " + out.str() + "not implemented"); 
    }

    mck2.setMu(SusyMI.getMuM());
    
    switch (mck2.getOrder()) {
        case NLO:
            mck2.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DLL * DLL, NLO);
            mck2.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mck2.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DRL * DRL, NLO);
            mck2.setCoeff(3, coNLO * C3NLOA(x,MuM2,Ms2) * DLR * DRL +
                          coNLO * C3NLOB(x,MuM2,Ms2) * DLL * DRR, NLO);
            mck2.setCoeff(4, coNLO * C4NLOA(x,MuM2,Ms2) * DLR * DRL + 
                          coNLO * C4NLOB(x,MuM2,Ms2) * DLL * DRR, NLO);
        case LO:
            mck2.setCoeff(0, coLO * C0LO(x) * DLL * DLL, LO);
            mck2.setCoeff(1, coLO * C1LO(x) * DRL * DRL, LO);
            mck2.setCoeff(2, coLO * C2LO(x) * DRL * DRL, LO);
            mck2.setCoeff(3, coLO * C3LOA(x) * DLR * DRL + coLO * coLO * C3LOB(x) * DLL * DRR, LO);
            mck2.setCoeff(4, coLO * C4LOA(x) * DLR * DRL + coLO * coLO * C4LOB(x) * DLL * DRR, LO);
            break;
        default:
            std::stringstream out;
            out << mck2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMk2(): order " + out.str() + "not implemented"); 
    }
    
    LRItoNDR(4);

    vmck2.push_back(mck2);
    
    switch (mck2.getOrder()) {
        case NLO:
            mck2.setCoeff(0, coNLO * C0NLO(x,MuM2,Ms2) * DRR * DRR, NLO);
            mck2.setCoeff(1, coNLO * C1NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mck2.setCoeff(2, coNLO * C2NLO(x,MuM2,Ms2) * DLR * DLR, NLO);
            mck2.setCoeff(3, 0., NLO);
            mck2.setCoeff(4, 0., NLO);
        case LO:
            mck2.setCoeff(0, coLO * C0LO(x) * DRR * DRR, LO);
            mck2.setCoeff(1, coLO * C1LO(x) * DLR * DLR, LO);
            mck2.setCoeff(2, coLO * C2LO(x) * DLR * DLR, LO);
            mck2.setCoeff(3, 0., LO);
            mck2.setCoeff(4, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mck2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMk2(): order " + out.str() + "not implemented"); 
    }

    LRItoNDR(4);
    
    vmck2.push_back(mck2);
    
    return(vmck2);
}

gslpp::matrix<double> SUSYMassInsertionMatching::RtoMisiak() const{
    
    gslpp::matrix<double> R(10,0.);
    
    R(0,0) = 2.;
    R(0,1) = 1./3.;
    R(1,1) = 1.;
    R(2,2) = -1./3.;
    R(2,4) = 1./12.;
    R(3,2) = -1./9.;
    R(3,3) = -2./3.;
    R(3,4) = 1./36.;
    R(3,5) = 1./6.;
    R(4,2) = 4./3.; 
    R(4,4) = -1./12.;
    R(5,2) = 4./9.;
    R(5,3) = 8./3.;
    R(5,4) = -1./36.;
    R(5,5) = -1./6.;
    R(6,6) = 1.;
    R(7,7) = 1.;
    R(8,8) = 1.;
    R(9,9) = 1.;
    return(R);
    
}

gslpp::matrix<double> SUSYMassInsertionMatching::EffectiveBase() const{

    gslpp::matrix<double> y(10, 0.);
    
    y(0,0) = 1.;
    y(1,1) = 1.;
    y(2,2) = 1.;
    y(3,3) = 1.;
    y(4,4) = 1.;
    y(5,5) = 1.;
    y(6,6) = 1.;
    y(7,7) = 1.;
    y(8,8) = 1.;
    y(9,9) = 1.;
    
    y(6,2) = -1./3.;
    y(6,3) = -4./9.;
    y(6,4) = -20./3.;
    y(6,5) = -80./9.;
    
    y(7,2) = 1.;
    y(7,3) = -1./6.;
    y(7,4) = 20.;
    y(7,5) = -10./3.;
 
    return(y);

}

void SUSYMassInsertionMatching::LRItoNDR(int i){
    switch (i){
        case 1:
            mcd2.setCoeff(*mcd2.getCoeff(NLO) + SusyMI.Als(mcd2.getMu()) / 4. / M_PI * 
                drNDRLRI.transpose() * (*mcd2.getCoeff(LO)), NLO);
            break;
        case 2:
            mcbd.setCoeff(*mcbd.getCoeff(NLO) + SusyMI.Als(mcbd.getMu()) / 4. / M_PI * 
                drNDRLRI.transpose() * (*mcbd.getCoeff(LO)), NLO);
            break;
        case 3:
            mcbs.setCoeff(*mcbs.getCoeff(NLO) + SusyMI.Als(mcbs.getMu()) / 4. / M_PI * 
                drNDRLRI.transpose() * (*mcbs.getCoeff(LO)), NLO);
            break;
        case 4:
            mck2.setCoeff(*mck2.getCoeff(NLO) + SusyMI.Als(mck2.getMu()) / 4. / M_PI * 
                drNDRLRI.transpose() * (*mck2.getCoeff(LO)), NLO);
            break;
        default:
            throw std::runtime_error("SUSYMassInsertionMatching::LRItoNDR : change of scheme not implemented"); 
    }
}
