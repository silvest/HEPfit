/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModelMatching.h"
#include "StandardModel.h"
#include "QCD.h"
#include <stdexcept>
#include <sstream>

#define ZEW_NUMERIC

StandardModelMatching::StandardModelMatching(const StandardModel & SM_i) 
: SM(SM_i),
        mcdbd2(5, NDR, NLO),
        mcdbs2(5, NDR, NLO),
        mcdd2(5, NDR, NLO),
        mcdk2(5, NDR, NLO),
        mck(10, NDR, NLO),
        mckcc(10, NDR, NLO),
        mcbsg(8, NDR, NNLO),
        mcprimebsg(8, NDR, NNLO),
        mcBMll(13, NDR, NLO),
        mcprimeBMll(13, NDR, NLO),
        mcbnlep(10, NDR, NLO, NLO_QED11),
        mcbnlepCC(10, NDR, NLO),
        mcd1(10, NDR, NLO),
        mcd1Buras(10, NDR, NLO),
        mckpnn(1, NDR, NLO, NLO_QED11),
        mckmm(1, NDR, NLO),
        mcbsnn(1, NDR, NLO),
        mcbdnn(1, NDR, NLO),
        mcbsmm(8, NDR, NNLO, NLO_QED22),
        mcbdmm(8, NDR, NNLO, NLO_QED22),
        mcbtaunu(3, NDR, LO),
        mcDLij(2, NDR, LO),
        mcDLi3j(20, NDR, LO),
        mcmueconv(8, NDR, LO),
        mcgminus2mu(2, NDR, NLO),
        mcC(2, NDR, NNLO, NLO_QED22),  // blocks have to have the same orders
        mcP(4, NDR, NNLO, NLO_QED22),
        mcM(2, NDR, NNLO, NLO_QED22),
        mcL(2, NDR, NNLO, NLO_QED22),
        mcQ(4, NDR, NNLO, NLO_QED22),
        mcB(1, NDR, NNLO, NLO_QED22),
        Vckm(3,3,0.)
{    
    swa = 0.;
    swb = 0.;
    swc = 0.;
    xcachea = 0.;
    xcacheb = 0.;
    xcachec = 0.;
    

    for (int j=0; j<10; j++) {
        CWD1ArrayLO[j] = 0.; 
        CWD1ArrayNLO[j] = 0.;
        CWbnlepArrayLOqcd[j] = 0.;
        CWbnlepArrayNLOqcd[j] = 0.;
        CWbnlepArrayLOew[j] = 0.;
        CWbnlepArrayNLOew[j] = 0.;
    };
    
    
    for(int j=0; j<19; j++){
        CWBMllArrayLO[j] = 0.;
        CWBMllArrayNLO[j] = 0.;
    }
    
    for(int j=0; j<8; j++){
        CWbsgArrayLO[j] = 0.;
        CWbsgArrayNLO[j] = 0.;
        CWbsgArrayNNLO[j] = 0.;
        CWprimebsgArrayLO[j] = 0.;
        CWprimebsgArrayNLO[j] = 0.;
    }
    
    for(int j=0; j<8; j++){
        CWBsmmArrayNNLOqcd[j] = 0.; 
        CWBsmmArrayNLOqcd[j] = 0.;
        CWBsmmArrayLOqcd[j] = 0.;
        CWBsmmArrayNLOewt4[j] = 0.;
        CWBsmmArrayNLOewt2[j] = 0.;
        CWBsmmArrayNLOew[j] = 0.;
    }
    
    for(int j=0; j<8; j++){
        CWBdmmArrayNNLOqcd[j] = 0.; 
        CWBdmmArrayNLOqcd[j] = 0.;
        CWBdmmArrayLOqcd[j] = 0.;
        CWBdmmArrayNLOewt4[j] = 0.;
        CWBdmmArrayNLOewt2[j] = 0.;
        CWBdmmArrayNLOew[j] = 0.;
    }
    
    
    Nc = SM.getNc();
    CF = SM.getCF();
    gamma0 = 6. * (Nc - 1.) / Nc;
    J5 = SM.Beta1(5) * gamma0 / 2. / SM.Beta0(5) / SM.Beta0(5) - ((Nc - 1.)/(2. * Nc) * (-21. + 57./Nc - 19./3. * Nc + 20./3.)) / 2. / SM.Beta0(5);
    BtNDR = 5. * (Nc - 1.) / 2. / Nc + 3. * CF;
}

StandardModelMatching::~StandardModelMatching()
{}

void StandardModelMatching::updateSMParameters()
{
    Mut = SM.getMut();
    Muw = SM.getMuw(); 
    Ale = SM.getAle();
    Mt_muw = SM.Mrun(Muw, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    Mt_mut = SM.Mrun(Mut, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    alstilde = SM.Als(Muw, FULLNNNLO, true) / 4. / M_PI; //  SM.Alstilde5(Muw) WHICH ONE TO USE?
    aletilde = SM.Ale(Muw, FULLNLO) / 4. / M_PI; // Ale / 4. / M_PI; // WHERE IS ale(mu)?
    GF = SM.getGF();
    Mw_tree = SM.Mw_tree();
    Mw = SM.Mw(); /* on-shell Mw */
    sW2 = SM.sW2(); /* on-shell sW2 */
    /* NP models should be added here after writing codes for Mw. */
    /* The following is commented out till further review.*/
//    if (SM.getModelName()=="StandardModel") {
//        Mw = SM.Mw(); /* on-shell Mw */
//        sW2 = SM.sW2(); /* on-shell sW2 */
//    } else {
//        Mw = Mw_tree;
//        sW2 = 1.0 - Mw*Mw/SM.getMz()/SM.getMz();
//    }
//    sW2 = (M_PI * Ale ) / ( sqrt(2.) * GF * Mw * Mw ); // WARNING: only for checking
    Vckm = SM.getVCKM();
    lam_t = SM.computelamt();
    L = 2*log(Muw/Mw);
    Lz = 2*log(Muw/SM.getMz());
    mu_b = SM.getMub();
}

double StandardModelMatching::x_c(const double mu, const orders order) const 
{
    double mc = SM.Mrun(mu, SM.getQuarks(QCD::CHARM).getMass_scale(), 
                        SM.getQuarks(QCD::CHARM).getMass(), order);
    return mc*mc/Mw/Mw;    
}

double StandardModelMatching::x_t(const double mu, const orders order) const 
{
    double mt;
    
    if(mu == Mut)
        mt = Mt_mut;
    else if (mu == Mw)
        mt = Mt_muw;
    else
        mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), order);

    //msbar mass here?
    //mt=163.836;
    return mt*mt/Mw/Mw;   
}

double StandardModelMatching::mt2omh2(const double mu, const orders order) const 
{
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), order);
    return (mt / SM.getMHl())*(mt / SM.getMHl());
}

double StandardModelMatching::S0(double x) const 
{
    return S0(x, x);
}

double StandardModelMatching::S0(double x, double y) const 
{ // Buras 2000 Appendix
    if (fabs(1. - y / x) < LEPS){
        return ((x * (-4. + 15. * x - 12. * x * x + x * x * x +
            6. * x * x * log(x))) / (4. * pow(-1. + x, 3.)));
    }
    else
        return (x * y * ((1./4. + 3./2. / (1. - x) - 3./4. / pow(1. - x, 2.)) *
            log(x) / (x - y) +
            (1./4. + 3./2. / (1. - y) - 3./4. / pow(1. - y, 2.)) *
            log(y) / (y - x) -
            3./4. / (1. - x) / (1. - y)));
}

double StandardModelMatching::S0p( double x) const
{
    double x2 = x * x;
    return (x * (4. - 22. * x + 15. * x2 + 2. * x2 * x + x2 * x2 - 18. * x2 * log(x)) / 4. / pow(x - 1., 4.));
}

double StandardModelMatching::S11(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm3 = (x - 1.) * (x - 1.) * (x - 1.);
    double xm4 = xm3 * (x - 1);
    
    return (x * (4. - 39. * x + 168. * x2 + 11. * x3) / 4. / xm3
            + 3. * x3 * gslpp_special_functions::dilog(1. - x) * (5. + x) / xm3
            + 3. * x * log(x)*(-4. + 24. * x - 36. * x2 - 7. * x3 - x4) / 2.
            / xm4 + 3. * x3 * pow(log(x), 2.) * (13. + 4. * x + x2) / 2.
            / xm4);
}

double StandardModelMatching::S18(double x) const 
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm2 = (x - 1.) * (x - 1.);
    double xm3 = xm2 * (x - 1.);
    return ((-64. + 68. * x + 17. * x2 - 11. * x3) / 4. / xm2
            + pow(M_PI, 2.) * 8./3. / x + 2. * gslpp_special_functions::dilog(1. - x) * (8. - 24. * x
            + 20. * x2 - x3 + 7. * x4 - x4 * x) / (x * xm3)
            + log(x)*(-32. + 68. * x - 32. * x2 + 28. * x3 - 3. * x4)
            / (2. * xm3) + x2 * pow(log(x), 2.) * (4. - 7. * x + 7. * x2
            - 2. * x3) / (2. * xm2 * xm2));
}

double StandardModelMatching::S1(double x) const 
{
    return (CF * S11(x) + (Nc - 1.) / 2. / Nc * S18(x));
}

/*******************************************************************************
 * loop functions misiak base for b -> s gamma                                 * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *         
 *                 - semileptonic                                              * 
 * ****************************************************************************/

double StandardModelMatching::A0t(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    
    return ((-3. * x3 + 2. * x2)/(2. * pow(1. - x, 4.)) * log(x) + 
            (22. * x3 - 153. * x2 + 159. * x - 46.)/(36. * pow(1. - x, 3.)));
}

double StandardModelMatching::B0t(double x) const
{
    return( x / (4.* (1. - x) * (1. - x)) * log(x) + 1. / (4. * (1. - x)) );
}

double StandardModelMatching::C0t(double x) const
{
    return( (3. * x * x + 2. * x) / (8. * (1. - x) * (1. - x)) * log(x) + (-x * x + 6. * x) / (8. * (1. - x)) );
}

double StandardModelMatching::D0t(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    
    return( (-3. * x4 + 30. * x3 - 54. * x2 + 32.* x - 8.) / (18.* pow(1. - x, 4)) * log(x)
            + (-47. * x3 + 237. * x2 - 312. * x + 104.) / (108. * pow(1.- x, 3.)) );
}

double StandardModelMatching::E0t(double x) const
{
//***// CHECK THIS FUNCTION    double x2 = x * x;
    
    return (-9. * x * x + 16. * x - 4.) / (6. * pow((1.- x),4) ) * log(x) + (-7. * x * x * x - 21. * x * x + 42. * x + 4.)/(36 * pow((1. - x), 3));
    //return (x * (18. - 11. * x - x * x) / (12. * pow(1. - x, 3.) + x * x * (15. - 16. * x + 4. * x * x) /(6. * pow(1. - x, 4.)) * log(x) - 2./3. * log(x)));
}

double StandardModelMatching::F0t(double x) const
{
    double x2 = x * x;
    double xm3 = (1. - x)*(1. - x)*(1. - x);
    
    return ((3. * x2) / (2. * xm3 * (1. - x)) * log(x) + ( 5. * x2 * x - 9. * x2 + 30. * x - 8.)/
            (12. * xm3));
}

double StandardModelMatching::A1t(double x, double mu) const
{
    double x2 = x * x;
    double x3 = x * x * x;
    double x4 = x * x * x * x;
    double xm2 = pow(1. - x, 2);
    double xm3 = xm2 * (1. - x);
    double xm4 = xm3 * (1. - x);
    double xm5 = xm4 * (1. - x);
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    
    return ((32. * x4 + 244. * x3 - 160. * x2 + 16. * x)/9./xm4 * gslpp_special_functions::dilog(1.- 1./x) +
            (-774. * x4 - 2826. * x3 + 1994. *x2 - 130. * x + 8.)/81./xm5 * log(x) + 
            (-94. * x4 - 18665. * x3 + 20682. * x2 - 9113. * x + 2006.)/243./xm4 +
            ((-12. * x4 - 92. * x3 + 56. * x2)/3./(1.-x)/xm4 * log(x) + 
            (-68. * x4 - 202. * x3 - 804. * x2 + 794. * x - 152.)/27./xm4) * 2. * log(mu/mt));
}

double StandardModelMatching::B1t(double x, double mu) const
{
    double x2 = x * x;
    double xm2 = pow(1. - x, 2);
    double xm3 = pow(1. - x, 3);
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    
    return (-2. * x)/xm2 * gslpp_special_functions::dilog(1.- 1./x) + (-x2 + 17. * x)/(3. * xm3)*log(x) + (13. * x + 3.)/(3. * xm2) + 
            ( (2. * x2 + 2. * x)/xm3*log(x) + (4. * x)/xm2 ) * 2. * log(mu / mt);
    /*(2. * x)/xm2 * gsl_sf_dilog(1.-x) + (3. * x2 + x)/xm3 * log(x) * log(x) + (-11. * x2 - 5. * x)/(3. * xm3) * log(x) +
            (-3. * x2 + 19. * x)/(3. * xm2) + 16. * x * (2. * (x - 1.) - (1. + x) * log(x))/(4. * xm3) * log(mu / Mw);*/
}

double StandardModelMatching::C1t(double x, double mu) const
{
    double x2 = x * x;
    double x3 = x * x2;
    double xm2 = pow(1. - x, 2);
    double xm3 = pow(1. - x, 3);
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
   
    return (-x3 - 4. * x)/xm2 * gslpp_special_functions::dilog(1.- 1./x) + (3. * x3 + 14. * x2 + 23 * x)/(3. * xm3)*log(x) + 
            (4. * x3 + 7. * x2 + 29. * x)/(3. * xm2) + ( (8. * x2 + 2. * x)/xm3*log(x) + (x3 + x2 + 8. * x)/xm2) * 2. * log(mu / mt);
    /*(x3 + 4. * x)/xm2 * gsl_sf_dilog(1.-x) + (x4 - x3  + 20. * x2)/(2. * xm3) * log(x) * log(x) +
            (-3. * x4 - 3. * x3 - 35. * x2 + x)/(3. * xm3) * log(x) + (4. * x3 + 7. * x2 + 29. * x)/(3. * xm2) +
            16. * x * (-8. + 7. * x + x3 - 2. * (1. + 4. * x) * log(x))/(8. * xm3) * log(mu / Mw);*/
}

double StandardModelMatching::D1t(double x, double mu) const
{
    double x2 = x * x;
    double x3 = x * x2;
    double x4 = x * x3;
    double xm4 = pow(1. - x, 4);
    double xm5 = pow(1. - x, 5);
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    
    return (380. * x4 - 1352. * x3 + 1656. * x2 - 784. * x + 256.)/(81. * xm4) * gslpp_special_functions::dilog(1.- 1./x) +
            (304. * x4 + 1716. * x3 - 4644. * x2 + 2768. * x - 720.)/(81. * xm5)*log(x) +
            (-6175. * x4 + 41608. * x3 - 66723. * x2 + 33106. * x - 7000.)/(729. * xm4) +
            ( (648. * x4 - 720. * x3 - 232. * x2 - 160. * x + 32.)/(81. * xm5)*log(x) + 
            (-352. * x4 + 4912. * x3 - 8280. * x2 + 3304. * x - 880.)/(243. * xm4) ) * 2. * log(mu / mt);
}

double StandardModelMatching::F1t(double x, double mu) const
{
    double x2 = x * x;
    double x3 = x * x2;
    double x4 = x * x3;
    double xm4 = pow(1. - x, 4);
    double xm5 = pow(1. - x, 5);
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
    
    return ((4. * x4 - 40. * x3 - 41. * x2 - x)/3./xm4 * gslpp_special_functions::dilog(1.- 1./x) +
            (-144. * x4 + 3177. * x3 + 3661. * x2 + 250. * x - 32.)/108./xm5 * log(x)
            + (-247. * x4 + 11890. * x3 + 31779. * x2 - 2966. * x + 1016.)/648./xm4
            + ((17. * x3 + 31. * x2)/xm5 * log(x) + (- 35. * x4 + 170. * x3 + 447. * x2  
            + 338. * x - 56.)/18./xm4)* 2. * log(mu/mt));
}

double StandardModelMatching::E1t(double x, double mu) const

{
double x2 = x * x;
double x3 = x * x2;
double x4 = x * x3;
double xm4 = pow(1. - x, 4);
double xm5 = pow(1. - x, 5);
double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);

 
return (515. * x4 - 614. * x3 - 81. * x2 - 190. * x + 40.)/(54. * xm4) * gslpp_special_functions::dilog(1.- 1./x) +
(-1030. * x4 + 435. * x3 + 1373. * x2 + 1950. * x - 424.)/(108. * xm5) * log(x) +
(-29467. * x4 + 45604. * x3 - 30237. * x2 + 66532. * x - 10960.)/(1944. * xm4) +
( (-1125. * x3 + 1685. * x2 + 380. * x - 76.)/(54. * xm5)*log(x) + 
(133. * x4 - 2758. * x3 - 2061. * x2 + 11522. * x - 1652.)/(324. * xm4) ) * 2. * log(mu / mt);
 }

double StandardModelMatching::G1t(double x, double mu) const

{
double x2 = x * x;
double x3 = x * x2;
double x4 = x * x3;
double xm3 = pow(1. - x, 3);
double xm4 = pow(1. - x, 4);

double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
  
return (10. * x4 - 100. * x3 + 30. * x2 + 160. * x - 40.)/(27. * xm4) * gslpp_special_functions::dilog(1.- 1./x) +
(30. * x3 - 42. * x2 - 332. * x + 68.)/(81. * xm4)*log(x) +
(-6. * x3 - 293. * x2 + 161. * x + 42.)/(81. * xm3) +
( (90. * x2 - 160. * x + 40.)/(27. * xm4)*log(x) + 
(35. * x3 + 105. * x2 - 210. * x - 20.)/(81. * xm3) ) * 2. * log(mu / mt);
 }

double StandardModelMatching::C7c_3L_at_mW(double x) const

{
   double z = 1./x;
   return (1.525 - 0.1165*z + 0.01975*z*log(z) + 0.06283*z*z + 0.005349*z*z*log(z)+ 0.01005*z*z*log(z)*log(z) 
           - 0.04202*z*z*z + 0.01535*z*z*z*log(z) - 0.00329*z*z*z*log(z)*log(z) + 0.002372*z*z*z*z - 0.0007910*z*z*z*z*log(z)); 
}

double StandardModelMatching::C7t_3L_at_mt(double x) const

{
   double z = 1./x;
   return (12.06 + 12.93*z + 3.013*z*log(z) + 96.71*z*z + 52.73*z*z*log(z) 
           + 147.9*z*z*z +187.7*z*z*z*log(z) - 144.9*z*z*z*z + 236.1*z*z*z*z*log(z)); 
}

double StandardModelMatching::C7t_3L_func(double x, double mu) const

{
double x2 = x * x;
double x3 = x * x2;
double x4 = x * x3;
double x5 = x * x3;
double xm1to5 = (x-1.)*(x-1.)*(x-1.)*(x-1.)*(x-1.);
double xm1to6 = xm1to5*(x-1.);

double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
  
return ( 2. * log(mu/mt) * (gslpp_special_functions::dilog(1.- 1./x) * (-592. * x5 - 22.* x4 + 12814. * x3 - 6376. * x2 + 512. * x )/27./xm1to5 
        + log(x) * (-26838. * x5 + 25938. * x4 + 627367. * x3 - 331956. * x2 + 16989. * x - 460.)/729./xm1to6 
        + (34400. * x5 + 276644.*x4 - 2668324. * x3 + 1694437.*x2 - 323354.*x + 53077.)/2187./xm1to5) 
        + 4.*log(mu/mt)*log(mu/mt) * (log(x)*(-63. * x5 + 532. * x4 + 2089. * x3 - 1118. * x2)/9./xm1to6  
        + (1186.*x5 - 2705.*x4 - 24791.*x3 - 16099.*x2 + 19229.*x - 2740.)/162./xm1to5) );    
    
}

double StandardModelMatching::C8c_3L_at_mW(double x) const

{
   double z = 1./x;
   return (- 1.870 + 0.1010*z - 0.1218*z*log(z) + 0.1045*z*z - 0.03748*z*z*log(z) 
           + 0.01151*z*z*log(z)*log(z) - 0.01023*z*z*z + 0.004342*z*z*z*log(z) 
           + 0.0003031*z*z*z*log(z)*log(z) - 0.001537*z*z*z*z + 0.0007532*z*z*z*z*log(z));    
}

double StandardModelMatching::C8t_3L_at_mt(double x) const

{
   double z = 1./x;
   return (- 0.8954 - 7.043*z - 98.34*z*z - 46.21*z*z*log(z) - 127.1*z*z*z 
           - 181.6*z*z*z*log(z) + 535.8*z*z*z*z - 76.76*z*z*z*z*log(z));    
}

double StandardModelMatching::C8t_3L_func(double x, double mu) const

{
double x2 = x * x;
double x3 = x * x2;
double x4 = x * x3;
double x5 = x * x3;
double xm1to5 = (x-1.)*(x-1.)*(x-1.)*(x-1.)*(x-1.);
double xm1to6 = xm1to5*(x-1.);

double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
  
return ( 2. * log(mu/mt) * (gslpp_special_functions::dilog(1.- 1./x) * (-148. * x5 + 1052. * x4 - 4811. * x3 - 3520. * x2 - 61. * x)/18./xm1to5 
        + log(x) * (-15984. * x5 + 152379. * x4 - 1358060. * x3 - 1201653. * x2 - 74190. * x + 9188.)/1944./xm1to6 
        + (109669. * x5 - 1112675. * x4 + 6239377. * x3 + 8967623. * x2 + 768722. * x - 42796.)/11664./xm1to5) 
        + 4. * log(mu/mt) * log(mu/mt) * (log(x) * (-139. * x4 - 2938. * x3 - 2683. * x2)/12./xm1to6 
        + (1295. * x5 - 7009. * x4 + 29495. * x3 + 64513. * x2 + 17458. * x - 2072.)/216./xm1to5) );    
    
}

double StandardModelMatching::Tt(double x) const
{
return ((-(16. * x + 8.) * sqrt(4. * x - 1.) * gslpp_special_functions::clausen(2. * asin(1./2./sqrt(x)))) +((16. * x + 20./3.) * log(x)) + (32. * x) + (112./9.)) ;
}

double StandardModelMatching::Wt(double x) const
{
 double x2 = x * x;
 double x3 = x * x * x;
 double x4 = x * x * x * x;
 double xm2 = pow(1. - x, 2);
 double xm3 = xm2 * (1. - x);
 double xm4 = xm3 * (1. - x);

return ((-32. * x4 + 38. * x3 + 15. * x2 - 18. * x)/18./xm4 * log(x) -
 (-18. * x4 + 163. * x3 - 259. *x2 + 108. * x)/36./xm3 );
}

double StandardModelMatching::Eet(double x) const
{
 double x2 = x * x;
 double xm2 = pow(1. - x, 2);
 double xm3 = xm2 * (1. - x);
 double xm4 = xm3 * (1. - x);

return ((x * (18. - 11. * x - x2))/(12. * xm3) +
 (log(x) * (x2 * (15. - 16. * x + 4. * x2))/(6. * xm4)) - 2. * log(x) /3.);
}

double StandardModelMatching::Rest(double x, double mu) const

{
   double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);

return (-37.01364013973161 + 7.870950908767437 * mt - 
    0.0015295355176062242 * mt * mt + 2.41071411865951 * Mw - 
    0.20348320015672194 * mt * Mw - 0.02858827491583899 * Mw * Mw + 
    0.001422822903303167 * mt * Mw * Mw + (4.257050684362808 + 0.17719711396626878 * mt - 
       0.8190947921716011 * Mw - 0.002315407459561656 * mt * Mw + 0.008797408866807221 * Mw * Mw) * log(
      mu) + (0.49627858125619595 - pow(5.784745743815408,-8) *mt + 
       0.031869225004473686 * Mw - 0.00041193393986696286 * Mw * Mw) * log(
       mu) * log(mu));
 }

double StandardModelMatching::Y0(double x) const
{
    return( x/8. * ((4 - 5 * x + x * x + 3 * x * log(x))/pow(x - 1., 2.)) );
}

double StandardModelMatching::Y1(double x, double mu) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double logx = log(x);
    double xm = 1. - x;
    double xm2 = xm * xm;
    double xm3 = xm2 * xm;
    
    return ((10. * x + 10. * x2 + 4. * x3)/(3 * xm2) - (2. * x - 8. * x2 - x3 - x4)/xm3 * logx 
            + (2. * x - 14. * x2 + x3 - x4)/(2. * xm3) * pow(logx, 2.) 
            + (2. * x + x3)/xm2 * gslpp_special_functions::dilog(1. - x) 
            + 16. * x * (-4. + 3. * x + x3 - 6. * x * logx)/(8. * -xm3) * log(mu / Mw));
}

double StandardModelMatching::C7LOeff(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    
    return( (3. * x3 - 2. * x2) / (4. * pow(x - 1., 4.)) * log(x) + (-8. * x3 - 5. * x2 +
              7. * x) / (24. * pow(x-1.,3.)));
}

double StandardModelMatching::C8LOeff(double x) const
{
    return( -3. * x * x / (4. * pow( x - 1., 4.)) * log(x) + (-x * x * x + 5. * x * x + 2. * x) / (8. * pow(x - 1., 3)) );
}

double StandardModelMatching::C7NLOeff(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm4 = pow(x-1.,4.);
    double xm5 = xm4 * (x - 1);
    double logx = log(x);
    
    double Li2 = gslpp_special_functions::dilog(1.-1./x);
    return( Li2 * ( -16. * x4 - 122. * x3 + 80. * x2 - 8. * x) / (9. * xm4) +
            (6. * x4 + 46. * x3 -28. * x2) / (3. * xm5) * logx * logx +
            (-102. * x4 * x - 588. * x4 - 2262. * x3 + 3244. * x2 - 1364. * x + 208.) / (81. * xm5) * logx +
            (1646. * x4 + 12205. * x3 - 10740. * x2 + 2509. * x - 436.) / (486. * xm4));
}

double StandardModelMatching::C8NLOeff(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm4 = pow(x-1.,4.);
    double xm5 = xm4 * (x - 1);
    double logx = log(x);
    
    double Li2 = gslpp_special_functions::dilog(1.-1./x);
    return(Li2 * ( -4. * x4 + 40. * x3 + 41. * x2 + x) / (6. * xm4) +
            (-17. * x3 - 31. * x2) / (2. * xm5) * logx * logx +
            (-210. * x * x4 + 1086. * x4 + 4893. * x3 + 2857. * x2 - 1994. * x + 280.)/(216. * xm5) * logx +
            (737. * x4 - 14102. * x3 - 28209. * x2 + 610. * x - 508.) / (1296. * xm4));
}


/******************************************************************************/


/*******************************************************************************
 * loop functions Buras base for nonlep. b decays                              * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - ew penguins                                               *
 * ****************************************************************************/

double StandardModelMatching::B0b(double x) const
{
    return ( 0.25 * ( x / (1. - x) + x / (x * x - 2. * x + 1.) * log(x) ) );
}

double StandardModelMatching::C0b(double x) const
{
    return ( x / 8. * ( (x - 6.) / (x - 1.) + (3. * x + 2.) / ( x * x - 2. * x + 1.) * log(x) ) );
}

double StandardModelMatching::D0b(double x) const
{
    double x2 = x * x;
    return ( -4. / 9. * log(x) + (-19. * x2 * x + 25. * x2) / (36. * (x2 * x - 3. * x2 + 3. * x - 1.))
            + (x2 * (5. * x2 - 2. * x - 6.) ) / (18. * pow(x - 1.,4.)) * log(x) );
}

double StandardModelMatching::D0b_tilde(double x) const
{
    return (D0b(x) - 4./9.);
}

double StandardModelMatching::E0b(double x) const
{
    double x2 = x * x;
    
    return ( -2./3. * log(x) + (18. * x - 11. * x2 - x2 * x) / (12.* (- x2 * x +3. * x2 -3. * x +1.)) +
            (x2 * (15. - 16. * x +4. * x2))/(6.*pow(1.-x,4.)) *log(x) );
}

//Loop functions for QED corrections
double StandardModelMatching::B1d(double x, double mu) const
{
    double xmo = x - 1.;
    double dilog1mx = gslpp_special_functions::dilog(1.-x);
    double xmuw = mu*mu/Mw/Mw;
    double mut = SM.getQuarks(QCD::TOP).getMass_scale();
    double xmut = mut*mut/Mw/Mw;

    return (-(8.-183.*x+47.*x*x)/24./xmo/xmo - (8.+27.*x+93.*x*x)/24./xmo/xmo/xmo*log(x) +
            (27.*x+71.*x*x-2.*x*x*x)/24./xmo/xmo/xmo*log(x)*log(x) - (2.-3.*x-9.*x*x+x*x*x)/6./x/xmo/xmo*dilog1mx +
            (2.+x)/36./x*M_PI*M_PI + 19./6.*B0b(x) - B0b(x)*log(xmuw) + 4.*x/xmo/xmo*log(xmut) -
            (2.*x+2.*x*x)/xmo/xmo/xmo*log(x)*log(xmut));
}

double StandardModelMatching::B1d_tilde(double x, double mu) const
{
    double xmo = x - 1.;
    double dilog1mx = gslpp_special_functions::dilog(1.-x);
    double xmuw = mu*mu/Mw/Mw;
    
    return ((-8.-23.*x)/8./xmo - (8.-5.*x)/8./xmo/xmo*log(x) + (3.*x+2.*x*x)/8./xmo/xmo*log(x)*log(x) +
            (2.-3.*x+3.*x*x+x*x*x)/2./x/xmo/xmo*dilog1mx - (2.+x)/12./x*M_PI*M_PI + 5./2.*B0b(x) +
            3.*B0b(x)*log(xmuw));
}

double StandardModelMatching::B1u(double x, double mu) const
{
    double xmo = x - 1.;
    double dilog1mx = gslpp_special_functions::dilog(1.-x);
    double xmuw = mu*mu/Mw/Mw;
    double mut = SM.getQuarks(QCD::TOP).getMass_scale();
    double xmut = mut*mut/Mw/Mw;

    return (-(46.*x+18.*x*x)/3./xmo/xmo - (16.*x-80.*x*x)/3./xmo/xmo/xmo*log(x) -
            (9.*x+23.*x*x)/2./xmo/xmo/xmo*log(x)*log(x) - 6.*x/xmo/xmo*dilog1mx - 38./3.*B0b(x) +
            4.*B0b(x)*log(xmuw) - 16.*x/xmo/xmo*log(xmut) + (8.*x+8.*x*x)/xmo/xmo/xmo*log(x)*log(xmut));
}

double StandardModelMatching::B1u_tilde(double x, double mu) const
{
    double xmo = x - 1.;
    double dilog1mx = gslpp_special_functions::dilog(1.-x);
    double xmuw = mu*mu/Mw/Mw;
    
    return ((-8.-23.*x)/8./xmo - (8.-5.*x)/8./xmo/xmo*log(x) + (3.*x+2.*x*x)/8./xmo/xmo*log(x)*log(x) +
            (2.-3.*x+3.*x*x+x*x*x)/2./x/xmo/xmo*dilog1mx - (2.+x)/12./x*M_PI*M_PI + 5./2.*B0b(x) +
            3.*B0b(x)*log(xmuw));
}

double StandardModelMatching::C1ew(double x) const
{
    double xmo = x - 1.;
    double dilog1mx = gslpp_special_functions::dilog(1.-x);
    double mut = SM.getQuarks(QCD::TOP).getMass_scale();
    double xmut = mut*mut/Mw/Mw;
    
    return ((29.*x+7.*x*x+4.*x*x*x)/3./xmo/xmo + (x-35.*x*x-3.*x*x*x-3.*x*x*x*x)/3./xmo/xmo/xmo*log(x) +
            (20.*x*x-x*x*x+x*x*x*x)/2./xmo/xmo/xmo*log(x)*log(x) + (4.*x+x*x*x)/xmo/xmo*dilog1mx +
            (8.*x+x*x+x*x*x)/xmo/xmo*log(xmut) + (2.*x+8.*x*x)/xmo/xmo/xmo*log(x)*log(xmut));
}

double StandardModelMatching::Zew(double xt, double xz) const
{
    double z0ew, z1ew;
#ifdef ZEW_NUMERIC
    z0ew = 5.1795 + 0.038*(Mt_muw-166.) + 0.015*(Mw-80.394);
    z1ew = 2.1095 + 0.0067*(Mt_muw-166.) + 0.026*(Mw-80.394);
#else    
    double xt2 = xt*xt;
    double xt3 = xt2*xt;
    double xz2 = xz*xz;
    double xz3 = xz2*xz;
    double xtmo = xt - 1.;
    double xzmo = xz - 1.;
    double dilog1mxt = gslpp_special_functions::dilog(1.-xt);
    double dilog1mxz = gslpp_special_functions::dilog(1.-xz);
    z0ew = -xt*(20.-20.*xt2-457.*xz+19.*xt*xz+8.*xz2)/32./xtmo/xz +
            xt*(10.*xt3-11.*xt2*xz-xt*(30.-16.*xz)+4.*(5.-17.*xz+xz2))/16./xtmo/xtmo/xz*log(xt) +
            xt*(10.-10.*xt2-17.*xz-xt*xz-4.*xz2)/16./xtmo/xz*log(xz) -
            xz*(10.*xt2-xt*(4.-xz)+8.*xz)/32./xtmo/xtmo*log(xt)*log(xt) - xz2/4.*log(xz)*log(xz) -
            ((8.+12.*xt+xt2)/4./xz - 5.*xtmo*xtmo*(2.+xt)/16./xz2 -
            (12.-3.*xt3-3.*xt2*(4.-xz)+4.*xt*(3.-xz)+4.*xz-xz2)/8./xtmo/xtmo)*log(xt)*log(xz) -
            ((8.+12.*xt+xt2)/2./xz - 5.*xtmo*xtmo*(2.+xt)/8./xz2 - 3.*(4.+8.*xt+2.*xt2-xt3)/4./xtmo/xtmo)*dilog1mxt +
            xzmo*xzmo*(5.-6.*xz-5.*xz2)/4./xz2*dilog1mxz - (5.-16.*xz+12.*xz2+2*xz3*xz)/24./xz2*M_PI*M_PI +
            xt*(4.-xz)(88.-30.*xz-25.*xz2-2.*xt*(44.-5.*xz-6.*xz2))/32./xtmo/xtmo/xz*phi_z(xz/4.) +
            (16.*xt3*xt-xt*(20.-xz)*xz2+8.*xz3-8.*xt3*(14.+5.*xz)+8.*xt2*(12.-7.*xz+xz2))/32./xtmo/xtmo/xz*phi_z(xz/4./xt) -
            ((22.+33.*xt-xt2)/16./xtmo/xz - 5.*xtmo*(2.+xt)/16./xz2 +
            (2.+5.*xt2+10.*xz+xt*(15.+xz))/16./xtmo/xtmo)*phi_xy(xt,xz);

    z1ew = xt*(20.-20.*xt2-265.*xz+67.*xt*xz+8.*xz2)/48./xtmo/xz -
            xt*(10.*xt3-15.*xt2*xz+4.*(5.-7.*xz+2.*xz2)-xt*(30.+20.*xz+4.*xz2))/24./xtmo/xtmo/xz*log(xt) -
            xt*(10.-10.*xt2-33.*xz+15.*xt*xz-4.*xz2)/24./xtmo/xz*log(xz) +
            xz*(8.-16.*xt+2.*xt2+10.*xz+7.*xt*xz)/48./xtmo/xtmo*log(xt)*log(xt) + xz*(4.+5.*xz)/24.*log(xz)*log(xz) +
            ((20.+6.*xt+xt2)/12./xz - 5.*xtmo*xtmo*(2.+xt)/24./xz2 +
            (3.*xt3+2.*xt2*(12.-xz)-xt*(18.-16.*xz+xz2)-2.*(9.+4.*xz-xz2))/12./xtmo/xtmo)*log(xt)*log(xz) +
            ((20.+6.*xt+xt2)/6./xz - 5.*xtmo*xtmo*(2.+xt)/12./xz2 - (6.+6.*xt-8.*xt2-xt3)/2./xtmo/xtmo)*dilog1mxt -
            xzmo*xzmo*(5.-10.*xz-7.*xz2)/6./xz2*dilog1mxz + (10.-40.*xz+36.*xz2+4.*xz3+5.*xz3*xz)/72./xz2*M_PI*M_PI +
            xt*(xz-4.)*(24.-26.*xz-13.*xz2-6.*xt*(4.-xz-xz2))/16./xtmo/xtmo/xz*phi_z(xz/4.) - (2.*xt2*(2.+xt)/3./xtmo/xz -
            (24.*xt3+12.*xt2*(14.+xz)-2.*xz*(4.+5.*xz)-xt*(80.-36.*xz+7.*xz2))/48./xtmo/xtmo)*phi_z(xz/4./xt) +
            ((10.-xt-xt2)/8./xtmo/xz - 5*xtmo*(2.+xt)/24./xz2 + (6.+3.*xt2+14.*xz+5.*xt*(7.-xz))/24./xtmo/xtmo)*phi_xy(xt,xz);
#endif
    return(z0ew+sW2*z1ew);
}


double StandardModelMatching::Gew(double xt, double xz, double mu) const
{
    double xmuw = mu*mu/Mw/Mw;

    return (Zew(xt,xz) + 5.*C0b(xt) + 6.*C0b(xt)*log(xmuw));
}

double StandardModelMatching::Hew(double xt, double xz, double mu) const
{
    double xmuw = mu*mu/Mw/Mw;
    
    return (Zew(xt,xz) - 7.*C0b(xt) + 6.*C0b(xt)*log(xmuw));
}

/******************************************************************************/
/* loop functions for rare K and B decays, K-> pi nu nu & B-> Xs nu nu        */
/******************************************************************************/

double StandardModelMatching::X0t(double x) const{
    return((x/8.)*((x+2.)/(x-1.)+(3.*x -6)/(x-1.)/(x-1.)*log(x)));
}

double StandardModelMatching::X1t(double x) const{
    
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm3 = pow(1.-x,3.);
    double logx = log(x);
    
    return (-(29. * x - x2 -4. * x3) / (3. * (1. - x) * (1. - x))
            - logx * (x + 9. * x2 - x3 - x4) / xm3
            + logx * logx * (8. * x + 4. * x2 + x3 - x4) / (2. * xm3)
            - gslpp_special_functions::dilog(1.-x) * (4. * x - x3) / ((1. - x) * (1. - x))
            - 8. * x * log(Mut*Mut/Muw/Muw) * (8. - 9. * x + x3 + 6. * logx)/8./xm3 );
  }

double StandardModelMatching::Xewt(double x, double a, double mu) const{
    double b = 0.;
    // WARNING: check consistency of EW scheme choice (see Gorbahn's NNLO papers)
    double swsq = (M_PI * Ale )/( sqrt(2) * GF * Mw * Mw);
    
    double A[17], C[17];
    
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double x5 = x4 * x;
    double x6 = x5 * x;
    double x7 = x6 * x;
    double x8 = x7 * x;
    double M_PI2 = M_PI * M_PI;
    double a2 = a * a;
    double a3 = a2 * a;
    double a4 = a3 * a;
    double a5 = a4 * a;
    double a6 = a5 * a;
    double xm2 = (x - 1.) * (x - 1.);
    double xm3 = xm2 * (x - 1.);
    double axm = a * x - 1.;
    double logx = log(x);
    double loga = log(a);
    
    A[0] = (16. - 48. * a) * M_PI2 + (288. * a - (32. - 88. * a) * M_PI2 ) * x
        + (2003. * a + 4. * (4. - 6. * a - a2 ) * M_PI2 )* x2
        + (9. * a * (93. + 28. * a) - 4. * a * (3. - 2. * a + 8. * a2) * M_PI2 ) * x3
        + (3. * a * (172. - 49. * a - 32. * a2) + 4. * a * (20. - a + 16. * a2 ) * M_PI2 ) * x4
        - (3. * a * (168. + 11. * a - 24. * a2 ) + 4. * a * (45. + 8. *a2) * M_PI2) * x5
        + 96. * a * M_PI2 * x6;
    
    A[1] = -768. * x - (525. - 867. * a)  * x2 + (303. + 318. * a) * x3 - 195. * a * x4;
    
    A[2] = -8.*(95. - 67. *  a + 11. * a2 ) * x2 + 2. * (662. - 78. *  a - 177. * a2 + 40. * a3 ) * x3
        - (608. + 476. *  a - 595. * a2 + 114. * a3 ) * x4 
        + (44. + 188. *  a - 321. * a2 + 103. * a3 - 8. * a4 ) * x5
        - a*(28. - 72. *  a + 33. * a2 - 4. * a3 ) * x6;
    
    A[3] = +48. - 10.*(57. + 4. * a) * x+ 51.*(29. + 10. * a) * x2 -
        (841. + 1265. * a) * x3 + (308. + 347. * a) * x4
        - (28. - 40. * a) * x5 + 12. * a * x6 ;
    
    A[4] = + 768. + (816. - 768. * a) * x+ (1240. - 1232. * a) * x2
        - 4.*(415. + 2.  * a) * x3 + (311. + 722.  * a) * x4
        + (145. - 267. * a) * x5 - (36. + 51. * a) * x6 + 20. * a * x7 ;
    
    A[5] = + 328. * x- (536. + 900. * a) * x2 + (208. + 1584. * a + 670. * a2 ) * x3
        - a * (668. + 1161. *  a + 225. * a2 ) * x4
        + a2 * (479. + 362. *  a + 28. * a2 ) * x5
        - a3 *(143. + 42. * a) * x6 + 16. * a4 * x7;
    
    A[6] = + 32. - 4.*(44. - 9. * a) * x + (384. - 322. *  a - 400. * a2 ) * x2
        - (400. - 869. *  a - 1126. * a2 - 696. * a3 ) * x3
        + 2.*(80. - 488. *  a - 517. * a2 - 631. * a3 - 264. * a4 ) * x4
        + (48. + 394. *  a + 269. * a2 + 190. * a3 + 882. * a4 + 196. * a5 ) * x5
        - (64. - 58. *  a - 89. * a2 - 95. * a3 + 34. * a4 + 296. * a5 + 32. * a6 ) * x6
        + (16. - 59. *  a - 79. * a2 + 256. * a3 - 239. * a4 
        + 57. * a5 + 48. * a6 ) * x7
        + (1. - a) * (1. - a) * (1. - a) * a2 * (29. + 16. * a) * x8 ;
    
    A[7] = + 28. * a2 * x2 - 32. * a3 * x3;
    
    A[8] = - 288. + 36.*(1. + 8. * a) * x + 6.*(647. + 87. * a) * x2 + 5.*(55. - 927. *  a - 132. * a2 ) * x3
        - (1233. + 98. *  a - 879. * a2 - 192. * a3 ) * x4 
        + (360. + 1371. *  a - 315. * a2 - 264. * a3 ) * x5
        - 24. * a * (17. - 4. * a2) * x6;
    
    A[9] = + 32. + 4.*(-44. + 29. * a) * x- 12.*(-32. + 77. *  a + 31. * a2 ) * x2
        + 2.*(-200. + 837. *  a + 767. * a2 + 182. * a3 ) * x3 
        - 2.*(-80. + 625. *  a + 905. * a2 + 520. * a3 + 82. * a4 ) * x4
        + (48. + 1079. *  a + 590. * a2 + 1002. * a3 + 462. * a4 + 32. * a5 ) * x5
        + (-64. - 1160. *  a - 501. * a2 - 364. * a3 - 486. * a4 - 72. * a5 ) * x6
        + (16. + 729. *  a + 1038. * a2 + 38. * a3 + 238. * a4 + 52. * a5 ) * x7
        - a*(192. + 743. *  a + 50. * a3 + 12. * a4 ) * x8 + 192. * a2 * x8 * x;
    
    A[10] = + 16. * x + 324. * x2 - 36. * x4;
    
    A[11] = + 216. * x - 672. * x2 + 152. * x3;
    
    A[12] = - 16. * x + (16. - 42. * a) * x2 + (16. + 21. *  a + 60. * a2 ) * x3
        - (16 - 21. *  a + 45. * a2 + 32. * a3 ) * x4 - a2 * (7. - 24. * a) * x5;
    
    A[13] = - 32. + (144. - 68. * a) * x + (-240. + 334. *  a + 332. * a2 ) * x2
        + (160. - 551. *  a - 660. * a2 - 364. * a3 ) * x3
        + a * (329. + 451. *  a + 650. * a2 + 164. * a3 ) * x4 
        + (-48. - a - 59. * a2 - 523. * a3 - 316. * a4 - 32. * a5 ) * x5
        + (16. - 43. *  a -93. * a2 + 255. * a3 + 287. * a4 + 32. * a5 ) * x6 
        - a2 * (-29. + 42. *  a + 103. * a2 + 8. * a3 ) * x7;

    A[14] = - 144.*(1. - a)*(1. - a) * x2 + 144. * (1. - a) * (1. - a) * x3 - 36. * (1. - a) * (1. - a) * x4;
    
    A[15] = - 32. + 96. *  a + (48. - 32. * a) * x - 176. * a * x2 - (16. - 74. * a) * x3 + 212. * a * x4;
    
    A[16] = - 32. + (64. - 100. * a) * x- 8.*(4. - 34. *  a -29. * a2 ) * x2
        - 4. * a * (34. + 170. *  a + 33. * a2 ) * x3
        + 8. * a2 * (47. + 51. *  a + 4. * a2) * x4 - 16. * a3 * (15. + 4. * a) * x5 
        + 32. * a4 * x6;
    
    C[0] = 1. / (3.* a * xm2 * x);
    
    C[1] = phi1(0.25) / (xm3 * axm);
    
    C[2] = phi1(0.25 * a) / (2. * xm3 * axm);
    
    C[3] = phi1(1. / 4. / x) / (2. * xm3 * axm);
    
    C[4] = phi1(0.25 * x) / (2. * xm3 * axm);
    
    C[5] = phi1(a * x * 0.25) / (xm3 * axm);
    
    C[6] = phi2(1. / a / x, 1. / a) / (2. * a2 * x2 * xm3 * axm);
    
    C[7] = loga * log(a) / axm;
    
    C[8] = logx / (xm3 * axm * 3.);
    
    C[9] = logx * logx / ((x-1.) * xm3 * axm * 2. * a * x);
    
    C[10] = 2. * log(mu/Mw) / xm2;
    
    C[11] = logx * 2. * log(mu/Mw) / xm3;
    
    C[12] = loga / (xm2 * axm);
    
    C[13] = logx * loga / (2. * a * xm3 * x * axm);
    
    C[14] = gslpp_special_functions::dilog(1. - a) / xm2;
    
    C[15] = gslpp_special_functions::dilog(1. - x) / a / x;
    
    C[16] = gslpp_special_functions::dilog(1. - a * x) / (a * x * xm2); 
    
    for (int i=0; i<10; i++){
        b += C[i]*A[i];
    }
    
    return (b/128./swsq);
}

double StandardModelMatching::phi1(double z) const{
    if (z >= 0.) {
        if (z < 1){
                return(4. * sqrt(z / (1. - z)) * gslpp_special_functions::clausen(2. * asin(sqrt(z))));
        }
        else{
            return((1. / sqrt(1. - 1. / z)) * (2. * log(0.5 * sqrt(1. - 1. / z)) * log(0.5 * sqrt(1. -1. / z))- 4. * gslpp_special_functions::dilog(0.5 * (1. - sqrt(z / (1. - z)))) - log(4. * z) * log(4. * z)
                    + M_PI * M_PI / 3.));
        }
    }
    else{
        std::stringstream out;
        out << z;
        throw std::runtime_error("StandardModelMatching::phi1(double z)" + out.str() + " <0");
    }
    return(0.);
}

double StandardModelMatching::phi2(double x, double y) const{
    double l = sqrt((1. - x - y) * (1. - x - y) - 4. * x * y);
    
    if ((l * l) >= 0. || (sqrt(x) + sqrt(y)) <= 1.){
        return( 1. / l * (M_PI * M_PI/3. + 2. * log(0.5 * (1. + x - y - l)) * log(0.5 * (1. - x + y - l)) - log(x) * log(y) - 2. * gslpp_special_functions::dilog(0.5 * (1. + x - y - l)) - 2. * gslpp_special_functions::dilog(0.5 * (1. -x + y - l))));
    }
    else if((l * l) < 0. || (sqrt(x) + sqrt(y)) > 1.){
        return(2./( -l * l) * (gslpp_special_functions::clausen(2. * acos((-1. + x + y)/(2. * sqrt(x * y))))
                +gslpp_special_functions::clausen(2. * acos((1. + x - y) / (2. * sqrt(x))))
                +gslpp_special_functions::clausen(2. * acos((1. - x + y) / (2. * sqrt(y))))));
    }
    else{
        std::stringstream out;
        out << x;
        throw std::runtime_error("StandardModelMatching::phi2(double x, double y) wrong" + out.str());
    }
    return(0.);
}

double StandardModelMatching::phi_z(double z) const
{
    double beta = sqrt(1.-1./z);
    double clausen = gslpp_special_functions::clausen(2.*asin(sqrt(z)));
    double dilog = gslpp_special_functions::dilog((1.-beta)/2.);
    
    if (z > 0.) {
        if (z <= 1.){
                return(4.*sqrt(z/(1.-z)) * clausen);
        }
        else{
            return(1./beta*(2.*log((1.-beta)/2.)*log((1.-beta)/2.) - 4.*dilog - log(4.*z)*log(4.*z) + M_PI*M_PI/3.));
        }
    }
    else{
        std::stringstream out;
        out << z;
        throw std::runtime_error("StandardModelMatching::phi_z(double z)" + out.str() + " <0");
    }
}

double StandardModelMatching::phi_xy(double x, double y) const
{
    double lambda = sqrt((1.-x-y)*(1.-x-y) - 4.*x*y);
    double diloga = gslpp_special_functions::dilog((1.+x-y-lambda)/2.);
    double dilogb = gslpp_special_functions::dilog((1.-x+y-lambda)/2.);
    double clausenxy = gslpp_special_functions::clausen(2.*acos((-1.+x+y)/2./sqrt(x*y)));
    double clausenx = gslpp_special_functions::clausen(2.*acos((1.+x-y)/2./sqrt(x)));
    double clauseny = gslpp_special_functions::clausen(2.*acos((1.-x+y)/2./sqrt(y)));
    
    if ((lambda*lambda) >= 0.){
        return(lambda*(2.*log((1.+x-y-lambda)/2.)*log((1.-x+y-lambda)/2.) - log(x)*log(y) -
                  2.*diloga - 2.*dilogb + M_PI*M_PI/3.));
    }
    else if((lambda*lambda) < 0.){
        return(-2.*sqrt(-lambda*lambda)*(clausenxy + clausenx + clauseny));
    }
    else{
        std::stringstream out;
        out << x;
        throw std::runtime_error("StandardModelMatching::phi_xy(double x, double y) wrong" + out.str());
    }
}

/*******************************************************************************
 * Wilson coefficients Buras base for Delta B = 2 observables                  *                                           
 * ****************************************************************************/

 std::vector<WilsonCoefficient>& StandardModelMatching::CMdbd2()  
{
    double gammam = 6. * CF;
    double Bt;  
    
    
    double xt = x_t(Mut);
    gslpp::complex co = GF / 4. / M_PI * Mw * SM.computelamt_d();
    
    vmcdb.clear();

    switch (mcdbd2.getScheme()) {
        case NDR:
            Bt = BtNDR;
            break;
        case HV:
        case LRI:
        default:
            std::stringstream out;
            out << mcdbd2.getScheme();
            throw std::runtime_error("StandardModel::CMdb2(): scheme " + out.str() + "not implemented"); 
    }

    mcdbd2.setMu(Mut);
    
    switch (mcdbd2.getOrder()) {
        case NNLO:
        case NLO:
            mcdbd2.setCoeff(0, co * co * 4. * (SM.Als(Mut, FULLNLO) / 4. / M_PI * (S1(xt) + //* CHECK ORDER *//
                    (Bt + gamma0 * log(Mut / Mw)) * S0(xt, xt) + 2. * gammam * S0p(xt) * log(Mut / Mw))), NLO);
#if SUSYFIT_DEBUG & 1
    std::cout << "Mw = " << Mw << " xt(muw=" << Muw << ")= " << xt << "matching of DB=2: S0(xt) = " << S0(xt) << 
                ", S1(xt) = " << S1(xt) +
                    (Bt + gamma0 * log(Muw / Mw)) * S0(xt, xt) + 2. * gammam * S0p(xt) * log(Muw / Mw) 
            << ", lambdat_d^2 = " << SM.getlamt_d()*SM.getlamt_d() << std::endl;
#endif
        case LO:
            mcdbd2.setCoeff(0, co * co * 4. * (S0(xt, xt)), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbd2(): order " + out.str() + "not implemented"); 
    }
    

    vmcdb.push_back(mcdbd2);
    return(vmcdb);
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMdbs2() 
{
    double gammam = 6. * CF;
    double Bt;
    double xt = x_t(Mut);
    gslpp::complex co = GF / 4. / M_PI * Mw * SM.computelamt_s();

    vmcds.clear();

    switch (mcdbs2.getScheme()) {
        case NDR:
            Bt = BtNDR;
            break;
        case HV:
        case LRI:
        default:
            std::stringstream out;
            out << mcdbs2.getScheme();
            throw std::runtime_error("StandardModel::CMdbs2(): scheme " + out.str() + "not implemented"); 
    }

    mcdbs2.setMu(Mut);
 
    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:          
            mcdbs2.setCoeff(0, co * co * 4. * (SM.Als(Mut, FULLNLO) / 4. / M_PI * (S1(xt) + //* CHECK ORDER *//
                    (Bt + gamma0 * log(Mut / Mw)) * S0(xt, xt) + 2. * gammam * S0p(xt) * log(Mut / Mw))), NLO);
         case LO:
            mcdbs2.setCoeff(0, co * co * 4. * (S0(xt, xt)), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbs2(): order " + out.str() + "not implemented"); 
    }    

    vmcds.push_back(mcdbs2);
    return(vmcds);
}

/*******************************************************************************
 * Wilson coefficients Buras base for Delta S = 2 observables                  *
 * Comment: they look empty because they are computed in Flavour/EvolDF2.cpp   *
 *          due to historical reasons.                                         *
 * ****************************************************************************/

 std::vector<WilsonCoefficient>& StandardModelMatching::CMdk2() 
{
    vmck2.clear();
    
    mcdk2.setMu(Mut);
 
    switch (mcdk2.getOrder()) {
        case NNLO:
        case NLO:
            mcdk2.setCoeff(0, 0., NLO);
        case LO:
            mcdk2.setCoeff(0, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcdk2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdk2(): order " + out.str() + "not implemented"); 
    }

    vmck2.push_back(mcdk2);
    return(vmck2);
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMd1Buras() 
{    
    vmcd1Buras.clear();
    
    switch (mcd1Buras.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcd1Buras.getScheme();
            throw std::runtime_error("StandardModel::CMd1Buras(): scheme " + out.str() + "not implemented"); 
    }

    mcd1Buras.setMu(Muw);
    
    switch (mcd1Buras.getOrder()) {
        case NNLO:
        case NLO:
            mcd1Buras.setCoeff(0, SM.Als(Muw, FULLNLO) / 4. / M_PI  * 11./2. , NLO); //* CHECK ORDER *//
            mcd1Buras.setCoeff(1, SM.Als(Muw, FULLNLO) / 4. / M_PI * (-11./6.) , NLO);
            for (int j=2; j<10; j++){
            mcd1Buras.setCoeff(j, 0., NLO);    
            }
            case LO:
            mcd1Buras.setCoeff(0, 0., LO);
            mcd1Buras.setCoeff(1, 1., LO);
            for (int j=2; j<10; j++){
            mcd1Buras.setCoeff(j, 0., LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcd1Buras.getOrder();
            throw std::runtime_error("StandardModelMatching::CMd1Buras(): order " + out.str() + "not implemented"); 
    }
        
    vmcd1Buras.push_back(mcd1Buras);
    
    return(vmcd1Buras);

}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMd1() 
{ 
    vmcd1.clear();
    
    switch (mcd1.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcd1.getScheme();
            throw std::runtime_error("StandardModel::CMd1(): scheme " + out.str() + "not implemented"); 
    }

    mcd1.setMu(Muw);
    
    switch (mcd1.getOrder()) {
        case NNLO:
        case NLO:
            mcd1.setCoeff(0, SM.Als(Muw, FULLNLO) / 4. / M_PI  * 15. , NLO); //* CHECK ORDER *//
            for (int j=1; j<10; j++){
            mcd1.setCoeff(j, 0., NLO);
            }
        case LO:
            mcd1.setCoeff(0,  0. , LO);
            mcd1.setCoeff(1,  1. , LO);
            for (int j=2; j<10; j++){
            mcd1.setCoeff(j, 0., LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcd1.getOrder();
            throw std::runtime_error("StandardModelMatching::CMd1(): order " + out.str() + "not implemented"); 
    }
      
    vmcd1.push_back(mcd1);
    
    return(vmcd1);
    
}

  std::vector<WilsonCoefficient>& StandardModelMatching::CMdd2() 
{
    vmcd2.clear();

    switch (mcdd2.getScheme()) {
        case NDR:
        case HV:
        case LRI:
          break; 
        default:
            std::stringstream out;
            out << mcdd2.getScheme();
            throw std::runtime_error("StandardModel::CMdd2(): scheme " + out.str() + "not implemented"); 
    }

    mcdd2.setMu(Muw);
 
    switch (mcdd2.getOrder()) {
        case NNLO:
        case NLO:
            for(int i=0; i<5; i++)
                mcdd2.setCoeff(i, 0., NLO);
        case LO:
            for(int j=0; j<5; j++)
                mcdd2.setCoeff(j, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcdd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdd2(): order " + out.str() + "not implemented"); 
    }

    vmcd2.push_back(mcdd2);
    return(vmcd2);
}

  std::vector<WilsonCoefficient>& StandardModelMatching::CMK()  
 {
    
    double xt = x_t(Muw);
    
    vmck.clear();
    
    switch (mck.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mck.getScheme();
            throw "StandardModel::CMK(): scheme " + out.str() + "not implemented";
    }

    mck.setMu(Muw);
    
    switch (mck.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<10; j++){
                mck.setCoeff(j, lam_t * SM.Als(Muw, FULLNLO) / 4. / M_PI * //* CHECK ORDER *//
                                setWCbnlep(j, xt,  NLO), NLO);
                mck.setCoeff(j, lam_t * Ale / 4. / M_PI *
                                setWCbnlepEW(j, xt), NLO_QED11);
                }
        case LO:
            for (int j=0; j<10; j++){
                mck.setCoeff(j, lam_t *  setWCbnlep(j, xt,  LO), LO);
                mck.setCoeff(j, 0., LO_QED); 
                }                   
            break;
        default:
            std::stringstream out;
            out << mck.getOrder();
            throw "StandardModelMatching::CMK(): order " + out.str() + "not implemented";
    }

    vmck.push_back(mck);
    return(vmck);
}

/*******************************************************************************
 * Wilson coefficients Buras base for K -> pi pi decays                        * 
 * operator basis: - current current                                           *
 * ****************************************************************************/
  std::vector<WilsonCoefficient>& StandardModelMatching::CMKCC() 
 {
    
    double xt = x_t(Muw);
    
    vmckcc.clear();
    
    switch (mckcc.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mckcc.getScheme();
            throw "StandardModel::CMKCC(): scheme " + out.str() + "not implemented";
    }

    mckcc.setMu(Muw);
    
    switch (mckcc.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<2; j++){
                mckcc.setCoeff(j, lam_t * setWCbnlep(j, xt, NLO), NLO); 
            }
            for (int j=2; j<10; j++){
                mckcc.setCoeff(j, 0. , NLO); 
            }
        case LO:
            for (int j=0; j<2; j++){
                mckcc.setCoeff(j, lam_t * setWCbnlep(j, xt, LO), LO); 
            }
            for (int j=2; j<10; j++){
                mckcc.setCoeff(j, 0. , LO); 
            }
            break;
        default:
            std::stringstream out;
            out << mckcc.getOrder();
            throw "StandardModelMatching::CMKCC(): order " + out.str() + "not implemented";
    }

    vmckcc.push_back(mckcc);
    return(vmckcc);
}

    
/*******************************************************************************
 * Wilson coefficients misiak base for b -> s g                            * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *
 * ****************************************************************************/
 std::vector<WilsonCoefficient>& StandardModelMatching::CMbsg() 
{    
    double xt = x_t(Muw);
    
    vmcbsg.clear();
    
    switch (mcbsg.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcbsg.getScheme();
            throw std::runtime_error("StandardModel::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbsg.setMu(Muw);
    double alSo4pi = SM.Alstilde5(Muw);
    
    switch (mcbsg.getOrder()) {
        case NNLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, alSo4pi * alSo4pi * setWCbsg(j, xt,  NNLO) , NNLO);
            }
        case NLO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, alSo4pi * setWCbsg(j, xt,  NLO) , NLO);
            }
        case LO:
            for (int j=0; j<8; j++){
            mcbsg.setCoeff(j, setWCbsg(j, xt,  LO), LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcbsg.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }
    
    vmcbsg.push_back(mcbsg);
    return(vmcbsg);
}

 
 std::vector<WilsonCoefficient>& StandardModelMatching::CMprimebsg() 
{    
    vmcprimebsg.clear();
    
    switch (mcprimebsg.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcprimebsg.getScheme();
            throw std::runtime_error("StandardModel::CMprimebsg(): scheme " + out.str() + "not implemented"); 
    }

    mcprimebsg.setMu(Muw);
    
    switch (mcprimebsg.getOrder()) {
        case NNLO:
            for (int j=0; j<8; j++){
            mcprimebsg.setCoeff(j, 0., NNLO);//* CHECK ORDER *//
            }
        case NLO:
            for (int j=0; j<8; j++){
            mcprimebsg.setCoeff(j, 0., NLO);//* CHECK ORDER *//
            }
        case LO:
            for (int j=0; j<8; j++){
            mcprimebsg.setCoeff(j, 0., LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcprimebsg.getOrder();
            throw std::runtime_error("StandardModelMatching::CMprimebsg(): order " + out.str() + "not implemented"); 
    }
    
    vmcprimebsg.push_back(mcprimebsg);
    return(vmcprimebsg);
}

/*******************************************************************************
 * Wilson coefficients calcoulus, misiak base for b -> s gamma                  *  
 * ****************************************************************************/

double StandardModelMatching::setWCbsg(int i, double x, orders order)
{    
    sw =  sqrt( sW2 );

    if ( swf == sw && xcachef == x){
        switch (order){
        case NNLO:
            return ( CWbsgArrayNNLO[i] );
            break;        
        case NLO:
            return ( CWbsgArrayNLO[i] );
            break;
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
    }
    
    swf = sw; xcachef = x;

    switch (order){
        case NNLO:
            // One and two loop matching from arXiv:9910220.
            CWbsgArrayNNLO[0] = -Tt(x)+7987./72.+17./3.*M_PI*M_PI+475./6.*L+17.*L*L;
            CWbsgArrayNNLO[1] = 127./18.+4./3.*M_PI*M_PI+46./3.*L+4.*L*L;
            CWbsgArrayNNLO[2] = G1t(x,Muw)-680./243.-20./81.*M_PI*M_PI-68./81.*L-20./27*L*L;
            CWbsgArrayNNLO[3] = E1t(x,Muw)+950./243.+10./81.*M_PI*M_PI+124./27.*L+10./27.*L*L;
            CWbsgArrayNNLO[4] = -0.1*G1t(x,Muw)+2./15.*E0t(x)+68./243.+2./81.*M_PI*M_PI+14./81.*L+2./27.*L*L;
            CWbsgArrayNNLO[5] = -3./16.*G1t(x,Muw)+0.25*E0t(x)+85./162.+5./108.*M_PI*M_PI+35./108.*L+5./36*L*L;
            // 3-loop matching from arXiv:0401041. Expansion around x = 1, on the basis of Fig.3 of the paper. 
            CWbsgArrayNNLO[6] = (C7t_3L_at_mt(x) + C7t_3L_func(x,Muw)-(C7c_3L_at_mW(x)+ 13763./2187.*L+814./729.*L*L)) - 1./3.*CWbsgArrayNNLO[2] - 4./9.*CWbsgArrayNNLO[3] - 20./3.*CWbsgArrayNNLO[4] - 80./9.*CWbsgArrayNNLO[5]; //effective C7
            CWbsgArrayNNLO[7] = (C8t_3L_at_mt(x) + C8t_3L_func(x,Muw)-(C8c_3L_at_mW(x) + 16607./5832.*L+397./486.*L*L)) + CWbsgArrayNNLO[2]-1./6.*CWbsgArrayNNLO[3]-20.*CWbsgArrayNNLO[4]-10./3.*CWbsgArrayNNLO[5]; //effective C8
        case NLO:
            CWbsgArrayNLO[0] = 15. + 6*L;
            CWbsgArrayNLO[3] = E0t(x) - 7./9. + 2./3.* L;
            CWbsgArrayNLO[6] = -0.5*A1t(x,Muw) + 713./243. + 4./81.*L - 4./9.*CWbsgArrayNLO[3]; //effective C7
            CWbsgArrayNLO[7] = -0.5*F1t(x,Muw) + 91./324. - 4./27.*L - 1./6.*CWbsgArrayNLO[3]; //effective C8
        case LO:
            CWbsgArrayLO[1] = 1.;
            CWbsgArrayLO[6] = -0.5*A0t(x) - 23./36.; //effective C7
            CWbsgArrayLO[7] = -0.5*F0t(x) - 1./3.; //effective C8
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    switch (order){
        case NNLO:
            return ( CWbsgArrayNNLO[i] );
            break;
        case NLO:
            return ( CWbsgArrayNLO[i] );
            break;
        case LO:
            return ( CWbsgArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
}




/*******************************************************************************
 * Wilson coefficients misiak base for B -> K^*ll                              * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *         
 *                 - semileptonic                                              * 
 * ****************************************************************************/
  std::vector<WilsonCoefficient>& StandardModelMatching::CMBMll(QCD::lepton lepton) 
    {    
    double xt = x_t(Muw); //* ORDER FULLNNLO*//
    
    vmcBMll.clear();
    
    switch (mcBMll.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcBMll.getScheme();
            throw std::runtime_error("StandardModel::CMBKstrall(): scheme " + out.str() + "not implemented"); 
    }

    mcBMll.setMu(Muw);
    
    switch (mcBMll.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<13; j++){
            mcBMll.setCoeff(j, SM.Als(Muw, FULLNNLO) / 4. / M_PI * setWCBMll(j, xt,  NLO) , NLO);
            }
        case LO:
            for (int j=0; j<13; j++){
            mcBMll.setCoeff(j, setWCBMll(j, xt,  LO), LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcBMll.getOrder();
            throw std::runtime_error("StandardModelMatching::CMBKstrall(): order " + out.str() + "not implemented"); 
    }
   
    vmcBMll.push_back(mcBMll);
    return(vmcBMll);
}
        


 /*******************************************************************************
 * Wilson coefficients calcoulus, misiak base for B -> K^*ll                    *  
 * *****************************************************************************/

double StandardModelMatching::setWCBMll(int i, double x, orders order)
{    
    sw =  sqrt( sW2 ) ;

    if ( swa == sw && xcachea == x){
        switch (order){
        case NNLO:
        case NLO:
            return ( CWBMllArrayNLO[i] );
            break;
        case LO:
            return ( CWBMllArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
    }
    
    swa = sw; xcachea = x;

    switch (order){
        case NNLO:
        case NLO:
            CWBMllArrayNLO[0] = 15. + 6*L;
            CWBMllArrayNLO[3] = E0t(x) - (7./9.) + (2./3.* L);
            CWBMllArrayNLO[6] = -0.5*A1t(x,Muw) + 713./243. + 4./81.*L - 4./9.*CWBMllArrayNLO[3];
            CWBMllArrayNLO[7] = -0.5*F1t(x,Muw) + 91./324. - 4./27.*L - 1./6.*CWBMllArrayNLO[3];
            CWBMllArrayNLO[8] = (1-4.*sW2) / (sW2) *C1t(x,Muw) - 1./(sW2) * B1t(x,Muw) - D1t(x,Muw) + 1./sW2 + 524./729. - 
                    128.*M_PI*M_PI/243. - 16.*L/3. -128.*L*L/81.;
            CWBMllArrayNLO[9] = (B1t(x,Muw) - C1t(x,Muw)) / sW2 - 1./sW2;
        case LO:
            CWBMllArrayLO[1] = 1.;
            CWBMllArrayLO[6] = -0.5*A0t(x) - 23./36.;
            CWBMllArrayLO[7] = -0.5*F0t(x) - 1./3.;
            CWBMllArrayLO[8] = (1-4.*sW2) / (sW2) *C0t(x) - 1./(sW2) * B0t(x) - D0t(x) + 38./27. + 1./(4.*sW2) - (4./9.)*L + 8./9. * log(SM.getMuw()/mu_b); 
            CWBMllArrayLO[9] = 1./(sW2) * (B0t(x) - C0t(x)) -1./(4.*sW2);
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    switch (order){
        case NNLO:
        case NLO:
            return ( CWBMllArrayNLO[i] );
            break;
        case LO:
            return ( CWBMllArrayLO[i] );
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
}
    


/*******************************************************************************
 * Wilson coefficients misiak primed base for B -> K^*ll                       * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *         
 *                 - semileptonic                                              * 
 * ****************************************************************************/
 std::vector<WilsonCoefficient>& StandardModelMatching::CMprimeBMll(QCD::lepton lepton) 
    {
        vmcprimeBMll.clear();
        mcprimeBMll.setMu(Muw);
        switch (mcprimeBMll.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<13; j++){
            mcprimeBMll.setCoeff(j, 0., NLO);
            }
        case LO:
            for (int j=0; j<13; j++){
            mcprimeBMll.setCoeff(j, 0., LO);
            }
            break;
        default:
            std::stringstream out;
            out << mcprimeBMll.getOrder();
            throw std::runtime_error("StandardModelMatching::CMBKstrall(): order " + out.str() + "not implemented"); 
    }
        vmcprimeBMll.push_back(mcprimeBMll);
        return(vmcprimeBMll);
    }


    /******************************************************************************/

/*******************************************************************************
 * Wilson coefficients Misiak base for bs -> mu mu. decays                      * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - semileptonic                                              *
 * ****************************************************************************/
    
 std::vector<WilsonCoefficient>& StandardModelMatching::CMbsmm() 
{
    
         // The couplings are not used here, but in the Bsmumu class.
            
        double xt = x_t(Muw);
         
        vmcbsmm.clear();
    
    
        
        switch (mcbsmm.getScheme()) {
        case NDR:
             //case HV:
             //case LRI:
             break;
             default:
            std::stringstream out;
            out << mcbsmm.getScheme();
            throw std::runtime_error("StandardModel::CMbsmm(): scheme " + out.str() + "not implemented"); 
    }
       
        mcbsmm.setMu(Muw);
        
        switch (mcbsmm.getOrder()) {
        case NNLO:
        
                for (int j=0; j<8; j++){
                   
                    mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) *   
                                setWCBsmm(j, xt,  NNLO), NNLO);
                    
                }
        case NLO:
           
            for (int j=0; j<8; j++){
             
               mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * 
                                setWCBsmm(j, xt,  NLO), NLO);
                }
        case LO:
      
            for (int j=0; j<8; j++){
              
                mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * setWCBsmm(j, xt,  LO), LO);
                }                   
            break;
            default:
            std::stringstream out;
            out << mcbsmm.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsmm(): order " + out.str() + "not implemented");
        }
    
        switch (mcbsmm.getOrder_qed()) {
        case NLO_QED22:
          ;
            for (int j=0; j<8; j++){
               
                mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * setWCBsmmEW(j, xt, NLO_QED22), NLO_QED22);
                
            }
        case NLO_QED12:  /*absent at high energy */
            for (int j=0; j<8; j++){
            
            mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * 0., NLO_QED12); 
            }
        case NLO_QED21: 
            
            for (int j=0; j<8; j++){
              
            mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * setWCBsmmEW(j, xt, NLO_QED21), NLO_QED21); 
            }
        case NLO_QED02:   /*absent at high energy */
            for (int j=0; j<8; j++){
            
              mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * 0., NLO_QED02);
            }
        case NLO_QED11:
            for (int j=0; j<8; j++){
            
               mcbsmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,1)) * setWCBsmmEW(j, xt, NLO_QED11), NLO_QED11);
             
            }
        case LO_QED:   /*absent at high energy */
            for (int j=0; j<8; j++){
              
             mcbsmm.setCoeff(j,(Vckm(2,2).conjugate() * Vckm(2,1)) * 0., LO_QED); 
            }                   
            break;
            default:
            std::stringstream out;
            out << mcbsmm.getOrder_qed();
            throw std::runtime_error("StandardModelMatching::CMbsmm(): order_qed " + out.str() + "not implemented");
    }
    vmcbsmm.push_back(mcbsmm);
    return(vmcbsmm);
    
}
    
/*******************************************************************************
 * Wilson coefficients Misiak base for bd -> mu mu. decays                      * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - semileptonic                                              *
 * ****************************************************************************/
    
 std::vector<WilsonCoefficient>& StandardModelMatching::CMbdmm() 
{
    
         // The couplings are not used here, but in the Bdmumu class.
           
        double xt = x_t(Muw);
         
        vmcbdmm.clear();
    
    
        
        switch (mcbdmm.getScheme()) {
        case NDR:
             //case HV:
             //case LRI:
             break;
             default:
            std::stringstream out;
            out << mcbdmm.getScheme();
            throw std::runtime_error("StandardModel::CMbdmm(): scheme " + out.str() + "not implemented"); 
    }
       
        mcbdmm.setMu(Muw);
        
        switch (mcbdmm.getOrder()) {
        case NNLO:
        
                for (int j=0; j<8; j++){
                   
                    mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) *   
                                setWCBdmm(j, xt,  NNLO), NNLO);
                }
        case NLO:
           
            for (int j=0; j<8; j++){
             
               mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * 
                                setWCBdmm(j, xt,  NLO), NLO);
                }
        case LO:
      
            for (int j=0; j<8; j++){
              
                mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * setWCBdmm(j, xt,  LO), LO);
                }                   
            break;
            default:
            std::stringstream out;
            out << mcbdmm.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbdmm(): order " + out.str() + "not implemented");
        }
    
        switch (mcbdmm.getOrder_qed()) {
        case NLO_QED22:
          ;
            for (int j=0; j<8; j++){
               
                mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * setWCBdmmEW(j, xt, NLO_QED22), NLO_QED22);
                
            }
        case NLO_QED12:  /*absent at high energy */
            for (int j=0; j<8; j++){
            
            mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * 0., NLO_QED12); 
            }
        case NLO_QED21: 
            
            for (int j=0; j<8; j++){
              
            mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * setWCBdmmEW(j, xt, NLO_QED21), NLO_QED21); 
            }
        case NLO_QED02:   /*absent at high energy */
            for (int j=0; j<8; j++){
            
              mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * 0., NLO_QED02);
            }
        case NLO_QED11:
            for (int j=0; j<8; j++){
            
               mcbdmm.setCoeff(j, (Vckm(2,2).conjugate() * Vckm(2,0)) * setWCBdmmEW(j, xt, NLO_QED11), NLO_QED11);
             
            }
        case LO_QED:   /*absent at high energy */
            for (int j=0; j<8; j++){
              
             mcbdmm.setCoeff(j,(Vckm(2,2).conjugate() * Vckm(2,0)) * 0., LO_QED); 
            }                   
            break;
            default:
            std::stringstream out;
            out << mcbdmm.getOrder_qed();
            throw std::runtime_error("StandardModelMatching::CMbdmm(): order_qed " + out.str() + "not implemented");
    }
    vmcbdmm.push_back(mcbdmm);
    return(vmcbdmm);
    
}  
  
/*******************************************************************************
 * Wilson coefficients calcoulus, misiak base for B -> tau nu                   *
 * ****************************************************************************/

 std::vector<WilsonCoefficient>& StandardModelMatching::CMbtaunu() 
{
    
    vmcbtaunu.clear();
    
    mcbtaunu.setMu(Muw);
 
    switch (mcbtaunu.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcbtaunu.setCoeff(0, 4.*GF * Vckm(0,2) / sqrt(2.) , LO);
            break;
        default:
            std::stringstream out;
            out << mcbtaunu.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbtaunu(): order " + out.str() + "not implemented");
    }
    
    vmcbtaunu.push_back(mcbtaunu);
    return(vmcbtaunu);
    
}     

    
/******************************************************************************/

/*******************************************************************************
 * Wilson coefficients Buras base for b -> nonlep. decays                      * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *         
 *                 - semileptonic                                              *
 * i=0 deltaS=0 deltaC=0;  i=1 1,0 ;                                           *
 * ****************************************************************************/
 std::vector<WilsonCoefficient>& StandardModelMatching::CMbnlep(const int a)
{
    gslpp::complex lambda;
    
    switch (a) {
        case 0: lambda = SM.computelamt_d();
        break;
        case 1: lambda = SM.computelamt_s();
        break;   
        default:
            std::stringstream out;
            out << a;
            throw std::runtime_error("case" + out.str() + "not implemented; implemented i=0,1,2,3"); 
    }
    
    double xt = x_t(Muw);
    double co = ( GF / sqrt(2));
    
    vmcbnlep.clear();
    
    switch (mcbnlep.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcbnlep.getScheme();
            throw std::runtime_error("StandardModel::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbnlep.setMu(Muw);
    
    switch (mcbnlep.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<10; j++){
                mcbnlep.setCoeff(j,co * lambda * SM.Als(Muw, FULLNLO) / 4. / M_PI * //* CHECK ORDER *//
                                setWCbnlep(j, xt,  NLO), NLO);
                mcbnlep.setCoeff(j, co * lambda * Ale / 4. / M_PI *
                                setWCbnlepEW(j, xt), NLO_QED11);
                }
        case LO:
            for (int j=0; j<10; j++){
                mcbnlep.setCoeff(j, co * lambda *  setWCbnlep(j, xt,  LO), LO);
                mcbnlep.setCoeff(j, 0., LO_QED); 
                }                   
            break;
        default:
            std::stringstream out;
            out << mcbnlep.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbnlep.push_back(mcbnlep);
    return(vmcbnlep);
}

/*******************************************************************************
 * Wilson coefficients Buras base for b -> nonlep. decays                      * 
 * operator basis: - current current opertors                                  *         
 * i=0 deltaS=0 deltaC=0;  i=1 1,0 ;  i=2 0,1 ; i=3 1,1                        *
 * ****************************************************************************/
 std::vector<WilsonCoefficient>& StandardModelMatching::CMbnlepCC(const int a) 
{    
    gslpp::complex lambda1 = 0.;
    //gslpp::complex lambda2 = 0.;
    //gslpp::matrix<gslpp::complex> ckm = SM.getVCKM();
    
    switch (a) {
        case 0: lambda1 = SM.computelamu_d();
                break;
        case 1: lambda1 = SM.computelamu_s();
                break;
        case 2: lambda1 = Vckm(0,2).conjugate()*Vckm(1,0);
                break;
        case 3: lambda1 = Vckm(1,2).conjugate()*Vckm(0,0); 
                break;
        case 4: lambda1 = Vckm(0,2).conjugate()*Vckm(1,1);
                break;
        case 5: lambda1 = Vckm(1,2).conjugate()*Vckm(2,1);
                break;
        default:
            std::stringstream out;
            out << a;
            throw std::runtime_error("case" + out.str() + " not existing; implemented i=0,1,2,3"); 
    }
    
    double xt = x_t(Muw);
    double co = ( GF / sqrt(2));
    
    vmcbnlepCC.clear();
    
    switch (mcbnlepCC.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mcbnlepCC.getScheme();
            throw std::runtime_error("StandardModel::CMbsg(): scheme " + out.str() + "not implemented"); 
    }

    mcbnlepCC.setMu(Muw);
    
    switch (mcbnlepCC.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<2; j++){
                mcbnlepCC.setCoeff(j, co * lambda1 * setWCbnlep(j, xt,  NLO), NLO); 
            }
            for (int j=2; j<10; j++){
                mcbnlepCC.setCoeff(j, 0. , NLO); 
            }
        case LO:
            for (int j=0; j<2; j++){
                mcbnlepCC.setCoeff(j, co * lambda1 *  setWCbnlep(j, xt,  LO), LO); }
            for (int j=2; j<10; j++){
                mcbnlepCC.setCoeff(j, 0. , LO); }
            break;
        default:
            std::stringstream out;
            out << mcbnlepCC.getOrder();
            throw std::runtime_error("StandardModelMatching::CMbsg(): order " + out.str() + "not implemented"); 
    }

    vmcbnlepCC.push_back(mcbnlepCC);
    return(vmcbnlepCC);
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMkpnn() {
    
    //scales assigned to xt, a and Xewt to be checked!
    
    double xt = x_t(SM.getMut());
    double a = 1./mt2omh2(Muw);
    double lambda5 = SM.getLambda()*SM.getLambda()*SM.getLambda()*SM.getLambda()*SM.getLambda();
    
    vmckpnn.clear();
    
    mckpnn.setMu(Mut);
 
    switch (mckpnn.getOrder()) {
        case NNLO:
        case NLO:
            mckpnn.setCoeff(0, SM.Als(SM.getMut(), FULLNLO)/4./M_PI*lam_t.imag()*X1t(xt)/lambda5, NLO);//* CHECK ORDER *//
        case LO:
            mckpnn.setCoeff(0, lam_t.imag()*X0t(xt)/lambda5, LO);
            break;
        default:
            std::stringstream out;
            out << mckpnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkpnn(): order " + out.str() + "not implemented"); 
    }
    
    switch (mckpnn.getOrder_qed()) {
        case NLO_QED11:
            mckpnn.setCoeff(0, Ale/4./M_PI*lam_t.imag()*Xewt(xt, a, Muw)/lambda5, NLO_QED11);
        case LO_QED:
            break; 
        default:
            std::stringstream out;
            out << mckpnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkpnn(): order " + out.str() + "not implemented"); 
    }

    vmckpnn.push_back(mckpnn);
    return(vmckpnn);
    
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMkmm() {
    
    //PROBLEMI: mu e sin(theta_weak) e la scala di als
    
    double xt = x_t(Muw);
    
    vmckmm.clear();
    
    mckmm.setMu(Mut);
 
    switch (mckmm.getOrder()) {
        case NNLO:
        case NLO:
            mckmm.setCoeff(0, SM.Als(Muw, FULLNLO)/4./M_PI*lam_t.real()*Y1(xt, Muw)/SM.getLambda(), NLO);//* CHECK ORDER *//
        case LO:
            mckmm.setCoeff(0, lam_t.real()*Y0(xt)/SM.getLambda(), LO);
            break;
        default:
            std::stringstream out;
            out << mckmm.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkmm(): order " + out.str() + "not implemented"); 
    }

    vmckmm.push_back(mckmm);
    return(vmckmm);
    
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMBXsnn() {
    
    double xt = x_t(Muw);
    
    vmcbsnn.clear();
    
    mcbsnn.setMu(Mut);
 
    switch (mcbsnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff(0, (Vckm(2,1).abs() / Vckm(1,2).abs())*
                                SM.Als(Muw, FULLNLO)/4./M_PI*X1t(xt), NLO);//* CHECK ORDER *//
        case LO:
            mcbsnn.setCoeff(0,  (Vckm(2,1).abs() / Vckm(1,2).abs())*
                                X0t(xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcbsnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMXsnn(): order " + out.str() + "not implemented"); 
    }

    vmcbsnn.push_back(mcbsnn);
    return(vmcbsnn);
    
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMBXdnn() {
    
    double xt = x_t(Muw);
    
    vmcbdnn.clear();
    
    mcbdnn.setMu(Mut);
 
    switch (mcbdnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff(0, (Vckm(2,2).abs() / Vckm(1,2).abs()) *
                                SM.Als(Muw, FULLNLO)/4./M_PI*X1t(xt), NLO);//* CHECK ORDER *//
        case LO:
            mcbsnn.setCoeff(0,  (Vckm(2,2).abs() / Vckm(1,2).abs()) *
                                X0t(xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcbdnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CXdnn(): order " + out.str() + "not implemented"); 
    }
    
    vmcbdnn.push_back(mcbdnn);
    return(vmcbdnn);
    
}

/*******************************************************************************
 * Wilson coefficients calculus, MISIAK base for Bs to mu mu  decay               *  
 * ****************************************************************************/
double StandardModelMatching::setWCBsmm(int i, double x, orders order) 	
{  
    
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2.) * GF * Mw * Mw) );
     
    if ( swd == sw && xcached == x){
        switch (order){
        case NNLO:
           return (CWBsmmArrayNNLOqcd[i]);
           break;                               
       case NLO:
            return (CWBsmmArrayNLOqcd[i]);
            break;
        case LO:
            return (CWBsmmArrayLOqcd[i]);
            break;
       default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted");      
         }
    }
  
    swd = sw;  xcached = x;
 
    switch (order){
    case NNLO:
        CWBsmmArrayNNLOqcd[0] = sw * sw * (-Tt(x) + 7987./72. + 17. * M_PI * M_PI/3. + 475. * L/6. + 17. * L * L);
        CWBsmmArrayNNLOqcd[1] = sw * sw * (127./18. + 4. * M_PI * M_PI /3. + 46. * L/3. + 4. * L * L);
        CWBsmmArrayNNLOqcd[2] = sw * sw * (G1t(x, Muw) - 680./243. - 20. * M_PI * M_PI /81. - 68. * L/81. - 20. * L* L/27.);
        CWBsmmArrayNNLOqcd[3] = sw * sw * (E1t(x, Muw) + 950./243. + 10.* M_PI * M_PI /81. + 124. * L/27. + 10. * L * L/27.);   
        CWBsmmArrayNNLOqcd[4] = sw * sw * (-G1t(x, Muw)/10. + 2. * E0t(x)/15. + 68./243. + 2. * M_PI * M_PI /81. + 14.* L/81. + 2. * L * L/27.);
        CWBsmmArrayNNLOqcd[5] = sw * sw * (-3. * G1t(x, Muw)/16. + E0t(x)/4. + 85./162. + 5. * M_PI * M_PI/108. + 35. * L/108. + 5. * L * L/36.);  
               
    case NLO:
        CWBsmmArrayNLOqcd[0] = sw * sw * (15. + 6. * L);
        CWBsmmArrayNLOqcd[3] = sw * sw * (Eet(x) - 2./3. + 2. * L/3.);
       
    case LO:
        CWBsmmArrayLOqcd[1] = sw * sw * 1.;
        
    break;
    default:
    std::stringstream out;
    out << order;
    throw std::runtime_error("order" + out.str() + "not implemeted"); 
    }
    switch (order){
    case NNLO:
        return (CWBsmmArrayNNLOqcd[i]);
       
        break;
    case NLO:
        return (CWBsmmArrayNLOqcd[i]);
       
        break;
    case LO:
        return (CWBsmmArrayLOqcd[i]);
        
        break;
        default:
        std::stringstream out;
        out << order;
        throw std::runtime_error("order" + out.str() + "not implemeted");      
    }
}

double StandardModelMatching::setWCBsmmEW(int i, double x, orders_qed order_qed) 	
{   
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2.) * GF * Mw * Mw) ) ;

    double mt = SM.Mrun(Muw, SM.getQuarks(QCD::TOP).getMass_scale(), 
                SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
            
    if ( swe == sw && xcachee == x){
        switch (order_qed){
        case NLO_QED22:
            return (CWBsmmArrayNLOewt4[i]);   
            break;
        case NLO_QED21:
            return (CWBsmmArrayNLOewt2[i]);
            break;
        case NLO_QED11:
            return (CWBsmmArrayNLOew[i]);
            break;
            default:
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("order_qed" + out.str() + "not implemeted");      
        }

    }
   
    
    swe = sw; xcachee = x;
    
///*non-student*/double xz = SM.getMz() * SM.getMz() / Mw / Mw;
///*non-student*/double xht = SM.getMHl()*SM.getMHl()/Mt_muw/Mt_muw;
 
    switch (order_qed){   
    case NLO_QED22:
///*non-student*/CWBsmmArrayNLOewt4[7] = sw * sw * (-x*x/32./sw/sw/sw/sw*(3.+taub2(xht)-Delta_t(Muw,xht))) ;
///*non-student*/CWBsmmArrayNLOewt4[6] = (4. * sw * sw - 1.) * CWBsmmArrayNLOewt4[7];
        CWBsmmArrayNLOewt4[7] = sw * sw * (1./(sw * sw)) * Rest(x, Muw) ;
//        std::cout << ">>>>>>> " << Rest(163.5*163.5/80.358/80.358, 80.358) << std::endl;
        
    case NLO_QED21:
///*non-student*/CWBsmmArrayNLOewt2[2] = sw * sw * (1. / sw / sw * (4. / 9. * B1d(x, Muw) + 4. / 27. * B1d_tilde(x, Muw) + 2. / 9. * B1u(x, Muw) +
//                    2. / 27. * B1u_tilde(x, Muw) - 2. / 9. * C1ew(x) + 320. / 27. * B0b(x) + 160. / 27. * C0b(x)));
///*non-student*/CWBsmmArrayNLOewt2[3] = sw * sw * (16. / 27. * C0b(x) + 1. / sw /sw * (8. / 9. * B1d_tilde(x, Muw) + 4. / 9. * B1u_tilde(x, Muw) -
//                    2. / 9. * Gew(x, xz, Muw) - 88. / 9. * B0b(x) - 184. / 27. * C0b(x)));
///*non-student*/CWBsmmArrayNLOewt2[4] = sw * sw * (1. / sw / sw * (-1. / 9. * B1d(x, Muw) - 1. / 27. * B1d_tilde(x, Muw) - 1. / 18. * B1u(x, Muw) -
//                    1. / 54. * B1u_tilde(x, Muw) + 1. / 18. * C1ew(x) - 32. / 27. * B0b(x) - 16. / 27. * C0b(x)));
///*non-student*/CWBsmmArrayNLOewt2[5] = sw * sw * (1. / sw / sw * (-2. / 9. * B1d_tilde(x, Muw) - 1. / 9. * B1u_tilde(x, Muw) + 1. / 18. * Gew(x, xz, Muw) +
//                    4. / 3. * B0b(x) + 2. / 3. * C0b(x)));
        CWBsmmArrayNLOewt2[6] = sw * sw * ((1. - 4. * sw * sw) * C1t(x, Muw) / (sw * sw) - B1t(x, Muw)/(sw * sw) 
                - D1t(x, Muw) + 1./ (sw * sw) + 524./729. - 128. * M_PI * M_PI / 243. 
                - 16. * L / 3. - 128. * L * L /81. ) ; 
        CWBsmmArrayNLOewt2[7] = sw * sw * ((1./(sw * sw)) * (B1t(x, Muw) - C1t(x, Muw)) - 1./(sw * sw)) ; 
     
    case NLO_QED11:
///*non-student*/CWBsmmArrayNLOew[1] = sw * sw * (-22. / 9. - 4. / 3. * Lz + 1. / 9.);
///*non-student*/CWBsmmArrayNLOew[2] = sw * sw * (-2. / 9. / sw /sw * (2. * B0b(x) + C0b(x)));
///*non-student*/CWBsmmArrayNLOew[4] = sw * sw * (1. / 9. / sw /sw * (B0b(x) + 1. / 2. * C0b(x)));
        CWBsmmArrayNLOew[6] = sw * sw * (Y0(x)/(sw * sw) + Wt(x) + 4./9. - 4. * 2 * log(Muw/mt)/9.);
        CWBsmmArrayNLOew[7] = sw * sw * (-Y0(x)/(sw * sw));   
            
        break;
        default:
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("order_qed" + out.str() + "not implemeted"); 
    }

    switch (order_qed){
    case NLO_QED22:
        return (CWBsmmArrayNLOewt4[i]);   
        break;  
    case NLO_QED21:
        return (CWBsmmArrayNLOewt2[i]);
        break;  
    case NLO_QED11:
        return (CWBsmmArrayNLOew[i]);
        break;
        default:
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("order_qed" + out.str() + "not implemeted");  
    } 
}




/*******************************************************************************
 * Wilson coefficients calculus, MISIAK base for Bd to mu mu  decay               *  
 * ****************************************************************************/
double StandardModelMatching::setWCBdmm(int i, double x, orders order) 	
{  
    
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2.) * GF * Mw * Mw) );
     
    if ( swb == sw && xcacheb == x){
        switch (order){
        case NNLO:
           return (CWBdmmArrayNNLOqcd[i]);
           break;                               
       case NLO:
            return (CWBdmmArrayNLOqcd[i]);
            break;
        case LO:
            return (CWBdmmArrayLOqcd[i]);
            break;
       default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted");      
         }
    }
  
    swb = sw;  xcacheb = x;
 
    switch (order){
    case NNLO:
        CWBdmmArrayNNLOqcd[0] = sw * sw * (-Tt(x) + 7987./72. + 17. * M_PI * M_PI/3. + 475. * L/6. + 17. * L * L);
        CWBdmmArrayNNLOqcd[1] = sw * sw * (127./18. + 4. * M_PI * M_PI /3. + 46. * L/3. + 4. * L * L);
        CWBdmmArrayNNLOqcd[2] = sw * sw * (G1t(x, Muw) - 680./243. - 20. * M_PI * M_PI /81. - 68. * L/81. - 20. * L* L/27.);
        CWBdmmArrayNNLOqcd[3] = sw * sw * (E1t(x, Muw) + 950./243. + 10.* M_PI * M_PI /81. + 124. * L/27. + 10. * L * L/27.);   
        CWBdmmArrayNNLOqcd[4] = sw * sw * (-G1t(x, Muw)/10. + 2. * E0t(x)/15. + 68./243. + 2. * M_PI * M_PI /81. + 14.* L/81. + 2. * L * L/27.);
        CWBdmmArrayNNLOqcd[5] = sw * sw * (-3. * G1t(x, Muw)/16. + E0t(x)/4. + 85./162. + 5. * M_PI * M_PI/108. + 35. * L/108. + 5. * L * L/36.);  
               
    case NLO:
        CWBdmmArrayNLOqcd[0] = sw * sw * (15. + 6. * L);
        CWBdmmArrayNLOqcd[3] = sw * sw * (Eet(x) - 2./3. + 2. * L/3.);
       
    case LO:
        CWBdmmArrayLOqcd[1] = sw * sw * 1.;
        
    break;
    default:
    std::stringstream out;
    out << order;
    throw std::runtime_error("order" + out.str() + "not implemeted"); 
    }
    switch (order){
    case NNLO:
        return (CWBdmmArrayNNLOqcd[i]);
       
        break;
    case NLO:
        return (CWBdmmArrayNLOqcd[i]);
       
        break;
    case LO:
        return (CWBdmmArrayLOqcd[i]);
        
        break;
        default:
        std::stringstream out;
        out << order;
        throw std::runtime_error("order" + out.str() + "not implemeted");      
    }
}

double StandardModelMatching::setWCBdmmEW(int i, double x, orders_qed order_qed) 	
{   
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2.) * GF * Mw * Mw) ) ;

    double mt = SM.Mrun(Muw, SM.getQuarks(QCD::TOP).getMass_scale(), 
                SM.getQuarks(QCD::TOP).getMass(), FULLNNLO);
            
    if ( swc == sw && xcachec == x){
        switch (order_qed){
        case NLO_QED22:
            return (CWBdmmArrayNLOewt4[i]);   
            break;
        case NLO_QED21:
            return (CWBdmmArrayNLOewt2[i]);
            break;
        case NLO_QED11:
            return (CWBdmmArrayNLOew[i]);
            break;
            default:
            std::stringstream out;
            out << order_qed;
            throw std::runtime_error("order_qed" + out.str() + "not implemeted");      
        }

    }
   
    
    swc = sw; xcachec = x;
 
    switch (order_qed){   
    case NLO_QED22: 
        CWBdmmArrayNLOewt4[7] = sw * sw * (1./(sw * sw)) * Rest(x, Muw) ;
        
    case NLO_QED21: 
        CWBdmmArrayNLOewt2[6] = sw * sw * ((1. - 4. * sw * sw) * C1t(x, Muw) / (sw * sw) - B1t(x, Muw)/(sw * sw) 
                - D1t(x, Muw) + 1./ (sw * sw) + 524./729. - 128. * M_PI * M_PI / 243. 
                - 16. * L / 3. - 128. * L * L /81. ) ; 
        CWBdmmArrayNLOewt2[7] = sw * sw * ((1./(sw * sw)) * (B1t(x, Muw) - C1t(x, Muw)) - 1./(sw * sw)) ; 
     
    case NLO_QED11:
        CWBdmmArrayNLOew[6] = sw * sw * (Y0(x)/(sw * sw) + Wt(x) + 4./9. - 4. * 2 * log(Muw/mt)/9.);
        CWBdmmArrayNLOew[7] = sw * sw * (-Y0(x)/(sw * sw));   
            
        break;
        default:
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("order_qed" + out.str() + "not implemeted"); 
    }

    switch (order_qed){
    case NLO_QED22:
        return (CWBdmmArrayNLOewt4[i]);   
        break;  
    case NLO_QED02:
        return (CWBdmmArrayNLOewt2[i]);
        break;  
    case NLO_QED11:
        return (CWBdmmArrayNLOew[i]);
        break;
        default:
        std::stringstream out;
        out << order_qed;
        throw std::runtime_error("order_qed" + out.str() + "not implemeted");  
    } 
}


/*******************************************************************************
 * Wilson coefficients calculus, Buras base for nonlep. b decays               *  
 * ****************************************************************************/
double StandardModelMatching::setWCbnlep(int i, double x, orders order) 
{    
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2) * GF * Mw * Mw) );
    
    if ( swb == sw && xcacheb == x){
        switch (order){
        case NNLO:
        case NLO:
            return (CWbnlepArrayNLOqcd[i]);
            break;
        case LO:
            return (CWbnlepArrayLOqcd[i]);
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted");      
        }
    }
    
    swb = sw; xcacheb = x;
    
    switch (order){
        case NNLO:
        case NLO:
            CWbnlepArrayNLOqcd[0] = 11./2.;
            CWbnlepArrayNLOqcd[1] = -11./6.;
            CWbnlepArrayNLOqcd[2] = -1./6. * (E0b(x) - 2./3.);
            CWbnlepArrayNLOqcd[3] = 0.5 * (E0b(x) - 2./3.);
            CWbnlepArrayNLOqcd[4] = -1./6. * (E0b(x) - 2./3.);
            CWbnlepArrayNLOqcd[5] = 0.5 * (E0b(x) - 2./3.);
        case LO:
            CWbnlepArrayLOqcd[1] = 1.;
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
        }
    
    switch (order){
        case NNLO:
        case NLO:
            return (CWbnlepArrayNLOqcd[i]);
            break;
        case LO:
            return (CWbnlepArrayLOqcd[i]);
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted");      
        }
}

double StandardModelMatching::setWCbnlepEW(int i, double x) 
{     
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2) * GF * Mw * Mw) ) ;
    
    if ( swb == sw && xcacheb == x){
        return (CWbnlepArrayNLOew[i]);
    }
    
    swc = sw; xcachec = x;
    
    CWbnlepArrayNLOew[1] = -35./18.;
    CWbnlepArrayNLOew[2] = 2. / (3. * sw * sw) * ( 2. * B0b(x) + C0b(x) );
    CWbnlepArrayNLOew[6] = 2./3. * (4. * C0b(x) + D0b(x) - 4./9.);
    CWbnlepArrayNLOew[8] = 2./3. * (4. * C0b(x) + D0b(x) - 4./9. + (1. / (sw * sw)) *
                           (10. * B0b(x) - 4. * C0b(x)) );
    
    return (CWbnlepArrayNLOew[i]);
}

gslpp::complex StandardModelMatching::S0c() const 
{
    double xc = x_c(SM.getMuc());
    gslpp::complex co = GF / 2. / M_PI * Mw_tree * SM.computelamc().conjugate(); /* Mw_tree...?? */
#if SUSYFIT_DEBUG & 2
    std::cout << "im lambdac = " << (SM.computelamc()*SM.computelamc()).imag() << std::endl;
#endif
    return(co * co * S0(xc, xc));
}

gslpp::complex StandardModelMatching::S0ct() const 
{
    double xc = SM.Mrun4(SM.getMuc(),SM.getQuarks(QCD::CHARM).getMass_scale(),SM.getQuarks(QCD::CHARM).getMass())/Mw;
    xc *= xc;
    double xt = SM.Mrun4(SM.getMuw(),SM.getQuarks(QCD::TOP).getMass_scale(),SM.getQuarks(QCD::TOP).getMass())/Mw;
    xt *= xt;
    double co = GF / 2. / M_PI * Mw;
#if SUSYFIT_DEBUG & 2
    std::cout << "im lamc lamt = " << (SM.computelamc()*SM.computelamt()).imag() << std::endl;
#endif
    
    return( co * co * 2. * SM.computelamc().conjugate() * lam_t.conjugate() * S0(xc, xt) );
}

gslpp::complex StandardModelMatching::S0tt() const
{
    double xt = x_t(Mut);
    gslpp::complex co = GF / 2. / M_PI * Mw * lam_t.conjugate();
#if SUSYFIT_DEBUG & 2
    double pino = SM.Mrun(Mut, SM.Mp2Mbar(SM.getMtpole(),FULLNLO), 
                        SM.Mp2Mbar(SM.getMtpole(),FULLNLO), FULLNLO);
    std::cout << "mt(" << Mut<< ")" << pino << std::endl;
    double poldo = pino*pino/SM.Mw()/SM.Mw() ;
    std::cout << "S0(" << poldo << ") = " << S0(poldo,poldo) << std::endl;
    std::cout << "S0(" << xt << ") = " << S0(xt,xt) << std::endl;
    std::cout << "im lamt = " << (SM.computelamt()*SM.computelamt()).imag() << std::endl;
#endif

    return ( co * co * S0(xt, xt) );
}

/*******************************************************************************
 * Wilson coefficients for Lepton Flavour Violation               *  
 * ****************************************************************************/

 std::vector<WilsonCoefficient>& StandardModelMatching::CMDLij(int li_lj) 
{
    
    vmcDLij.clear();
    
    mcDLij.setMu(Muw);
    
    switch (mcDLij.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcDLij.setCoeff(0, 0., LO);
            mcDLij.setCoeff(1, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcDLij.getOrder();
            throw std::runtime_error("StandardModelMatching::CMDLij(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }
    
    vmcDLij.push_back(mcDLij);
    return(vmcDLij);
    
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMDLi3j(int li_lj) 
{

    vmcDLi3j.clear();

    mcDLi3j.setMu(Muw);

    switch (mcDLi3j.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcDLi3j.setCoeff(0, 0., LO);
            mcDLi3j.setCoeff(1, 0., LO);
            mcDLi3j.setCoeff(2, 0., LO);
            mcDLi3j.setCoeff(3, 0., LO);
            mcDLi3j.setCoeff(4, 0., LO);
            mcDLi3j.setCoeff(5, 0., LO);
            mcDLi3j.setCoeff(6, 0., LO);
            mcDLi3j.setCoeff(7, 0., LO);
            mcDLi3j.setCoeff(8, 0., LO);
            mcDLi3j.setCoeff(9, 0., LO);
            mcDLi3j.setCoeff(10, 0., LO);
            mcDLi3j.setCoeff(11, 0., LO);
            mcDLi3j.setCoeff(12, 0., LO);
            mcDLi3j.setCoeff(13, 0., LO);
            mcDLi3j.setCoeff(14, 0., LO);
            mcDLi3j.setCoeff(15, 0., LO);
            mcDLi3j.setCoeff(16, 0., LO);
            mcDLi3j.setCoeff(17, 0., LO);
            mcDLi3j.setCoeff(18, 0., LO);
            mcDLi3j.setCoeff(19, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcDLi3j.getOrder();
            throw std::runtime_error("StandardModelMatching::CMDLi3j(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }
    
    vmcDLi3j.push_back(mcDLi3j);
    return(vmcDLi3j);
    
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMmueconv() 
{
    
    vmcmueconv.clear();
    
    mcmueconv.setMu(Muw);
    
    switch (mcmueconv.getOrder()) {
        case NNLO:
        case NLO:
        case LO:
            mcmueconv.setCoeff(0, 0., LO);
            mcmueconv.setCoeff(1, 0., LO);
            mcmueconv.setCoeff(2, 0., LO);
            mcmueconv.setCoeff(3, 0., LO);
            mcmueconv.setCoeff(4, 0., LO);
            mcmueconv.setCoeff(5, 0., LO);
            mcmueconv.setCoeff(6, 0., LO);
            mcmueconv.setCoeff(7, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcmueconv.getOrder();
            throw std::runtime_error("StandardModelMatching::CMmueconv(): order " + out.str() + " not implemented.\nFor lepton flavour violating observables only Leading Order (LO) necessary.");
    }
    
    vmcmueconv.push_back(mcmueconv);
    return(vmcmueconv);
    
}

 std::vector<WilsonCoefficient>& StandardModelMatching::CMgminus2mu() 
{
    
    vmcgminus2mu.clear();
    
    mcgminus2mu.setMu(Muw);
    
    switch (mcgminus2mu.getOrder()) {
        case NNLO:
        case NLO:
            mcgminus2mu.setCoeff(0, 0., NLO);
            mcgminus2mu.setCoeff(1, 0., NLO);
            break;
        case LO:
            mcgminus2mu.setCoeff(0, 0., LO);
            mcgminus2mu.setCoeff(1, 0., LO);
            break;
        default:
            std::stringstream out;
            out << mcgminus2mu.getOrder();
            throw std::runtime_error("StandardModelMatching::CMgminus2mu(): order " + out.str() + " not implemented.");
    }
    
    vmcgminus2mu.push_back(mcgminus2mu);
    return(vmcgminus2mu);
    
}


WilsonCoefficient& StandardModelMatching::mc_C()
{
    switch (mcC.getScheme()) {
        case NDR:
            //case HV:
            //case LRI:
            break;
        default:
            std::stringstream out;
            out << mcC.getScheme();
            throw "StandardModel::mc_C(): scheme " + out.str() + "not implemented";
    }

    mcC.setMu(Muw); // cleared too

    switch (mcC.getOrder()) {
        case NNLO:
            mcC.setCoeff(0, alstilde * alstilde * (-Tt(x_t(Muw)) + 7987. / 72. + 17. / 3. * M_PI * M_PI +
                    L * 475. / 6. + 17. * L * L), NNLO);
            mcC.setCoeff(1, alstilde * alstilde * (127. / 18. +  4. / 3. * M_PI * M_PI  + 46. / 3. * L + 4. * L * L),
                    NNLO);
        case NLO:
            mcC.setCoeff(0, alstilde * (15. + 6. * L), NLO);
        case LO:
            mcC.setCoeff(1, 1., LO);
            break;
        default:
            std::stringstream out;
            out << mcC.getOrder();
            throw "StandardModelMatching::mc_C(): order " + out.str() + "not implemented";
    }

    switch (mcC.getOrder_qed()) {
        case NLO_QED22:
        case NLO_QED12:
        case NLO_QED21:
        case NLO_QED02:            
        case NLO_QED11:
            mcC.setCoeff(1, aletilde * (-22. / 9. - 4. / 3. * Lz + 1. / 9.), NLO_QED11);
        case LO_QED:
            break;
        default:
            std::stringstream out;
            out << mcC.getOrder_qed();
            throw "StandardModelMatching::mc_C(): order " + out.str() + "not implemented";
    }

    return(mcC);
}

double StandardModelMatching::C3funNNLO(double x)
{
    return(G1t(x,Muw)-680./243.-20./81.*M_PI*M_PI-68./81.*L-20./27*L*L);
}

double StandardModelMatching::C4fun(double x, orders ord)
{
    switch (ord) {
        case NNLO:
            return(E1t(x,Muw)+950./243.+10./81.*M_PI*M_PI+124./27.*L+10./27.*L*L);
        case NLO:
            return(E0t(x)- 7./9.+2./3.* L);
        default:
            std::stringstream out;
            out << ord;
            throw "StandardModelMatching::C4fun(): order " + out.str() + "not implemented";

    }
}

double StandardModelMatching::C5funNNLO(double x)
{
    return(-0.1*G1t(x,Muw)+2./15.*E0t(x)+68./243.+2./81.*M_PI*M_PI+14./81.*L+2./27.*L*L);
}

double StandardModelMatching::C6funNNLO(double x)
{
    return(-3./16.*G1t(x,Muw)+0.25*E0t(x)+85./162.+5./108.*M_PI*M_PI+35./108.*L+5./36*L*L);
}

WilsonCoefficient& StandardModelMatching::mc_P()
{
//    double sW2 = (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw() ); // ******* FOR TEST *********
    double xt = x_t(Muw);
    double xz = SM.getMz() * SM.getMz() / Mw / Mw;

    switch (mcP.getScheme()) {
        case NDR:
            //case HV:
            //case LRI:
            break;
        default:
            std::stringstream out;
            out << mcP.getScheme();
            throw "StandardModel::mc_P(): scheme " + out.str() + "not implemented";
    }

    mcP.setMu(Muw); // cleared too

    switch (mcP.getOrder()) {
        case NNLO:
            mcP.setCoeff(0, alstilde * alstilde * C3funNNLO(xt), NNLO);
            mcP.setCoeff(1, alstilde * alstilde * C4fun(xt, NNLO), NNLO);
            mcP.setCoeff(2, alstilde * alstilde * C5funNNLO(xt), NNLO);
            mcP.setCoeff(3, alstilde * alstilde * C6funNNLO(xt), NNLO);
        case NLO:
            mcP.setCoeff(1, alstilde * C4fun(xt, NLO), NLO);
        case LO:
            break;
        default:
            std::stringstream out;
            out << mcP.getOrder();
            throw "StandardModelMatching::mc_P(): order " + out.str() + "not implemented";
    }

    switch (mcP.getOrder_qed()) {
        case NLO_QED22:
        case NLO_QED12:
        case NLO_QED02:            
        case NLO_QED21:
            mcP.setCoeff(0, aletilde * alstilde * (1. / sW2 * (4. / 9. * B1d(xt, Muw) + 4. / 27. * B1d_tilde(xt, Muw) + 2. / 9. * B1u(xt, Muw) +
                    2. / 27. * B1u_tilde(xt, Muw) - 2. / 9. * C1ew(xt) + 320. / 27. * B0b(xt) +
                    160. / 27. * C0b(xt))), NLO_QED21);
            mcP.setCoeff(1, aletilde * alstilde * (16. / 27. * C0b(xt) + 1. / sW2 * (8. / 9. * B1d_tilde(xt, Muw) + 4. / 9. * B1u_tilde(xt, Muw) -
                    2. / 9. * Gew(xt, xz, Muw) - 88. / 9. * B0b(xt) - 184. / 27. * C0b(xt))), NLO_QED21);
            mcP.setCoeff(2, aletilde * alstilde * (1. / sW2 * (-1. / 9. * B1d(xt, Muw) - 1. / 27. * B1d_tilde(xt, Muw) - 1. / 18. * B1u(xt, Muw) -
                    1. / 54. * B1u_tilde(xt, Muw) + 1. / 18. * C1ew(xt) - 32. / 27. * B0b(xt) -
                    16. / 27. * C0b(xt))), NLO_QED21);
            mcP.setCoeff(3, aletilde * alstilde * (1. / sW2 * (-2. / 9. * B1d_tilde(xt, Muw) - 1. / 9. * B1u_tilde(xt, Muw) + 1. / 18. * Gew(xt, xz, Muw) +
                    4. / 3. * B0b(xt) + 2. / 3. * C0b(xt))), NLO_QED21);
        case NLO_QED11:
            mcP.setCoeff(0, aletilde * (-2. / 9. / sW2 * (2. * B0b(xt) + C0b(xt))), NLO_QED11);
//            mcP.setCoeff(0, aletilde*(-2./9./sW2*(2.*Y0(xt) - X0t(xt))), NLO_QED11);
            mcP.setCoeff(2, aletilde * (1. / 9. / sW2 * (B0b(xt) + 1. / 2. * C0b(xt))), NLO_QED11);
//            mcP.setCoeff(2, aletilde*(1./9./sW2*(Y0(xt) - 1./2.*X0t(xt))), NLO_QED11);
        case LO_QED:
            break;
        default:
            std::stringstream out;
            out << mcP.getOrder_qed();
            throw "StandardModelMatching::mc_P(): order " + out.str() + "not implemented";

    }
    
    return (mcP);
}


double StandardModelMatching::C7funLO(double x)
{
    return(-0.5*A0t(x) - 23./36.);
}

double StandardModelMatching::C8funLO(double x)
{
    return(-0.5*F0t(x) - 1./3.);
}

WilsonCoefficient& StandardModelMatching::mc_M()
{
    double xt = x_t(Muw);
    double mH = SM.getMHl();

    switch (mcM.getScheme()) {
        case NDR:
            //case HV:
            //case LRI:
            break;
        default:
            std::stringstream out;
            out << mcM.getScheme();
            throw "StandardModel::mc_M(): scheme " + out.str() + "not implemented";
    }

    mcM.setMu(Muw); // cleared too

    switch (mcM.getOrder()) {
        case NNLO:
            mcM.setCoeff(0, alstilde * alstilde * (C7t_3L_at_mt(xt) + C7t_3L_func(xt, Muw)-(C7c_3L_at_mW(xt) + 13763. / 2187. * L + 814. / 729. * L * L)
                    - 1. / 3. * C3funNNLO(xt) - 4. / 9. * C4fun(xt, NNLO) - 20. / 3. * C5funNNLO(xt) - 80. / 9. * C6funNNLO(xt)), NNLO);
            mcM.setCoeff(1, alstilde * alstilde * (C8t_3L_at_mt(xt) + C8t_3L_func(xt, Muw)-(C8c_3L_at_mW(xt) + 16607. / 5832. * L + 397. / 486. * L * L)
                    + C3funNNLO(xt) - 1. / 6. * C4fun(xt, NNLO) - 20. * C5funNNLO(xt) - 10. / 3. * C6funNNLO(xt)), NNLO);
        case NLO:
            mcM.setCoeff(0, alstilde * (-0.5 * A1t(xt, Muw) + 713. / 243. + 4. / 81. * L - 4. / 9. * C4fun(xt, NLO)), NLO);
            mcM.setCoeff(1, alstilde * (-0.5 * F1t(xt, Muw) + 91. / 324. - 4. / 27. * L - 1. / 6. * C4fun(xt, NLO)), NLO);
        case LO:
            mcM.setCoeff(0, C7funLO(xt), LO);
            mcM.setCoeff(1, C8funLO(xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcM.getOrder();
            throw "StandardModelMatching::mc_M(): order " + out.str() + "not implemented";
    }

    switch (mcM.getOrder_qed()) {
        case NLO_QED22:
        case NLO_QED12:
        case NLO_QED21:
        case NLO_QED02:            
        case NLO_QED11:
            mcM.setCoeff(0, aletilde * (1. / sW2 * (1.11 - 1.15 * (1. - Mt_muw * Mt_muw / 170. / 170.) - 0.444 * log(mH / 100.) -
                    0.21 * log(mH / 100.) * log(mH / 100.) - 0.513 * log(mH / 100.) * log(Mt_muw / 170.)) +
                    (8. / 9. * C7funLO(xt) - 104. / 243.) * L), NLO_QED11);
            mcM.setCoeff(1, aletilde * (1. / sW2 * (-0.143 + 0.156 * (1. - Mt_muw * Mt_muw / 170. / 170.) - 0.129 * log(mH / 100.) -
                    0.0244 * log(mH / 100.) * log(mH / 100.) - 0.037 * log(mH / 100.) * log(Mt_muw / 170.)) +
                    (4. / 9. * C8funLO(xt) - 4. / 3. * C7funLO(xt) - 58. / 81.) * L), NLO_QED11);
        case LO_QED:
            break;
        default:
            std::stringstream out;
            out << mcM.getOrder();
            throw "StandardModelMatching::mc_M(): order " + out.str() + "not implemented";

    }

    return (mcM);
}

double StandardModelMatching::fbb(double x) // from hep-ph/9707243 
{
    return(-0.0380386-2.24502*x+3.8472*x*x-7.80586*x*x*x+10.0763*x*x*x*x-6.9751*x*x*x*x*x
            +1.95163*x*x*x*x*x*x-0.00550335*log(x)); // fitted between 0 and 1
}

double StandardModelMatching::gbb(double x) // from hep-ph/9707243 
{
    if(x > 4.)
        return(sqrt(x-4.)*log((1.-sqrt(1.-4./x))/(1.+sqrt(1.-4./x))));
    else if(x >= 0.)
        return(2.*sqrt(4.-x)*acos(sqrt(x/4.)));
    else
        throw "StandardModelMatching::gbb(): defined for non-negative argument only.";
}

double StandardModelMatching::taub2(double x) // from hep-ph/9707243 
{
    double ll=log(x);
    return( 9.-13./4.*x-2.*x*x-x/4.*(19.+6.*x)*ll-x*x/4.*(7.-6.*x)*ll*ll-(1./4.+7./2.*x*x-3.*x*x*x)*M_PI*M_PI/6.+
            (x/2.-2.)*sqrt(x)*gbb(x)+(x-1.)*(x-1.)*(4.*x-7./4.)*gslpp_special_functions::dilog(1.-x)-(x*x*x-33./4.*x*x+
            18.*x-7.)*fbb(x) );
}

double StandardModelMatching::Delta_t(double mu, double x) // from hep-ph/9707243
{
    return(18.*log(mu/Mt_muw)+11.-x/2.+x*(x-6.)/2.*log(x)+(x-4.)/2.*sqrt(x)*gbb(x));
}

WilsonCoefficient& StandardModelMatching::mc_L()
{
//    double sW2 = (M_PI * SM.getAle() ) / ( sqrt(2.) * SM.getGF() * SM.Mw() * SM.Mw() ); // ******* FOR TEST *********
    double xt = x_t(Muw);
    double xht = SM.getMHl()*SM.getMHl()/Mt_muw/Mt_muw;

    
    switch (mcL.getScheme()) {
        case NDR:
            //case HV:
            //case LRI:
            break;
        default:
            std::stringstream out;
            out << mcL.getScheme();
            throw "StandardModel::mc_L(): scheme " + out.str() + "not implemented";
    }

    mcL.setMu(Muw); // cleared too
    
    switch (mcL.getOrder_qed()) {
        case NLO_QED12:
        case NLO_QED02:            
        case NLO_QED22:
            //Eqs. (32-33) of ref. Huber et al.
            //Delta_t and tau_b in hep-ph/9707243
//             mcL.setCoeff(0, aletilde*aletilde*(-xt*xt/32/sW2/sW2*(4.*sW2-1.)*(3.+taub2(xht)-Delta_t(Muw,xht))), NLO_QED22);
             mcL.setCoeff(1, aletilde*aletilde*(-xt*xt/32./sW2/sW2*(3.+taub2(xht)-Delta_t(Muw,xht))), NLO_QED22);
             mcL.setCoeff(0, (4. * sW2 - 1.) * (*(mcL.getCoeff(NLO_QED22)))(1), NLO_QED22);
//             mcL.setCoeff(1, aletilde*aletilde*Rest(xt, Muw)/sW2, NLO_QED22);
        case NLO_QED21:
            mcL.setCoeff(0, aletilde*alstilde*((1.-4.*sW2)/sW2*C1t(xt,Muw) - 1./sW2*B1t(xt,Muw) - D1t(xt,Muw) + 1./sW2 +
                                                524./729. - 128./243.*M_PI*M_PI - 16./3.*L - 128./81.*L*L), NLO_QED21);
            mcL.setCoeff(1, aletilde*alstilde*(1./sW2*(B1t(xt,Muw) - C1t(xt,Muw)) - 1./sW2), NLO_QED21);
        case NLO_QED11:
            mcL.setCoeff(0, aletilde*(1./sW2*Y0(xt) + Wt(xt) + 4./9. - 4./9.* 2.* log(Muw/Mt_muw)), NLO_QED11);
            mcL.setCoeff(1, aletilde*(-1./sW2*Y0(xt)), NLO_QED11);
        case LO_QED:
            break;
        default:
            std::stringstream out;
            out << mcL.getOrder_qed();
            throw "StandardModelMatching::mc_L(): order " + out.str() + "not implemented";
    }

    return (mcL);
}

WilsonCoefficient& StandardModelMatching::mc_Q()
{
    double xt = x_t(Muw);
    double xz = SM.getMz()*SM.getMz()/Mw/Mw;
    
    switch (mcQ.getScheme()) {
        case NDR:
            //case HV:
            //case LRI:
            break;
        default:
            std::stringstream out;
            out << mcQ.getScheme();
            throw "StandardModel::mc_Q(): scheme " + out.str() + "not implemented";
    }

    mcQ.setMu(Muw); // cleared too
    
    switch (mcQ.getOrder_qed()) {
        case NLO_QED22:
        case NLO_QED12:
        case NLO_QED02:            
        case NLO_QED21:
            mcQ.setCoeff(0, aletilde*alstilde*(4.*C1ew(xt) + 4.*D1t(xt,Muw) + 320./9.*C0b(xt) +
                                                1./sW2*(-2./3.*B1d(xt,Muw) - 2./9.*B1d_tilde(xt,Muw) +
                                                        2./3.*B1u(xt,Muw) + 2./9.*B1u_tilde(xt,Muw) +
                                                        4./3.*C1ew(xt) + 800./9.*B0b(xt) - 640./9.*C0b(xt))), NLO_QED21);
            mcQ.setCoeff(1, aletilde*alstilde*(-4./3.*Gew(xt,xz,Muw) - 16./3.*Hew(xt,xz,Muw) - 32.*C0b(xt) +
                                                1./sW2*(-4./3.*B1d_tilde(xt,Muw) + 4./3.*B1u_tilde(xt,Muw) +
                                                        4./3.*Gew(xt,xz,Muw) - 80.*B0b(xt) + 112./3.*C0b(xt))), NLO_QED21);
            mcQ.setCoeff(2, aletilde*alstilde*(-32./9.*C0b(xt) +
                                                1./sW2*(1./6.*B1d(xt,Muw) + 1./18.*B1d_tilde(xt,Muw) -
                                                        1./6.*B1u(xt,Muw) - 1./18.*B1u_tilde(xt,Muw) -
                                                        1./3.*C1ew(xt) - 80./9.*B0b(xt) + 64./9.*C0b(xt))), NLO_QED21);
            mcQ.setCoeff(3, aletilde*alstilde*(1./3.*Gew(xt,xz,Muw) + 1./3.*Hew(xt,xz,Muw) + 4.*C0b(xt) +
                                                1./sW2*(1./3.*B1d_tilde(xt,Muw) - 1./3.*B1u_tilde(xt,Muw) -
                                                        1./3.*Gew(xt,xz,Muw) +10.*B0b(xt) - 16./3.*C0b(xt))), NLO_QED21);
        case NLO_QED11:
            mcQ.setCoeff(0, aletilde*(4.*C0b(xt) + D0b_tilde(xt) + 4./9.*L - 1./sW2*(10./3.*B0b(xt)-4./3.*C0b(xt))), NLO_QED11); // log from Misiak's notes
            mcQ.setCoeff(2, aletilde*(1./sW2*(5./6.*B0b(xt)-1./3.*C0b(xt))), NLO_QED11);
        case LO_QED:
            break;
        default:
            std::stringstream out;
            out << mcQ.getOrder_qed();
            throw "StandardModelMatching::mc_Q(): order " + out.str() + "not implemented";
    }

    return (mcQ);
}

WilsonCoefficient& StandardModelMatching::mc_B()
{
    double xt = x_t(Muw);
    
    switch (mcB.getScheme()) {
        case NDR:
            //case HV:
            //case LRI:
            break;
        default:
            std::stringstream out;
            out << mcB.getScheme();
            throw "StandardModel::mc_B(): scheme " + out.str() + "not implemented";
    }

    mcB.setMu(Muw); // cleared too
    
    switch (mcB.getOrder_qed()) {
        case NLO_QED22:
        case NLO_QED12:
        case NLO_QED21:
        case NLO_QED02:            
        case NLO_QED11:
            mcB.setCoeff(0, aletilde*(-1./2./sW2*S0(xt)), NLO_QED11);
        case LO_QED:
            break;
        default:
            std::stringstream out;
            out << mcB.getOrder_qed();
            throw "StandardModelMatching::mc_B(): order " + out.str() + "not implemented";
    }

    return (mcB);
}

unsigned int StandardModelMatching::setCMDF1(WilsonCoefficient& CMDF1, WilsonCoefficient& DF1block,
        unsigned int tot, schemes scheme, orders order, orders_qed order_qed)
{
    unsigned int j, nops = DF1block.getSize();
    int ord;

    for (ord = LO; ord <= order; ord++)
        for (j = 0; j < nops; j++)
            CMDF1.setCoeff(j + tot, (*(DF1block.getCoeff((orders) ord)))(j), (orders) ord);

    for (ord = LO_QED; ord <= order_qed; ord++)
        for (j = 0; j < nops; j++)
            CMDF1.setCoeff(j + tot, (*(DF1block.getCoeff((orders_qed) ord)))(j), (orders_qed) ord);

    return (nops + tot);
}

std::vector<WilsonCoefficient>& StandardModelMatching::CMDF1(std::string blocks, unsigned int nops)
{
    unsigned int tot = 0;
    schemes scheme = mcC.getScheme();
    orders order = mcC.getOrder();
    orders_qed order_qed = mcC.getOrder_qed();
    WilsonCoefficient mcDF1(nops, scheme, order, order_qed);

    typedef WilsonCoefficient& (StandardModelMatching::*BlockM)();
    std::vector< std::pair<std::string, BlockM> > Methods = {
        std::make_pair("C", &StandardModelMatching::mc_C),
        std::make_pair("P", &StandardModelMatching::mc_P),
        std::make_pair("M", &StandardModelMatching::mc_M),
        std::make_pair("L", &StandardModelMatching::mc_L),
        std::make_pair("Q", &StandardModelMatching::mc_Q),
        std::make_pair("B", &StandardModelMatching::mc_B)
    };

    mcDF1.setMu(Muw);
    for (std::vector< std::pair<std::string, BlockM> >::iterator it = Methods.begin(); it != Methods.end(); it++)
        if (blocks.find(it->first) != std::string::npos)
            tot = setCMDF1(mcDF1, (this->*(it->second))(), tot, scheme, order, order_qed);

    vmcDF1.push_back(mcDF1);
    return (vmcDF1);
}
