/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "StandardModelMatching.h"
#include "StandardModel.h"
#include "QCD.h"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>
#include <stdexcept>


StandardModelMatching::StandardModelMatching(const StandardModel & SM_i) 
: ModelMatching(), SM(SM_i),
        mcdbd2(5, NDR, NLO),
        mcdbs2(5, NDR, NLO),
        mcdd2(5, NDR, NLO),
        mcdk2(5, NDR, NLO),
        mck(10, NDR, NLO),
        mckcc(10, NDR, NLO),
        mcbsg(10, NDR, NLO),
        mcbnlep(10, NDR, NLO, NLO_ew),
        mcbnlepCC(10, NDR, NLO),
        mcd1(10, NDR, NLO),
        mcd1Buras(10, NDR, NLO),
        mckpnn(1, NDR, NLO, NLO_ew),
        mckppnn(1, NDR, NLO, NLO_ew),
        mckmm(1, NDR, NLO),
        mcbsnn(1, NDR, NLO),
        mcbdnn(1, NDR, NLO),
        mcbsmm(1, NDR, NLO),
        mcbdmm(1, NDR, NLO),
        Vckm(3, 3, 0)
{
    swa = 0.;
    swb = 0.;
    swc = 0.;
    xcachea = 0.;
    xcacheb = 0.;
    xcachec = 0.;
    

    for (int j=0; j<10; j++) {
        CWbsgArrayLO[j] = 0.;
        CWbsgArrayNLO[j] = 0.;
        CWD1ArrayLO[j] = 0.; 
        CWD1ArrayNLO[j] = 0.;
        CWbnlepArrayLOqcd[j] = 0.;
        CWbnlepArrayNLOqcd[j] = 0.;
        CWbnlepArrayLOew[j] = 0.;
        CWbnlepArrayNLOew[j] = 0.;
    };
    
}

void StandardModelMatching::updateSMParameters()
{
    Mut = SM.getMut();
    Muw = SM.getMuw();
    Ale = SM.getAle();
    GF = SM.getGF();
    MW_tree = SM.Mw_tree();
    Vckm = SM.getVCKM();
    lam_t = SM.getlamt();
//   std::cout << "UPDATED SM" << std::endl;

}

double StandardModelMatching::x_c(const double mu, const orders order) const 
{
    double mc = SM.Mrun(mu, SM.getQuarks(QCD::CHARM).getMass_scale(), 
                        SM.getQuarks(QCD::CHARM).getMass(), order);
    return pow(mc / MW_tree, 2.); 
}

double StandardModelMatching::x_t(const double mu, const orders order) const 
{
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), order);
    return pow(mt / MW_tree, 2.); 
}

double StandardModelMatching::mt2omh2(const double mu, const orders order) const 
{
    double mt = SM.Mrun(mu, SM.getQuarks(QCD::TOP).getMass_scale(), 
                        SM.getQuarks(QCD::TOP).getMass(), order);
    return pow(mt / SM.getMHl(), 2.);
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
    return (x * (-4. + 18. * x + 3. * x * x + x * x * x) / 4. / pow(x - 1., 3.)
            - 9. * x * x * x / 2. / pow(x - 1., 4.) * log(x));
}

double StandardModelMatching::S11(double x) const
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm3 = pow(x - 1., 3.);
    double xm4 = xm3 * (x - 1);
    
    return (x * (4. - 39. * x + 168. * x2 + 11. * x3) / 4. / xm3
            + 3. * x3 * gsl_sf_dilog(1. - x) * (5. + x) / xm3
            + 3. * x * log(x)*(-4. + 24. * x - 36. * x2 - 7. * x3 - x4) / 2.
            / xm4 + 3. * x3 * pow(log(x), 2.) * (13. + 4. * x + x2) / 2.
            / xm4);
}

double StandardModelMatching::S18(double x) const 
{
    double x2 = x * x;
    double x3 = x2 * x;
    double x4 = x3 * x;
    double xm2 = pow(x - 1., 2.);
    double xm3 = xm2 * (x - 1);
    
    return ((-64. + 68. * x + 17. * x2 - 11. * x3) / 4. / xm2
            + pow(M_PI, 2.) * 8./3. / x + 2. * gsl_sf_dilog(1. - x) * (8. - 24. * x
            + 20. * x2 - x3 + 7. * x4 - x4 * x) / (x * xm3)
            + log(x)*(-32. + 68. * x - 32. * x2 + 28. * x3 - 3. * x4)
            / (2. * xm3) + x2 * pow(log(x), 2.) * (4. - 7. * x + 7. * x2
            - 2. * x3) / (2. * x4));
}

double StandardModelMatching::S1(double x) const 
{
    return (SM.getCF() * S11(x) + (SM.getNc() - 1.) / 2. / SM.getNc() * S18(x));
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
    return ( (-3. * x * x * x + 2. * x * x)/(2. * pow(1. - x, 4.)) * log(x) + (22. * x * x * x - 153. * x * x +
            159. * x - 46.)/(36. * pow(1. - x, 3.)) );
}

double StandardModelMatching::B0t(double x) const
{
    return( x / (4.* (1. - x) * (1. - x)) * log(x) + 1. / 4. * (1. - x) );
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
    
    return (x * (18. - 11. * x - x * x) / (12. * pow(1. - x, 3.) + x * x * (15. - 16. * x + 4. * x * x) /
            (6. * pow(1. - x, 4.)) * log(x) - 2./3. * log(x)));
}

double StandardModelMatching::F0t(double x) const
{
    double x2 = x * x;
    double xm3 = pow(1 - x, 3);
    
    return ( (3. * x2) / (2. * xm3 * (1 - x)) * log(x) + ( 5. * x2 * x - 9. * x2 + 30. * x - 8)/
            (12. * xm3) );
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
    
    double Li2 = gsl_sf_dilog(1.-1./x);
    return( Li2 * ( -16. * x4 - 122. * x3 + 80. * x2 - 8. * x) / (9. * xm4) +
            (6. * x4 + 46. * x3 -28. * x2) / (3. * xm5) * logx * logx +
            (-102. * x4 * x - 588. * x4 + 3244. * x2 - 1364. * x + 208.) / (81. * xm5) * logx +
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
    
    double Li2 = gsl_sf_dilog(1.-1./x);
    return(Li2 * ( -4. * x4 + 40. * x3 + 41. * x2 + x) / (6. * xm4) +
            (-17. * x3 - 31. * x2) / (2. * xm5) * logx * logx +
            (-210. * x * x4 + 1086. * x4 + 4893. * x3 + 2857. * x2 - 1994. * x + 280.)/(216. * xm5) * logx+
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

double StandardModelMatching::E0b(double x) const
{
    double x2 = x * x;
    
    return ( -2./3. * log(x) + (18. * x - 11. * x2 - x2 * x) / (12.* (- x2 * x +3. * x2 -3. * x +1.)) +
            (x2 * (15. - 16. * x +4. * x2))/(6.*pow(1.-x,4.)) *log(x) );
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
            - gsl_sf_dilog(1.-x) * (4. * x - x3) / ((1. - x) * (1. - x))
            - 8. * x * log(Mut*Mut/Muw/Muw) * (8. - 9. * x + x3 + 6. * logx)/8./xm3 );
  }

double StandardModelMatching::Xewt(double x, double a, double mu) const{
    double b = 0.;
    // WARNING: sin^{2}(theta_w) needs to be implemented at NLO in the MSbar scheme, according to hep-ph/1009.0947v2
    double swsq = (M_PI * Ale )/( sqrt(2) * GF * MW_tree * MW_tree);
    
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
    
    C[10] = 2. * log(mu/MW_tree) / xm2;
    
    C[11] = logx * 2. * log(mu/MW_tree) / xm3;
    
    C[12] = loga / (xm2 * axm);
    
    C[13] = logx * loga / (2. * a * xm3 * x * axm);
    
    C[14] = gsl_sf_dilog(1. - a) / xm2;
    
    C[15] = gsl_sf_dilog(1. - x) / a / x;
    
    C[16] = gsl_sf_dilog(1. - a * x) / (a * x * xm2); 
    
    for (int i=0; i<10; i++){
        b += C[i]*A[i];
    }
    
    return (b/128./swsq);
}

double StandardModelMatching::phi1(double z) const{
    if (z >= 0.) {
        if (z < 1){
                return(4. * sqrt(z / (1. - z)) * gsl_sf_clausen(2. * asin(sqrt(z))));
        }
        else{
            return((1. / sqrt(1. - 1. / z)) * (2. * log(0.5 * sqrt(1. - 1. / z)) * log(0.5 * sqrt(1. -1. / z))- 4. * gsl_sf_dilog(0.5 * (1. - sqrt(z / (1. - z)))) - log(4. * z) * log(4. * z)
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
        return( 1. / l * (M_PI * M_PI/3. + 2. * log(0.5 * (1. + x - y - l)) * log(0.5 * (1. - x + y - l)) - log(x) * log(y) - 2. * gsl_sf_dilog(0.5 * (1. + x - y - l)) - 2. * gsl_sf_dilog(0.5 * (1. -x + y - l))));
    }
    else if((l * l) < 0. || (sqrt(x) + sqrt(y)) > 1.){
        return(2./( -l * l) * (gsl_sf_clausen(2. * acos((-1. + x + y)/(2. * sqrt(x * y))))
                +gsl_sf_clausen(2. * acos((1. + x - y) / (2. * sqrt(x))))
                +gsl_sf_clausen(2. * acos((1. - x + y) / (2. * sqrt(y))))));
    }
    else{
        std::stringstream out;
        out << x;
        throw std::runtime_error("StandardModelMatching::phi2(double x, double y) wrong" + out.str());
    }
    return(0.);
}

double StandardModelMatching::Y0t(double x) const{
    return(x / 8. * ((4. - x)/(1. - x) + 3. * x * log(x) / (1. - x) / (1. - x)) );
}

double StandardModelMatching::Y1t(double x) const{
    double x2 = x * x;
    double x4 = x2 * x2;
    double logx = log(x);
    double logMutMuW = log(Mut * Mut / Muw / Muw);
    
    return( (4. * x + 16. * x2 + 4. * x2 * x) / (3. * (1. - x) * (1. - x))
           -(4. * x - 10. * x2 - x2 * x - x4) * logx / (1.-x) / (1.-x) / (1.-x)
           +(2. * x - 14. * x4 * x - x4) * logx * logx / 2. / (1.-x) / (1.-x) / (1.-x)
           +(2. * x4) / (1.-x) / (1.-x) * gsl_sf_dilog(1.-x)
           +x * ((4.-x) / (1.-x) + 3. * x * logx * (1.-x) * (1.-x)) * logMutMuW
           +x2 * (-1. / (1.-x) + (4.-x) / (1.-x) / (1.-x) + 3. * logx / (1.-x) / (1.-x)
                 + 6. * x * logx / (1.-x) / (1.-x) / (1.-x) + 3. / (1.-x) / (1.-x)) * logMutMuW
            );
}

/******************************************************************************/
const std::vector<WilsonCoefficient>& StandardModelMatching::CMdbd2() 
{   
//    if(SM_i == SM)
//        return(vmc);
    
    double gammam = 8.;                                                         
    double Bt;  
    
    
    double xt = x_t(Muw);
    complex co = GF / 4. / M_PI * MW_tree * SM.getlamt_d();
    double Nc = SM.getNc();

    vmcdb.clear();

    switch (mcdbd2.getScheme()) {
        case NDR:
            Bt = 5. * (Nc - 1.) / 2. / Nc + 3. * SM.getCF();
            break;
        case HV:
        case LRI:
        default:
            std::stringstream out;
            out << mcdbd2.getScheme();
            throw std::runtime_error("StandardModel::CMdb2(): scheme " + out.str() + "not implemented"); 
    }

    mcdbd2.setMu(Muw);
    
    switch (mcdbd2.getOrder()) {
        case NNLO:
        case NLO:
            mcdbd2.setCoeff(0, co * co * 4. * (SM.Als(Muw, FULLNLO) / 4. / M_PI * (S1(xt) +
                    Bt * S0(xt, xt) + 2. * gammam * S0p(xt) * log(Muw / MW_tree))), NLO);   
        case LO:
            mcdbd2.setCoeff(0, co * co * 4. * S0(xt, xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbd2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbd2(): order " + out.str() + "not implemented"); 
    }
    

    vmcdb.push_back(mcdbd2);
    return(vmcdb);
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMdbs2() 
{   
//    if(SM_i == SM)
//        return(vmc);
    
   
    double gammam = 8.;
    double Bt;
    double xt = x_t(Muw);
    complex co = GF / 4. / M_PI * MW_tree * SM.getlamt_s();
    double Nc = SM.getNc();

    vmcds.clear();

    switch (mcdbs2.getScheme()) {
        case NDR:
            Bt = 5. * (Nc - 1.) / 2. / Nc + 3. * SM.getCF();
            break;
        case HV:
        case LRI:
        default:
            std::stringstream out;
            out << mcdbs2.getScheme();
            throw std::runtime_error("StandardModel::CMdbs2(): scheme " + out.str() + "not implemented"); 
    }

    mcdbs2.setMu(Muw);
 
    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:          
            mcdbs2.setCoeff(0, co * co * 4. * (SM.Als(Muw, FULLNLO) / 4. / M_PI * (S1(xt) +
                    Bt * S0(xt, xt) + 2. * gammam * S0p(xt) * log(Muw / MW_tree))), NLO);      
         case LO:
            mcdbs2.setCoeff(0, co * co * 4. * S0(xt, xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcdbs2.getOrder();
            throw std::runtime_error("StandardModelMatching::CMdbs2(): order " + out.str() + "not implemented"); 
    }    

    vmcds.push_back(mcdbs2);
    return(vmcds);
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMdk2() 
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMd1Buras()
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
            mcd1Buras.setCoeff(0, SM.Als(Muw, FULLNLO) / 4. / M_PI  * 11./2. , NLO);
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMd1()
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
            mcd1.setCoeff(0, SM.Als(Muw, FULLNLO) / 4. / M_PI  * 15. , NLO);
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMdd2() 
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMK(){
    
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
                mck.setCoeff(j, lam_t * SM.Als(Muw, FULLNLO) / 4. / M_PI * 
                                setWCbnlep(j, xt,  NLO), NLO);
                mck.setCoeff(j, lam_t * Ale / 4. / M_PI *
                                setWCbnlepEW(j, xt), NLO_ew);
                }
        case LO:
            for (int j=0; j<10; j++){
                mck.setCoeff(j, lam_t *  setWCbnlep(j, xt,  LO), LO);
                mck.setCoeff(j, 0., LO_ew); 
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
const std::vector<WilsonCoefficient>& StandardModelMatching::CMKCC(){
    
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
 * Wilson coefficients misiak base for b -> s gamma                            * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *         
 *                 - semileptonic                                              * 
 * ****************************************************************************/
const std::vector<WilsonCoefficient>& StandardModelMatching::CMbsg() 
{    
    double xt = x_t(Muw);
    complex co = (- 4. * GF / sqrt(2)) * SM.getlamt_s();
    
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
    
    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<10; j++){
            mcbsg.setCoeff(j, co * SM.Als(Muw, FULLNLO) / 4. / M_PI * setWCbsg(j, xt,  NLO) , NLO);
            }
            std::cout<<std::endl;
        case LO:
            for (int j=0; j<10; j++){
            mcbsg.setCoeff(j, co * setWCbsg(j, xt,  LO), LO);
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


/*******************************************************************************
 * Wilson coefficients calcoulus, misiak base for b -> s gamma                  *  
 * ****************************************************************************/

double StandardModelMatching::setWCbsg(int i, double x, orders order)
{    
    sw =  sqrt( (M_PI * Ale )/( sqrt(2) * GF * MW_tree * MW_tree) ) ;

    if ( swa == sw && xcachea == x){
        switch (order){
        case NNLO:
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
    
    swa = sw; xcachea = x;
    // this function returns the effective Wilson coefficients if  CWbsgArrayNLO[7] = C8NLOeff(x);
    //                                                             CWbsgArrayNLO[6] = C7NLOeff(x);
    //                                                             CWbsgArrayLO[6] = C7LOeff(x);
    //                                                             CWbsgArrayLO[7] = C8LOeff(x);
    // or the standard one if CWbsgArrayNLO[7] = -0.5 * A0t(x)- 23./36.;
    //                        CWbsgArrayNLO[6] = -0.5 * F0t(x)- 1./3.;
    //                        CWbsgArrayLO[6] = 0.;
    //                        CWbsgArrayLO[7] = 0.;
    switch (order){
        case NNLO:
        case NLO:
            CWbsgArrayNLO[0] = 15.;
            CWbsgArrayNLO[3] = E0t(x)-(2./3.);
            CWbsgArrayNLO[6] = C7NLOeff(x);//-0.5 * A0t(x)- 23./36.;
            CWbsgArrayNLO[7] = C8NLOeff(x);//-0.5 * F0t(x)- 1./3.;
            CWbsgArrayNLO[8] = (1-4.*sw*sw) / sw *C0t(x) - 1./(sw*sw) *
                                B0t(x) - D0t(x) + 38./27. + 1./(4.*sw*sw);
            CWbsgArrayNLO[9] = 1./(sw*sw) * (B0t(x) - C0t(x)) -1./(4.*sw*sw);
        case LO:
            CWbsgArrayLO[1] = 1.;
            CWbsgArrayLO[6] = C7LOeff(x);//0.;
            CWbsgArrayLO[7] = C8LOeff(x);//0.;
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("order" + out.str() + "not implemeted"); 
            }
    
    switch (order){
        case NNLO:
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

/******************************************************************************/

/*******************************************************************************
 * Wilson coefficients Buras base for b -> nonlep. decays                      * 
 * operator basis: - current current                                           *         
 *                 - qcd penguins                                              * 
 *                 - magnetic and chromomagnetic penguins                      *         
 *                 - semileptonic                                              *
 * i=0 deltaS=0 deltaC=0;  i=1 1,0 ;                                           *
 * ****************************************************************************/
const std::vector<WilsonCoefficient>& StandardModelMatching::CMbnlep(const int& a)
{
    complex lambda;
    
    switch (a) {
        case 0: lambda = SM.getlamt_d();
        break;
        case 1: lambda = SM.getlamt_s();
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
                mcbnlep.setCoeff(j,co * lambda * SM.Als(Muw, FULLNLO) / 4. / M_PI * 
                                setWCbnlep(j, xt,  NLO), NLO);
                mcbnlep.setCoeff(j, co * lambda * Ale / 4. / M_PI *
                                setWCbnlepEW(j, xt), NLO_ew);
                }
        case LO:
            for (int j=0; j<10; j++){
                mcbnlep.setCoeff(j, co * lambda *  setWCbnlep(j, xt,  LO), LO);
                mcbnlep.setCoeff(j, 0., LO_ew); 
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
const std::vector<WilsonCoefficient>& StandardModelMatching::CMbnlepCC(const int& a) 
{    
    complex lambda1 = 0.;
    //complex lambda2 = 0.;
    //matrix<complex> ckm = SM.getVCKM();
    
    switch (a) {
        case 0: lambda1 = SM.getlamu_d();
                break;
        case 1: lambda1 = SM.getlamu_s();
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
            throw std::runtime_error("case" + out.str() + "unexsting; implemented i=0,1,2,3"); 
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMkp0nn() {
    
    //PROBLEMI: mu e sin(theta_weak) e la scala di als
    
    double xt = x_t(Muw);
    double a = mt2omh2(Muw);
    
    vmckpnn.clear();
    
    mckpnn.setMu(Mut);
 
    switch (mckpnn.getOrder()) {
        case NNLO:
        case NLO:
            mckpnn.setCoeff(0, SM.Als(Muw, FULLNLO)/4./M_PI*lam_t.imag()*X1t(xt)/SM.GetLambda(), NLO);
        case LO:
            mckpnn.setCoeff(0, lam_t.imag()*X0t(xt)/SM.GetLambda(), LO);
            break;
        default:
            std::stringstream out;
            out << mckpnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkp0nn(): order " + out.str() + "not implemented"); 
    }
    
    switch (mckpnn.getOrder_ew()) {
        case NLO_ew:
            mckpnn.setCoeff(0, Ale/4./M_PI*lam_t.imag()*Xewt(xt, a, Muw)/SM.GetLambda(), NLO_ew);
        case LO_ew:
            break; 
        default:
            std::stringstream out;
            out << mckpnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkp0nn(): order " + out.str() + "not implemented"); 
    }

    vmckpnn.push_back(mckpnn);
    return(vmckpnn);
    
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMkppnn() {
    
    //PROBLEMI: mu e sin(theta_weak) e la scala di als
    
    double xt = x_t(Muw);
    double a = mt2omh2(Muw);
    
    vmckppnn.clear();
    
    mckppnn.setMu(Mut);
 
    switch (mckppnn.getOrder()) {
        case NNLO:
        case NLO:
            mckppnn.setCoeff(0, SM.Als(Muw, FULLNLO)/4./M_PI*lam_t.imag()*X1t(xt)/SM.GetLambda(), NLO);
        case LO:
            mckppnn.setCoeff(0, lam_t.imag()*X0t(xt)/SM.GetLambda(), LO);
            break;
        default:
            std::stringstream out;
            out << mckppnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkppnn(): order " + out.str() + "not implemented"); 
    }
    
    switch (mckppnn.getOrder_ew()) {
        case NLO_ew:
            mckppnn.setCoeff(0, Ale/4./M_PI*lam_t.imag()*Xewt(xt, a, Muw)/SM.GetLambda(), NLO_ew);
        case LO_ew:
            break; 
        default:
            std::stringstream out;
            out << mckppnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkppnn(): order " + out.str() + "not implemented"); 
    }

    vmckppnn.push_back(mckppnn);
    return(vmckppnn);
    
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMkmm() {
    
    //PROBLEMI: mu e sin(theta_weak) e la scala di als
    
    double xt = x_t(Muw);
    
    vmckmm.clear();
    
    mckmm.setMu(Mut);
 
    switch (mckmm.getOrder()) {
        case NNLO:
        case NLO:
            mckmm.setCoeff(0, SM.Als(Muw, FULLNLO)/4./M_PI*lam_t.real()*Y1t(xt)/SM.GetLambda(), NLO);
        case LO:
            mckmm.setCoeff(0, lam_t.real()*Y0t(xt)/SM.GetLambda(), LO);
            break;
        default:
            std::stringstream out;
            out << mckmm.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkmm(): order " + out.str() + "not implemented"); 
    }

    vmckmm.push_back(mckmm);
    return(vmckmm);
    
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMBXsnn() {
    
    double xt = x_t(Muw);
    
    vmcbsnn.clear();
    
    mcbsnn.setMu(Mut);
 
    switch (mcbsnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (Vckm(2,1).abs() / Vckm(1,2).abs())*
                                SM.Als(Muw, FULLNLO)/4./M_PI*X1t(xt), NLO);
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMBXdnn() {
    
    double xt = x_t(Muw);
    
    vmcbdnn.clear();
    
    mcbdnn.setMu(Mut);
 
    switch (mcbdnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (Vckm(2,2).abs() / Vckm(1,2).abs()) *
                                SM.Als(Muw, FULLNLO)/4./M_PI*X1t(xt), NLO);
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

const std::vector<WilsonCoefficient>& StandardModelMatching::CMbsmm() {
    
    double xt = x_t(Muw);
    
    vmcbsmm.clear();
    
    mcbsmm.setMu(Mut);
 
    switch (mcbsmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (Vckm(2,2).conjugate() * Vckm(2,1)).abs() *
                                SM.Als(Muw, FULLNLO) / 4. / M_PI*Y1t(xt), NLO);
        case LO:
            mcbsnn.setCoeff(0,  (Vckm(2,2).conjugate() * Vckm(2,1)).abs() *
                                Y0t(xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcbsmm.getOrder();
            throw std::runtime_error("StandardModelMatching::CXdnn(): order " + out.str() + "not implemented"); 
    }
    
    vmcbsmm.push_back(mcbsmm);
    return(vmcbsmm);
    
}

const std::vector<WilsonCoefficient>& StandardModelMatching::CMbdmm() {
    
    double xt = x_t(Muw);
    
    vmcbdmm.clear();
    
    mcbdmm.setMu(Mut);
 
    switch (mcbdmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (Vckm(2,2).conjugate() * Vckm(2,0)).abs() *
                                SM.Als(Muw, FULLNLO) / 4. /M_PI * Y1t(xt), NLO);
        case LO:
            mcbsnn.setCoeff(0,  (Vckm(2,2).conjugate() * Vckm(2,0)).abs() *
                                Y0t(xt), LO);
            break;
        default:
            std::stringstream out;
            out << mcbdmm.getOrder();
            throw std::runtime_error("StandardModelMatching::CXdnn(): order " + out.str() + "not implemented"); 
    }
    
    vmcbdmm.push_back(mcbdmm);
    return(vmcbdmm);
    
}


/*******************************************************************************
 * Wilson coefficients calculus, Buras base for nonlep. b decays               *  
 * ****************************************************************************/
double StandardModelMatching::setWCbnlep(int i, double x, orders order) 
{    
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2) * GF * MW_tree * MW_tree) );
    
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
    sw =  sqrt( (M_PI * Ale ) / ( sqrt(2) * GF * MW_tree * MW_tree) ) ;
    
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

complex StandardModelMatching::S0c() const 
{
    double xc = x_c(SM.getMuc());
    complex co = GF / 2. / M_PI * MW_tree * SM.getlamc().conjugate();
    
    return(co * co * S0(xc, xc));
}

complex StandardModelMatching::S0ct() const 
{
    double xc = x_c(SM.getMuc());
    double xt = x_t(Mut);
    double co = GF / 2. / M_PI * MW_tree;
    
    return( co * co * 2. * SM.getlamc().conjugate() * lam_t.conjugate() * S0(xc, xt) );
}

complex StandardModelMatching::S0tt() const
{
    double xt = x_t(Mut);
    complex co = GF / 2. / M_PI * MW_tree * lam_t.conjugate();
    
    return ( co * co * S0(xt, xt) );
}