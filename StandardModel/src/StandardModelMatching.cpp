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
        mcdbd2(5, NDR, NLO), mcdbs2(5, NDR, NLO),
        mcdd2(5, NDR, NLO), mcdk2(5, NDR, NLO),
        mck(10, NDR, NLO), mckcc(10, NDR, NLO), 
        mcbsg(10, NDR, NLO), mcbnlep(10, NDR, NLO, NLO_ew), 
        mcbnlepCC(10, NDR, NLO), mcd1(10, NDR, NLO), 
        mcd1Buras(10, NDR, NLO) , mckpnn(1, NDR, NLO, NLO_ew), mckppnn(1, NDR, NLO, NLO_ew),
        mckmm(1, NDR, NLO), mcbsnn(1, NDR, NLO), mcbdnn(1, NDR, NLO), 
        mcbsmm(1, NDR, NLO), mcbdmm(1, NDR, NLO) {   
    
    swa = 0.;  swb = 0.; swc = 0.;
    xcachea = 0.; xcacheb = 0.; xcachec = 0.;
    
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

double StandardModelMatching::S0(double x) const 
{
    return S0(x, x);
}

double StandardModelMatching::S0(double x, double y) const 
{ // Buras 2000 Appendix
    if (fabs(1. - y / x) < LEPS){
        return ((x * (-4. + 15. * x - 12. * x * x + x*x*x +
            6. * x * x * log(x))) / (4. * pow(-1. + x, 3.)));
    }
    else
        return (x * y * ((1. / 4. + 3. / 2. / (1. - x) - 3. / 4. / pow(1. - x, 2.)) *
            log(x) / (x - y) +
            (1. / 4. + 3. / 2. / (1. - y) - 3. / 4. / pow(1. - y, 2.)) *
            log(y) / (y - x) -
            3. / 4. / (1. - x) / (1. - y)));
}

double StandardModelMatching::S0p( double x) const 
{
    return (x * (-4. + 18. * x + 3. * x*x + x*x*x) / 4. / pow(x - 1., 3.)
            - 9. * x*x*x / 2. / pow(x - 1., 4.) * log(x));
}

double StandardModelMatching::S11(double x) const 
{
    return (x * (4. - 39. * x + 168. * x*x + 11. * x*x*x) / 4. / pow(x - 1., 3.)
            + 3. * x*x*x * gsl_sf_dilog(1. - x)*(5. + x) / pow(x - 1., 3.)
            + 3. * x * log(x)*(-4. + 24. * x - 36. * x*x - 7. * x*x*x - x*x*x*x) / 2.
            / pow(x - 1., 4.) + 3. * x*x*x * pow(log(x), 2.)*(13. + 4. * x + x*x) / 2.
            / pow(x - 1., 4.));
}

double StandardModelMatching::S18(double x) const 
{
    return ((-64. + 68. * x + 17. * x*x - 11. * x*x*x) / 4. / pow(x - 1., 2.)
            + pow(M_PI, 2.)*8. / 3. / x + 2. * gsl_sf_dilog(1. - x)*(8. - 24. * x
            + 20. * x*x - x*x*x + 7. * x*x*x*x - pow(x, 5.)) / x / pow(x - 1., 3.)
            + log(x)*(-32. + 68. * x - 32. * x*x + 28. * x*x*x - 3. * x*x*x*x)
            / 2. / pow(x - 1., 3.) + x*x * pow(log(x), 2.)*(4. - 7. * x + 7. * x*x
            - 2. * x*x*x) / 2. / pow(x - 1., 4.));
}

double StandardModelMatching::S1(double x) const 
{
    return (SM.getCF() * S11(x) + (SM.getNc()-1.)/2./SM.getNc() * S18(x));
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
    return ( (-3.*x*x*x + 2.*x*x)/(2.*pow(1.-x,4.)) *log(x) + (22.*x*x*x - 153.*x*x +
            159.*x - 46.)/(36.*pow(1.-x,3.)) );
}

double StandardModelMatching::B0t(double x) const
{
    return( x / (4.* (1.-x)*(1.-x)) *log(x) + 1. / 4.*(1.-x) );
}

double StandardModelMatching::C0t(double x) const
{
    return( (3.*x*x + 2.*x)/(8.*(1.-x)*(1.-x)) *log(x) + (-x*x + 6.*x)/(8.*(1.-x)) );
}

double StandardModelMatching::D0t(double x) const
{
    return( (-3.*x*x*x*x + 30.*x*x*x - 54.*x*x + 32.*x -8.) / (18.*pow(1. -x, 4)) *log(x)
            + (-47.*x*x*x + 237.*x*x -312.*x +104.) / (108.*pow (1.-x,3.)) );
}

double StandardModelMatching::E0t(double x) const
{
    return (x*(18. - 11.*x - x*x)/(12. *pow(1.-x,3.) + x*x*(15. - 16.*x + 4.*x*x)/
            (6. *pow(1.-x,4.))*log(x) - 2./3.*log(x)));
}

double StandardModelMatching::C7LOeff(double x) const
{
    return( (3.*x*x*x - 2.*x*x)/(4.*pow(x-1.,4.))*log(x) + (-8.*x*x*x - 5.*x*x +
              7.*x)/(24.*pow(x-1.,3.)) );
}

double StandardModelMatching::C8LOeff(double x) const
{
    return( -3.*x*x/(4.*pow(x-1.,4.))*log(x) + (-x*x*x + 5.*x*x + 2.*x)/(8.*pow(x-1.,3)) );
}

double StandardModelMatching::C7NLOeff(double x) const
{
    double Li2 = gsl_sf_dilog(1.-1./x);
    return( Li2*( -16.*x*x*x*x - 122.*x*x*x + 80.*x*x - 8.*x)/(9.*pow(x-1.,4.)) +
            (6.*x*x*x*x + 46.*x*x*x -28.*x*x)/(3.*pow(x-1.,5.))*log(x)*log(x) + 
            (-102.*x*x*x*x*x - 588.*x*x*x*x + 3244.*x*x - 1364.*x + 208.)/(81.*pow(x-1.,5.))*log(x) + 
            (1646.*x*x*x*x + 12205.*x*x*x - 10740.*x*x + 2509.*x - 436.)/(486.*pow(x-1.,4.)));
}

double StandardModelMatching::C8NLOeff(double x) const
{
    double Li2 = gsl_sf_dilog(1.-1./x);
    return(Li2*(-4.*x*x*x*x + 40.*x*x*x + 41.*x*x + x)/(6.*pow(x-1.,4.)) +
            (-17.*x*x*x - 31.*x*x)/(2.*pow(x-1.,5.))*log(x)*log(x) +
            (-210.*x*x*x*x*x + 1086.*x*x*x*x + 4893.*x*x*x + 2857.*x*x - 1994.*x + 280.)/(216.*pow(x-1.,5.))*log(x)+
            (737.*x*x*x*x - 14102.*x*x*x - 28209.*x*x + 610.*x -508.)/(1296.*pow(x-1.,4.)));
}

double StandardModelMatching::F0t(double x) const
{
    return ( (3.*x*x)/(2.*pow((1.-x), 4.)) * log(x) + ( 5.*x*x*x - 9.*x*x + 30.*x -8)/
            (12.*(1.-x)*(1.-x)*(1.-x)) );
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
    return ( 0.25*( x/(1.-x) + x/(x*x-2.*x+1.) * log(x) ) );
}

double StandardModelMatching::C0b(double x) const
{
    return ( x/8. * ( (x-6.)/(x-1.) + (3.*x+2.)/(x*x-2.*x+1.) * log(x) ) );
}

double StandardModelMatching::D0b(double x) const
{
    return ( -4./9. * log(x) + (-19.*x*x*x + 25.*x*x)/(36.*(x*x*x -3.*x*x +3.*x -1.)) 
            + ( x*x*(5.*x*x -2.*x -6.) )/(18.*pow(x-1.,4.)) * log(x) );
}

double StandardModelMatching::E0b(double x) const
{
    return ( -2./3. *log(x) + (x*(18. -11.*x -x*x))/(12.*(-x*x*x +3.*x*x -3.*x +1.)) +
            (x*x*(15. -16.*x +4.*x*x))/(6.*pow(1.-x,4.)) *log(x) );
}

/******************************************************************************/
/* loop functions for rare K and B decays, K-> pi nu nu & B-> Xs nu nu        */
/******************************************************************************/

double StandardModelMatching::X0t(double x) const{
    return ( (x/8.)*( (x+2.)/(x-1.) - (6.-3.*x)/(1.-x)/(1.-x)*log(x) ) );
}

double StandardModelMatching::X1t(double x) const{
    return ( 
            -(29.*x -x*x -4.*x*x*x)/(3.*(1.-x)*(1.-x))
            -log(x)*(x +9.*x*x -x*x*x -x*x*x*x)/((1.-x)*(1.-x)*(1.-x)) 
            +log(x)*log(x)*(8.*x +4.*x*x +x*x*x -x*x*x*x)/(2.*(1.-x)*(1.-x)*(1.-x))
            -gsl_sf_dilog(1.-x)*(4.*x -x*x*x)/((1.-x)*(1.-x))
            +8.*x* (1/8.)*( (x+2.)/(x-1.) + (3.*x-6.)/(x-1.)/(x-1.)*log(x) )
            *log(SM.getMut()*SM.getMut()/SM.getMuw()/SM.getMuw())
            +8.*x* (x/8.)*( 1./(x-1.) -(x+2.)/(x*x -2.*x +1.) + (3.*x -6.)/(x*(x*x -2.*x +1.))
            +3.*log(x)/(x*x -2.*x +1.) -2.*log(x)*(3.*x-6.)/(x-1.)/(x-1.)/(x-1.) ) );
}

double StandardModelMatching::Xewt(double x, double a, double mu) const{
    double b = 0.;
    
    double swsq = (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * SM.Mw_tree() * SM.Mw_tree());
    
    double A[17], C[17];
    
    A[0]=(16. - 48.*a)*M_PI*M_PI + (288.*a - (32. - 88.*a)*M_PI*M_PI )*x 
        + (2003.*a + 4.*(4. - 6.*a - a*a )*M_PI*M_PI )*x*x
        + (9.*a*(93. + 28.*a) - 4.*a*(3. - 2.*a + 8.*a*a )*M_PI*M_PI )*x*x*x
        + (3.*a*(172. - 49.*a - 32.*a*a ) + 4.*a*(20. - a + 16.*a*a )*M_PI*M_PI )*x*x*x*x
        - (3.*a*(168. + 11.*a - 24.*a*a ) + 4.*a*(45. + 8.*a*a )*M_PI*M_PI)*x*x*x*x*x
        + 96.*a*M_PI*M_PI*x*x*x*x*x*x ;
    
    A[1]=- 768.*x - (525. - 867.*a)*x*x + (303. + 318.*a)*x*x*x - 195.*a*x*x*x*x;
    
    A[2]=- 8.*(95. - 67.*a + 11.*a*a )*x*x + 2.*(662. - 78.*a - 177.*a*a + 40.*a*a*a )*x*x*x
        - (608. + 476.*a - 595.*a*a + 114.*a*a*a )*x*x*x*x 
        + (44. + 188.*a - 321.*a*a + 103.*a*a*a - 8.*a*a*a*a )*x*x*x*x*x
        - a*(28. - 72.*a + 33.*a*a - 4.*a*a*a )*x*x*x*x*x*x;
    
    A[3]=+ 48. - 10.*(57. + 4.*a)*x + 51.*(29. + 10.*a)*x*x - 
        (841. + 1265.*a)*x*x*x + (308. + 347.*a)*x*x*x*x
        - (28. - 40.*a )*x*x*x*x*x + 12.*a*x*x*x*x*x*x ;
    
    A[4]= + 768. + (816. - 768.*a)*x + (1240. - 1232.*a)*x*x 
        - 4.*(415. + 2.*a)*x*x*x + (311. + 722.*a)*x*x*x*x
        + (145. - 267.*a)*x*x*x*x*x - (36. + 51.*a)*x*x*x*x*x*x + 20.*a*x*x*x*x*x*x*x ;
    
    A[5]=+ 328.*x - (536. + 900.*a)*x*x + (208. + 1584.*a + 670.*a*a )*x*x*x 
        - a*(668. + 1161.*a + 225.*a*a )*x*x*x*x
        + a*a*(479. + 362.*a + 28.*a*a )*x*x*x*x*x 
        - a*a*a*(143. + 42.*a)*x*x*x*x*x*x + 16.*a*a*a*a*x*x*x*x*x*x*x;
    
    A[6] =+ 32. - 4.*(44. - 9.*a)*x + (384. - 322.*a - 400.*a*a )*x*x 
        - (400. - 869.*a - 1126.*a*a - 696.*a*a*a )*x*x*x
        + 2.*(80. - 488.*a - 517.*a*a - 631.*a*a*a - 264.*a*a*a*a )*x*x*x*x
        + (48. + 394.*a + 269.*a*a + 190.*a*a*a + 882.*a*a*a*a + 196.*a*a*a*a*a )*x*x*x*x*x
        - (64. - 58.*a - 89.*a*a - 95.*a*a*a + 34.*a*a*a*a + 296.*a*a*a*a*a + 32.*a*a*a*a*a*a )*x*x*x*x*x*x
        + (16. - 59.*a - 79.*a*a + 256.*a*a*a - 239.*a*a*a*a 
        + 57.*a*a*a*a*a + 48.*a*a*a*a*a*a )*x*x*x*x*x*x*x
        + (1. - a)*(1. - a)*(1. - a)*a*a*(29. + 16.*a)*x*x*x*x*x*x*x*x ;
    
    A[7]= + 28.*a*a*x*x - 32.*a*a*a*x*x*x;
    
    A[8]=- 288. + 36.*(1. + 8.*a)*x + 6.*(647. + 87.*a)*x*x + 5.*(55. - 927.*a - 132.*a*a )*x*x*x
        - (1233. + 98.*a - 879.*a*a - 192.*a*a*a )*x*x*x*x 
        + (360. + 1371.*a - 315.*a*a - 264.*a*a*a )*x*x*x*x*x
        - 24.*a*(17. - 4.*a*a)*x*x*x*x*x*x;
    
    A[9]=+ 32. + 4.*(-44. + 29.*a)*x - 12.*(-32. + 77.*a + 31.*a*a )*x*x
        + 2.*(-200. + 837.*a + 767.*a*a + 182.*a*a*a )*x*x*x 
        - 2.*(-80. + 625.*a + 905.*a*a + 520.*a*a*a + 82.*a*a*a*a )*x*x*x*x
        + (48. + 1079.*a + 590.*a*a + 1002.*a*a*a + 462.*a*a*a*a + 32.*a*a*a*a*a )*x*x*x*x*x
        + (-64. - 1160.*a - 501.*a*a - 364.*a*a*a - 486.*a*a*a*a - 72.*a*a*a*a*a )*x*x*x*x*x*x
        + (16. + 729.*a + 1038.*a*a + 38.*a*a*a + 238.*a*a*a*a + 52.*a*a*a*a*a )*x*x*x*x*x*x*x
        - a*(192. + 743.*a + 50.*a*a*a + 12.*a*a*a*a )*x*x*x*x*x*x*x*x + 192.*a*a*x*x*x*x*x*x*x*x*x;
    
    A[10]=+ 16.*x + 324.*x*x - 36.*x*x*x*x;
    
    A[11]=+ 216.*x - 672.*x*x + 152.*x*x*x;
    
    A[12]=- 16.*x + (16. - 42.*a)*x*x + (16. + 21.*a + 60.*a*a )*x*x*x
        - (16 - 21.*a + 45.*a*a + 32.*a*a*a )*x*x*x*x - a*a*(7. - 24.*a)*x*x*x*x*x; 
    
    A[13]=- 32. + (144. - 68.*a)*x + (-240. + 334.*a + 332.*a*a )*x*x 
        + (160. - 551.*a - 660.*a*a - 364.*a*a*a )*x*x*x
        + a*(329. + 451.*a + 650.*a*a + 164.*a*a*a )*x*x*x*x 
        + (-48. - a - 59.*a*a - 523.*a*a*a - 316.*a*a*a*a - 32.*a*a*a*a*a )*x*x*x*x*x
        + (16. - 43.*a - 93.*a*a + 255.*a*a*a + 287.*a*a*a*a + 32.*a*a*a*a*a )*x*x*x*x*x*x 
        - a*a*(-29. + 42.*a + 103.*a*a + 8.*a*a*a )*x*x*x*x*x*x*x;

    A[14]=- 144.*(1. - a)*(1. - a)*x*x + 144.*(1. - a)*(1. - a)*x*x*x - 36.*(1. - a)*(1. - a)*x*x*x*x;
    
    A[15]=- 32. + 96.*a + (48. - 32.*a)*x - 176.*a*x*x - (16. - 74.*a)*x*x*x + 212.*a*x*x*x*x;
    
    A[16]=- 32. + (64. - 100.*a)*x - 8.*(4. - 34.*a - 29.*a*a )*x*x 
        - 4.*a*(34. + 170.*a + 33.*a*a )*x*x*x
        + 8.*a*a*(47. + 51.*a + 4.*a*a)*x*x*x*x - 16.*a*a*a*(15. + 4.*a)*x*x*x*x*x 
        + 32.*a*a*a*a*x*x*x*x*x*x;
    
    ////////////////////////////////////////////////////////////////////////////
    
    C[0] = 1./(3.*a*(x-1.)*(x-1.)*x);
    
    C[1] = phi1(0.25)/((x-1.)*(x-1.)*(x-1.)*(a*x-1.));
    
    C[2] = phi1(0.25*a)/(2.*(x-1.)*(x-1.)*(x-1.)*(a*x-1.));
    
    C[3] = phi1(1./4./x)/(2.*(x-1.)*(x-1.)*(x-1.)*(a*x-1.));
    
    C[4] = phi1(0.25*x)/(2.*(x-1.)*(x-1.)*(x-1.)*(a*x-1.));
    
    C[5] = phi1(a*x*0.25)/((x-1.)*(x-1.)*(x-1.)*(a*x-1.));
    
    C[6] = phi2(1/a/x, 1/a)/(2.*a*a*x*x*(x-1.)*(x-1.)*(x-1.)*(a*x-1.));
    
    C[7] = log(a)*log(a)/(a*x-1.);
    
    C[8] = log(x)/((x-1.)*(x-1.)*(x-1.)*(a*x-1.)*3.);
    
    C[9] = log(x)*log(x)/((x-1.)*(x-1.)*(x-1.)*(x-1.)*(a*x-1.)*2.*a*x);
    
    C[10] = 2.*log(mu/SM.Mw_tree())/((x-1.)*(x-1.));
    
    C[11] = log(x)*2.*log(mu/SM.Mw_tree())/((x-1.)*(x-1.)*(x-1.));
    
    C[12] = log(a)/((x-1.)*(x-1.)*(a*x-1.));
    
    C[13] = log(x)*log(a)/(2.*a*(x-1.)*(x-1.)*(x-1.)*x*(a*x-1.));
    
    C[14] = gsl_sf_dilog(1.-a)/(x-1.)/(x-1.);
    
    C[15] = gsl_sf_dilog(1.-x)/a*x;
    
    C[16] = gsl_sf_dilog(1.-a*x)/(a*x*(x-1.)*(x-1.));
    
    for (int i=0; i<10; i++){
        b += C[i]*A[i];
    }
    
    return (b/128./swsq);
}

double StandardModelMatching::phi1(double z) const{
    if (z >= 0.) {
        if (z < 1){
                return(4.*sqrt(z/(1.-z))*gsl_sf_clausen(2.*asin(sqrt(z))));
        }
        else{
            return((1./sqrt(1.-1./z))*(2.*log(0.5*sqrt(1.-1./z))*log(0.5*sqrt(1.-1./z))-
                    4.*gsl_sf_dilog(0.5*(1.-sqrt(z/(1.-z)))) - log(4.*z)*log(4.*z)
                    + M_PI*M_PI/3.));
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
    double l = sqrt((1.-x-y)*(1.-x-y)-4.*x*y);
    
    if ((l*l) >= 0. || (sqrt(x) + sqrt(y)) <= 1.){
        return( 1./l*(M_PI*M_PI/3. +2.*log(0.5*(1.+x-y-l))*log(0.5*(1.-x+y-l)) -log(x)*log(y)
                -2.*gsl_sf_dilog(0.5*(1.+x-y-l)) -2.*gsl_sf_dilog(0.5*(1.-x+y-l))));        
    }
    else if((l*l) < 0. || (sqrt(x) + sqrt(y)) > 1.){
        return(2./(-l*l)*(gsl_sf_clausen(2.*acos((-1.+x+y)/(2.*sqrt(x*y))))
                +gsl_sf_clausen(2.*acos((1.+x-y)/(2.*sqrt(x))))
                +gsl_sf_clausen(2.*acos((1.-x+y)/(2.*sqrt(y))))));
    }
    else{
        std::stringstream out;
        out << x;
        throw std::runtime_error("StandardModelMatching::phi2(double x, double y) wrong" + out.str());
    }
    return(0.);
}

double StandardModelMatching::Y0t(double x) const{
    return( x/8.*((4.-x)/(1.-x) + 3.*x*log(x)*(1.-x)*(1.-x)) );
}

double StandardModelMatching::Y1t(double x) const{
    return( (4.*x+16.*x*x+4.*x*x*x)/(3.*(1.-x)*(1.-x))
           -(4.*x-10.*x*x+x*x*x-x*x*x*x)*log(x)/(1.-x)/(1.-x)/(1.-x)
           +(2.*x-14.*x*x+x*x*x-x*x*x*x)*log(x)*log(x)/2./(1.-x)/(1.-x)
           +(2.*x+x*x*x)/(1.-x)/(1.-x)*gsl_sf_dilog(1.-x) 
           +x*((4.-x)/(1.-x) + 3.*x*log(x)*(1.-x)*(1.-x))*0.//log(mu^2/Muw^2)
           +x*x*(-1./(1.-x) + (4.-x)/(1.-x)/(1.-x) + 3.*log(x)/(1.-x)/(1.-x)
                 +6.*x*log(x)/(1.-x)/(1.-x)/(1.-x) +3./(1.-x)/(1.-x))*0.////log(mu^2/Muw^2)
            );
}

/******************************************************************************/
const std::vector<WilsonCoefficient>& StandardModelMatching::CMdbd2() 
{   
//    if(SM_i == SM)
//        return(vmc);
    
    double gammam = 8.;                                                         
    double Bt;  
    
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)  
            / SM.Mw_tree(), 2.); // always FULLNLO
    complex co = SM.getGF() / 4. / M_PI * SM.Mw_tree() * SM.getlamt_d();
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

    mcdbd2.setMu(SM.getMuw());
    
    switch (mcdbd2.getOrder()) {
        case NNLO:
        case NLO:
            mcdbd2.setCoeff(0, co * co * 4. * (SM.Als(SM.getMuw()) / 4. / M_PI * (S1(xt) +
                    Bt * S0(xt, xt) + 2. * gammam * S0p(xt) * log(SM.getMuw() / SM.Mw_tree()))), NLO);   
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
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); // always FULLNLO
    complex co = SM.getGF() / 4. / M_PI * SM.Mw_tree() * SM.getlamt_s();
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

    mcdbs2.setMu(SM.getMuw());
 
    switch (mcdbs2.getOrder()) {
        case NNLO:
        case NLO:          
            mcdbs2.setCoeff(0, co * co * 4. * (SM.Als(SM.getMuw()) / 4. / M_PI * (S1(xt) +
                    Bt * S0(xt, xt) + 2. * gammam * S0p(xt) * log(SM.getMuw() / SM.Mw_tree()))), NLO);      
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
    
    mcdk2.setMu(SM.getMut());
 
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

    mcd1Buras.setMu(SM.getMuw());
    
    switch (mcd1Buras.getOrder()) {
        case NNLO:
        case NLO:
            mcd1Buras.setCoeff(0, SM.Als(SM.getMuw()) / 4. / M_PI  * 11./2. , NLO);
            mcd1Buras.setCoeff(1, SM.Als(SM.getMuw()) / 4. / M_PI * (-11./6.) , NLO);
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

    mcd1.setMu(SM.getMuw());
    
    switch (mcd1.getOrder()) {
        case NNLO:
        case NLO:
            mcd1.setCoeff(0, SM.Als(SM.getMuw()) / 4. / M_PI  * 15. , NLO);
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

    mcdd2.setMu(SM.getMuw());
 
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), 
            SM.getQuarks(QCD::TOP).getMass(), 5.) / SM.Mw_tree(), 2.); // always FULLNLO
    
    vmck.clear();
    
    switch (mck.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mck.getScheme();
            throw "StandardModel::CMbsg(): scheme " + out.str() + "not implemented";
    }

    mck.setMu(SM.getMuw());
    
    switch (mck.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<10; j++){
                mck.setCoeff(j, SM.getlamt() * SM.Als(SM.getMuw()) / 4. / M_PI * 
                                setWCbnlep(j, xt,  NLO), NLO);
                mck.setCoeff(j, SM.getlamt() * SM.getAle() / 4. / M_PI *
                                setWCbnlepEW(j, xt), NLO_ew);
                }
        case LO:
            for (int j=0; j<10; j++){
                mck.setCoeff(j, SM.getlamt() *  setWCbnlep(j, xt,  LO), LO);
                mck.setCoeff(j, 0., LO_ew); 
                }                   
            break;
        default:
            std::stringstream out;
            out << mck.getOrder();
            throw "StandardModelMatching::CMbsg(): order " + out.str() + "not implemented";
    }

    vmck.push_back(mck);
    return(vmck);
}

/*******************************************************************************
 * Wilson coefficients Buras base for K -> pi pi decays                        * 
 * operator basis: - current current                                           *
 * ****************************************************************************/
const std::vector<WilsonCoefficient>& StandardModelMatching::CMKCC(){
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); // always FULLNLO
    
    vmckcc.clear();
    
    switch (mckcc.getScheme()) {
        case NDR:
        //case HV:
        //case LRI:
        break;
        default:
            std::stringstream out;
            out << mckcc.getScheme();
            throw "StandardModel::CMbsg(): scheme " + out.str() + "not implemented";
    }

    mckcc.setMu(SM.getMuw());
    
    switch (mckcc.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<2; j++){
                mckcc.setCoeff(j, SM.getlamt() * setWCbnlep(j, xt, NLO), NLO); 
            }
            for (int j=2; j<10; j++){
                mckcc.setCoeff(j, 0. , NLO); 
            }
        case LO:
            for (int j=0; j<2; j++){
                mckcc.setCoeff(j, SM.getlamt() * setWCbnlep(j, xt, LO), LO); 
            }
            for (int j=2; j<10; j++){
                mckcc.setCoeff(j, 0. , LO); 
            }
            break;
        default:
            std::stringstream out;
            out << mckcc.getOrder();
            throw "StandardModelMatching::CMbsg(): order " + out.str() + "not implemented";
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
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass())
            / SM.Mw_tree(), 2.);
    complex co = (- 4. * SM.getGF() / sqrt(2)) * SM.getlamt_s();
    
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

    mcbsg.setMu(SM.getMuw());
    
    switch (mcbsg.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<10; j++){
            mcbsg.setCoeff(j, co * SM.Als(SM.getMuw()) / 4. / M_PI * setWCbsg(j, xt,  NLO) , NLO);
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
    sw =  sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * SM.Mw_tree() * SM.Mw_tree()) ) ;

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
                                B0t(x) - D0t(x) + 38./27. + 1/(4.*sw*sw);
            CWbsgArrayNLO[9] = 1./(sw*sw) * (B0t(x) - C0t(x)) -1/(4.*sw*sw);
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
                / SM.Mw_tree(), 2.);
    double co = ( SM.getGF() / sqrt(2));
    
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

    mcbnlep.setMu(SM.getMuw());
    
    switch (mcbnlep.getOrder()) {
        case NNLO:
        case NLO:
            for (int j=0; j<10; j++){
                mcbnlep.setCoeff(j,co * lambda * SM.Als(SM.getMuw()) / 4. / M_PI * 
                                setWCbnlep(j, xt,  NLO), NLO);
                mcbnlep.setCoeff(j, co * lambda * SM.getAle() / 4. / M_PI *
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
    matrix<complex> ckm = SM.getVCKM();
    
    switch (a) {
        case 0: lambda1 = SM.getlamu_d();
                break;
        case 1: lambda1 = SM.getlamu_s();
                break;
        case 2: lambda1 = ckm(0,2).conjugate()*ckm(1,0);
                break;
        case 3: lambda1 = ckm(1,2).conjugate()*ckm(0,0); 
                break;
        case 4: lambda1 = ckm(0,2).conjugate()*ckm(1,1);
                break;
        case 5: lambda1 = ckm(1,2).conjugate()*ckm(2,1);
                break;
        default:
            std::stringstream out;
            out << a;
            throw std::runtime_error("case" + out.str() + "unexsting; implemented i=0,1,2,3"); 
    }
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
                / SM.Mw_tree(), 2.);
    double co = ( SM.getGF() / sqrt(2));
    
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

    mcbnlepCC.setMu(SM.getMuw());
    
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); // always FULLNLO
    double a = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.getMHl(), 2.);
    
    vmckpnn.clear();
    
    mckpnn.setMu(SM.getMut());
 
    switch (mckpnn.getOrder()) {
        case NNLO:
        case NLO:
            mckpnn.setCoeff(0, SM.Als(SM.getMuw())/4./M_PI*SM.getlamt().imag()*X1t(xt)/SM.GetLambda(), NLO);
        case LO:
            mckpnn.setCoeff(0, SM.getlamt().imag()*X0t(xt)/SM.GetLambda(), LO);
            break;
        default:
            std::stringstream out;
            out << mckpnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkp0nn(): order " + out.str() + "not implemented"); 
    }
    
    switch (mckpnn.getOrder_ew()) {
        case NLO_ew:
            mckpnn.setCoeff(0, SM.getAle()/4./M_PI*SM.getlamt().imag()*Xewt(xt, a, SM.getMuw())/SM.GetLambda(), NLO_ew);
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); // always FULLNLO
    double a = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.getMHl(), 2.);
    
    vmckppnn.clear();
    
    mckppnn.setMu(SM.getMut());
 
    switch (mckppnn.getOrder()) {
        case NNLO:
        case NLO:
            mckppnn.setCoeff(0, SM.Als(SM.getMuw())/4./M_PI*SM.getlamt().imag()*X1t(xt)/SM.GetLambda(), NLO);
        case LO:
            mckppnn.setCoeff(0, SM.getlamt().imag()*X0t(xt)/SM.GetLambda(), LO);
            break;
        default:
            std::stringstream out;
            out << mckppnn.getOrder();
            throw std::runtime_error("StandardModelMatching::CMkppnn(): order " + out.str() + "not implemented"); 
    }
    
    switch (mckppnn.getOrder_ew()) {
        case NLO_ew:
            mckppnn.setCoeff(0, SM.getAle()/4./M_PI*SM.getlamt().imag()*Xewt(xt, a, SM.getMuw())/SM.GetLambda(), NLO_ew);
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); // always FULLNLO
    
    
    vmckmm.clear();
    
    mckmm.setMu(SM.getMut());
 
    switch (mckmm.getOrder()) {
        case NNLO:
        case NLO:
            mckmm.setCoeff(0, SM.Als(SM.getMuw())/4./M_PI*SM.getlamt().real()*Y1t(xt)/SM.GetLambda(), NLO);
        case LO:
            mckmm.setCoeff(0, SM.getlamt().real()*Y0t(xt)/SM.GetLambda(), LO);
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); // always FULLNLO
    
    vmcbsnn.clear();
    
    mcbsnn.setMu(SM.getMut());
 
    switch (mcbsnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (SM.getVCKM()(2,1).abs() / SM.getVCKM()(1,2).abs())*
                                SM.Als(SM.getMuw())/4./M_PI*X1t(xt), NLO);
        case LO:
            mcbsnn.setCoeff(0,  (SM.getVCKM()(2,1).abs() / SM.getVCKM()(1,2).abs())*
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); 
    
    vmcbdnn.clear();
    
    mcbdnn.setMu(SM.getMut());
 
    switch (mcbdnn.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (SM.getVCKM()(2,2).abs() / SM.getVCKM()(1,2).abs())*
                                SM.Als(SM.getMuw())/4./M_PI*X1t(xt), NLO);
        case LO:
            mcbsnn.setCoeff(0,  (SM.getVCKM()(2,2).abs() / SM.getVCKM()(1,2).abs())*
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); 
    
    vmcbsmm.clear();
    
    mcbsmm.setMu(SM.getMut());
 
    switch (mcbsmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (SM.getVCKM()(2,2).conjugate() * SM.getVCKM()(2,1)).abs()*
                                SM.Als(SM.getMuw())/4./M_PI*Y1t(xt), NLO);
        case LO:
            mcbsnn.setCoeff(0,  (SM.getVCKM()(2,2).conjugate() * SM.getVCKM()(2,1)).abs()*
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
    
    double xt = pow(SM.Mrun(SM.getMuw(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
            / SM.Mw_tree(), 2.); 
    
    vmcbdmm.clear();
    
    mcbdmm.setMu(SM.getMut());
 
    switch (mcbdmm.getOrder()) {
        case NNLO:
        case NLO:
            mcbsnn.setCoeff( 0, (SM.getVCKM()(2,2).conjugate() * SM.getVCKM()(2,0)).abs()*
                                SM.Als(SM.getMuw())/4./M_PI*Y1t(xt), NLO);
        case LO:
            mcbsnn.setCoeff(0,  (SM.getVCKM()(2,2).conjugate() * SM.getVCKM()(2,0)).abs()*
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
    sw =  sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * SM.Mw_tree() * SM.Mw_tree()) );
    
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
    sw =  sqrt( (M_PI * SM.getAle() )/( sqrt(2) * SM.getGF() * SM.Mw_tree() * SM.Mw_tree()) ) ;
    
    if ( swb == sw && xcacheb == x){
        return (CWbnlepArrayNLOew[i]);
    }
    
    swc = sw; xcachec = x;
    
    CWbnlepArrayNLOew[1] = - 35./18.;
    CWbnlepArrayNLOew[2] = 2./(3. *sw*sw) * ( 2. *B0b(x) + C0b(x) );
    CWbnlepArrayNLOew[6] = 2./3. * (4. *C0b(x) + D0b(x) - 4./9.);
    CWbnlepArrayNLOew[8] = 2./3. * (4.*C0b(x) + D0b(x) - 4./9. + (1./(sw*sw)) *
                           (10.*B0b(x) - 4.*C0b(x)) );
    
    return (CWbnlepArrayNLOew[i]);
}

complex StandardModelMatching::S0c() const 
{
    double xc = pow(SM.Mrun(SM.getMuc(), SM.getQuarks(QCD::CHARM).getMass(), SM.getQuarks(QCD::CHARM).getMass(), 4.)
                / SM.Mw_tree(), 2.);
    complex co = SM.getGF() / 2. / M_PI * SM.Mw_tree() * SM.getlamc().conjugate();
    
    return(co*co*xc);    
}

complex StandardModelMatching::S0ct() const 
{
    double xc = pow(SM.Mrun(SM.getMuc(), SM.getQuarks(QCD::CHARM).getMass(), SM.getQuarks(QCD::CHARM).getMass(), 4.)
                / SM.Mw_tree(), 2.);
    double xt = pow(SM.Mrun(SM.getMut(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
                / SM.Mw_tree(), 2.);
    double co = SM.getGF() / 2. / M_PI * SM.Mw_tree();
    
    return( co*co*2.*SM.getlamc().conjugate()*SM.getlamt().conjugate()* xc *
            (log(xt/xc) - 3.*xt/4./(1.-xt) - 3.*xt*xt*log(xt)/4./(1.-xt)/(1.-xt)) );
}

complex StandardModelMatching::S0tt() const
{
    double x = pow(SM.Mrun(SM.getMut(), SM.getQuarks(QCD::TOP).getMass(), SM.getQuarks(QCD::TOP).getMass(), 5.)
                / SM.Mw_tree(), 2.);
    complex co = SM.getGF() / 2. / M_PI * SM.Mw_tree() * SM.getlamt().conjugate();
    
    return ( co * co * ((4.*x - 11.*x*x + x*x*x)/4./(1.-x)/(1.-x) - 
            3.*x*x*x/2./(1.-x)/(1.-x)/(1.-x) * log(x)) );
}
