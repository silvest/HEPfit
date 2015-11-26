/* 
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Flavour.h"
#include "bsgamma.h"
#include <gslpp_complex.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>
#include <boost/bind.hpp>

Bsgamma::Bsgamma(const StandardModel& SM_i, int obsFlag)
: ThObservable(SM_i),
Intbc_cache(2, 0.)
{
    if (SM.ModelName().compare("StandardModel") != 0 && SM.ModelName().compare("FlavourWilsonCoefficient") != 0) std::cout << "\nWARNING: b to s gamma not implemented in: " + SM.ModelName() + " model, returning Standard Model value.\n" << std::endl;
    
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR) or 2 (BR_CPodd)");
    
    Intb1Cached = 0;
    Intb2Cached = 0;
    Intb3Cached = 0;
    Intb4Cached = 0;
    Intbb1Cached = 0;
    Intbb2Cached = 0;
    Intbb4Cached = 0;
    Intbc1Cached = 0;
    Intbc2Cached = 0;
    Intc1Cached = 0;
    Intc2Cached = 0;
    Intc3Cached = 0;
    Intcc1Cached = 0;
    
    w_INT = gsl_integration_cquad_workspace_alloc (100);
}

void Bsgamma::checkCache()
{
    if (Mb_kin == Intb_cache)
        Intb_updated = 1;
    else {
        Intb_cache = Mb_kin;
        Intb_updated = 0;
    }
    
    if (Mb_kin == Intbc_cache(0) && Mc == Intbc_cache(1)) 
        Intbc_updated = 1;
    else {
        Intbc_cache(0) = Mb_kin;
        Intbc_cache(1) = Mc;
        Intbc_updated = 0;
    }
}

double Bsgamma::delta(double E0)
{
    return 1. - 2.*E0/Mb_kin;
}

double Bsgamma::rho(double E0)
{
    double d=delta(E0);
    double d4=d*d*d*d;
    
    return d + d4/6. + log(1. - d);
}

double Bsgamma::omega(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return 3./2. * d2 - 2. * d3 + d4;
}

double Bsgamma::T1(double E0,double t)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    double Li2 = gsl_sf_dilog(d);
    
    return 109./18. * d + 17./18. * d2 - 191./108. * d3 + 23./16. * d4 
            + 79./18. * log(1. - d) - 5./3. * Li2 
            - (5./3. * rho(E0) + 2./9. * omega(E0)) * log(t*d) 
            /* + rho(E0)/9.*log(ms^5/(mu^4*md)) */;
}

double Bsgamma::T2(double E0,double t)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    double Li2 = gsl_sf_dilog(d);
    
    return 187./108. * d + 7./18. * d2 - 395./648. * d3 + 1181./2592. * d4 
            + 133./108. * log(1. - d) - Li2/2. 
            - (rho(E0)/2. + 2./27. * omega(E0)) * log(t*d) 
            /* + rho(E0)/9.*log(ms/mu) */;
}

double Bsgamma::T3(double E0,double t)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    double Li2 = gsl_sf_dilog(d);
    
    return 35./162. * d + 1./72. * d2 - 89./1944. * d3 + 341./7776. * d4 
            + 13./81. * log(1. - d) - Li2/18. 
            - (rho(E0)/18. + omega(E0)/162.) * log(t*d);
}

double Bsgamma::P0_4body(double E0, double t)
{
    gslpp::complex A1 =-C1_0*CKMu;
    gslpp::complex A2 =-C2_0*CKMu;
    
    return (C3_0.abs2() + 20. * (C3_0*C5_0).real() + 2./9. * C4_0.abs2() 
            + 40./9. * (C4_0*C6_0).real() + 136. * C5_0.abs2() 
            + 272./9. * C6_0.abs2()) * T1(E0,t)
            + (2./9. * A1.abs2() + A2.abs2() 
            + (8./9. * C3_0.real() 
            - 4./27. * C4_0.real() + 128./9. * C5_0.real() 
            - 64./27. * C6_0.real()) * A1.real() 
            +(2./3. * C3_0.real() + 8./9. * C4_0.real() + 32./3. * C5_0.real() 
            + 128./9. * C6_0.real()) * A2.real()) * T2(E0,t)
            + (C3_0.abs2() + 8./3. * (C3_0*C4_0).real() + 32. * (C3_0*C5_0).real() 
            + 128./3. * (C3_0*C6_0).real() - 2./9. * C4_0.abs2() 
            + 128./3. * (C4_0*C5_0).real() - 64./9. * (C4_0*C6_0).real()
            + 256. * C5_0.abs2() + 2048./3 * (C5_0*C6_0).real() 
            - 512./9. * C6_0.abs2()) * T3(E0,t);
}

double Bsgamma::zeta()
{
    return Mc*Mc/Mb_kin/Mb_kin;
}

gslpp::complex Bsgamma::a(double z)
{
    double zeta3 = gsl_sf_zeta_int(3);
    
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double z6=z5*z;
    
    double L=log(z);
    double L2=L*L;
    double L3=L2*L;
    
    double pi2=M_PI*M_PI;
    
    if (z == 1.) return 4.0859 + 4./9. * M_PI * gslpp::complex::i();
    else return 16./9. * ( ( 5./2. - pi2/3. - 3.*zeta3 
            + ( 5./2. - 3./4.*pi2 )*L + L2/4. + L3/12. )*z 
            + ( 7./4. + 2./3.*pi2 - pi2*L/2. - L2/4. + L3/12. )*z2 
            + ( -7./6. - pi2/4. + 2*L - 3./4.*L2 )*z3 
            + ( 457./216. - 5./18*pi2 - L/72. - 5./6.*L2 )*z4 
            + ( 35101./8640. - 35./72.*pi2 - 185./144.*L - 35./24.*L2 )*z5
            + ( 67801./8000. - 21./20.*pi2 - 3303./800.*L - 63./20.*L2 )*z6 + 
            gslpp::complex::i()*M_PI*( ( 2. - pi2/6. + L/2. + L2/2. )*z 
            + ( 1./2. - pi2/6. - L + L2/2. )*z2
            + z3 + 5./9.*z4 + 49./72.*z5 + 231./200.*z6) );
}

gslpp::complex Bsgamma::b(double z)
{
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double z6=z5*z;
    
    double L=log(z);
    double L2=L*L;
    
    double pi2=M_PI*M_PI;
    
    if (z == 1.) return 0.0316 + 4./81. * M_PI * gslpp::complex::i();
    else return -8./9. * ( ( -3. + pi2/6. - L )*z - 2./3.*pi2*pow(z,3./2.) 
            + ( 1./2. + pi2 -2.*L - L2/2. )*z2 
            + ( -25./12. - pi2/9. - 19./18.*L + 2.*L2 )*z3 
            + ( -1376./225. + 137./30.*L + 2.*L2 + 2./3.*pi2 )*z4 
            + ( -131317./11760. + 887./84.*L + 5.*L2 + 5./3.*pi2 )*z5 
            + ( -2807617./97200. + 16597./540.*L + 14.*L2 + 14./3.*pi2 )*z6 + 
            gslpp::complex::i()*M_PI*( -z + ( 1 - 2.*L )*z2 
            + ( -10./9. + 4./3.*L )*z3 + z4 + 2./3.*z5 + 7./9.*z6) );
}

gslpp::complex Bsgamma::r1(int i, double z)
{
    double Xb = -0.16844083981858157;
    
    switch(i){
        case 1:
            return 833./729. - (a(z) + b(z))/3. + 40./243.*gslpp::complex::i()*M_PI;
        case 2:
            return - 1666./243. + 2.*(a(z) + b(z)) - 80./81.*gslpp::complex::i()*M_PI;
        case 3:
            return 2392./243. + 8.*M_PI/3./sqrt(3.) + 32./9.*Xb - a(1.) + 2.*b(1.) + 56./81.*gslpp::complex::i()*M_PI;
        case 4:
            return -761./729. - 4.*M_PI/9./sqrt(3.) - 16./27.*Xb + a(1.)/6. + 5.*b(1.)/3. + 2.*b(z) - 148./243.*gslpp::complex::i()*M_PI;
        case 5:
            return 56680./243. + 32.*M_PI/3./sqrt(3.) + 128./9.*Xb - 16.*a(1.) + 32.*b(1.) + 896./81.*gslpp::complex::i()*M_PI;
        case 6:
            return 5710./729. - 16.*M_PI/9./sqrt(3.) - 64./27.*Xb - 10./3.*a(1.) + 44./3.*b(1.) + 12.*a(z) + 20.*b(z) 
                    - 2296./243.*gslpp::complex::i()*M_PI;
        case 8:
            return 44./9. - 8./27.*M_PI*M_PI + 8./9.*gslpp::complex::i()*M_PI;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("Bsgamma::r1(): index " + out.str() + " not implemented");
    }
}

gslpp::complex Bsgamma::Gamma_t(double t)
{
    if (t<4) return -2. * atan( sqrt(t/(4.-t)) ) * atan( sqrt(t/(4.-t)) );
    else return -M_PI*M_PI/2. + 2.*log( ( sqrt(t) + sqrt(t-4.) ) / 2. )*log( ( sqrt(t) + sqrt(t-4.) ) / 2. )
            - 2.*gslpp::complex::i()*M_PI*log( ( sqrt(t) + sqrt(t-4.) ) / 2. );
}


gslpp::complex Bsgamma::kappa(double Mq,double t)
{
    double s = t * Mb_kin*Mb_kin/Mq/Mq;
    return 1./2. + Gamma_t(s)/s;
}

double Bsgamma::Int_b1(double E0)
{
    if (Intb1Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntb1 = delta(E0)*mt + mt2;
        Intb1Cached = 1;
    }

    return CacheIntb1;
}

double Bsgamma::Int_b2(double E0)
{
    if (Intb2Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntb2 = delta(E0)*mt + mt2;
        Intb2Cached = 1;
    }

    return CacheIntb2;
}

double Bsgamma::Int_b3(double E0)
{
    if (Intb3Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_t, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double t = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        CacheIntb3 = delta(E0)*t + mt;
        Intb3Cached = 1;
    }

    return CacheIntb3;
}

double Bsgamma::Int_b4(double E0)
{
    if (Intb4Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_t2_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_re_t2_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntb4 = delta(E0)*mt + mt2;
        Intb4Cached = 1;
    }

    return CacheIntb4;
}

double Bsgamma::Int_bb1(double E0)
{
    if (Intbb1Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_abs2_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_abs2_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntbb1 = delta(E0)*mt + mt2;
        Intbb1Cached = 1;
    }

    return CacheIntbb1;
}

double Bsgamma::Int_bb2(double E0)
{
    if (Intbb2Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_abs2_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_abs2_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntbb2 = delta(E0)*mt + mt2;
        Intbb2Cached = 1;
    }

    return CacheIntbb2;
}

double Bsgamma::Int_bb4(double E0)
{
    if (Intbb4Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_abs2_t2_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKb_abs2_t2_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntbb4 = delta(E0)*mt + mt2;
        Intbb4Cached = 1;
    }

    return CacheIntbb4;
}

double Bsgamma::Int_bc1(double E0)
{
    if (Intbc1Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntbc1 = delta(E0)*mt + mt2;
        Intbc1Cached = 1;
    }

    return CacheIntbc1;
}

double Bsgamma::Int_bc2(double E0)
{
    if (Intbc2Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntbc2 = delta(E0)*mt + mt2;
        Intbc2Cached = 1;
    }

    return CacheIntbc2;
}

double Bsgamma::Int_c1(double E0)
{
    if (Intc1Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntc1 = delta(E0)*mt + mt2;
        Intc1Cached = 1;
    }

    return CacheIntc1;
}

double Bsgamma::Int_c1_im(double E0)
{
    if (Intc1imCached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntc1im = delta(E0)*mt + mt2;
        Intc1imCached = 1;
    }

    return CacheIntc1im;
}

double Bsgamma::Int_c2(double E0)
{
    if (Intc2Cached == 0) {
        double t1 = (1. - delta(E0));

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;

        CacheIntc2 = delta(E0)*mt + mt2;
        Intc2Cached = 1;
    }

    return CacheIntc2;
}

double Bsgamma::Int_c3(double E0)
{
    if (Intc3Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double t = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;
        
        CacheIntc3 = delta(E0)*t + mt;
        Intc3Cached = 1;
    }

    return CacheIntc3;
}

double Bsgamma::Int_cc(double E0)
{
    if (IntccCached == 0) {
        double t1 = (1. - delta(E0));

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_abs2_t, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_abs2_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;
        
        CacheIntcc = delta(E0)*mt + 2. * mt2;
        IntccCached = 1;
    }
    
    return CacheIntcc;
}

double Bsgamma::Int_cc1(double E0)
{
    if (Intcc1Cached == 0) {
        double t1 = (1. - delta(E0));

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_abs2_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_abs2_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        double mt2 = avaINT;
        
        CacheIntcc1 = delta(E0)*mt + mt2;
        Intcc1Cached = 1;
    }
    
    return CacheIntcc1;
}

double Bsgamma::Int_cc1_part1(double E0)
{
    if (Intcc1p1Cached == 0) {
        double t1 = (1. - delta(E0));

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_abs2_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        
        CacheIntcc1p1 = delta(E0)*avaINT;
        Intcc1p1Cached = 1;
    }
    
    return CacheIntcc1p1;
}

double Bsgamma::ff7_dMP(double E0)
{
    if (FOUR_BODY){
        double d=delta(E0);
        double d2=d*d;
        double d3=d2*d;

        return 4. * d * (18. - 33.*d + 2.*d2 + 13.*d3 - 6.* d2 * (2. + d) * log(d)) 
                / (81. * (d - 1.));   
    }
    else 
        return 0.;
}

double Bsgamma::ff7_sMP(double E0)
{
    if (FOUR_BODY){
        double d=delta(E0);
        double d2=d*d;
        double d3=d2*d;

        return (-2. * d * (72. + 39.*d - 76.*d2 - 35.*d3 
                + 6.*d*(18. + 13.*d + 2.*d2)*log(d))) / (243.*(d - 1.));
    }
    else 
        return 0.;
}

double Bsgamma::ff8_dMP(double E0)
{
    if (FOUR_BODY){
        double d=delta(E0);
        double d2=d*d;
        double d3=d2*d;
        double ld = log(d);
        double l1d = log(1. - d);
        double Li2 = gsl_sf_dilog(d);

        return -136./27. * d - 724./81. * d2 + 20./27. * d3 
                + (-8./9. + 16./9. * d - 8./9. * d2) * l1d* l1d 
                + (32./27. * d + 76./27. * d2 - 16./81. * d3) * ld 
                + (-104./27. - 80./9. * d + 40./9. * d2 + (32./27. 
                + 32./9. * d - 16./9. * d2) * ld) * l1d 
                + (-64./27. * d - 152./27. * d2 + 32./81. * d3 
                + (-64./27. - 64./9. * d + 32./9. * d2) * l1d) * log(Ms/Mb_kin) 
                + (32./27. + 32./9. * d - 16./9. * d2) * Li2;
    }
    else 
        return 0.;
}

double Bsgamma::ff8_sMP(double E0)
{
    if (FOUR_BODY){
        double d=delta(E0);
        double d2=d*d;
        double d3=d2*d;
        double ld = log(d);
        double l1d = log(1. - d);
        double Li2 = gsl_sf_dilog(d);

        return -340./243. * d - 104./81. * d2 + 16./729. * d3 
                + (-4./27. + 8./27. * d - 4./27. * d2) * l1d* l1d 
                + (8./27. * d + 4./9. * d2) * ld 
                + (-16./27. * d - 8./9. * d2) * log(Ms/Mb_kin)
                + (-268./243. - 40./27. * d + 20./27. * d2 + (8./27. 
                + 16./27. * d - 8./27. * d2) * ld
                + (-16./27. - 32./27. * d + 16./27. * d2) * log(Ms/Mb_kin)) * l1d 
                + (8./27. + 16./27. * d - 8./27. * d2) * Li2;
    }
    else 
        return 0.;
}

double Bsgamma::Phi11_1(double E0)
{
    return Phi22_1(E0)/36.;
}

double Bsgamma::Phi12_1(double E0)
{
    return -Phi22_1(E0)/3.;
}

double Bsgamma::Phi13_1(double E0)
{
    return -Phi23_1(E0)/6.;
}

double Bsgamma::Phi14_1(double E0)
{
    return -Phi24_1(E0)/6.;
}

double Bsgamma::Phi15_1(double E0)
{
    return -Phi25_1(E0)/6.;
}

double Bsgamma::Phi16_1(double E0)
{
    return -Phi26_1(E0)/6.;
}

double Bsgamma::Phi17_1(double E0, double z)
{
    return -Phi27_1(E0,z)/6.;
}

double Bsgamma::Phi18_1(double E0, double z)
{
    return Phi27_1(E0,z)/18.;
}

double Bsgamma::Phi22_1(double E0)
{
    return 16./27. * Int_cc1(E0);
}

double Bsgamma::Phi23_1_4body(double E0)
{
    if (FOUR_BODY)
        return 0.0039849625073434735;
    else
        return 0.;
}

double Bsgamma::Phi23_1(double E0)
{
    return -8./27. * (Int_c1(E0) + Int_c2(E0) + 2.*Int_bc1(E0) - 2.*Int_bc2(E0))
            - Phi23_1_4body(E0);
}

double Bsgamma::Phi24_1_4body(double E0)
{
    if (FOUR_BODY)
        return 0.012330977673588935;
    else
        return 0.;
}

double Bsgamma::Phi24_1(double E0)
{
    return -1./6. * (Phi23_1(E0) + Phi23_1_4body(E0))
            - Phi24_1_4body(E0);
}

double Bsgamma::Phi25_1_4body(double E0)
{
    if (FOUR_BODY)
        return 0.06375940011749558;
    else
        return 0.;
}

double Bsgamma::Phi25_1(double E0)
{
    return -32./27. * (4.*Int_c1(E0) + Int_c2(E0) + 8.*Int_bc1(E0) - 2.*Int_bc2(E0))
            - Phi25_1_4body(E0);
}

double Bsgamma::Phi26_1_4body(double E0)
{
    if (FOUR_BODY)
        return 0.11932481422855279;
    else
        return 0.;
}

double Bsgamma::Phi26_1(double E0)
{
    return 16./81. * (4.*Int_c1(E0) + Int_c2(E0) - 10.*Int_bc1(E0) - 2.*Int_bc2(E0) + 36.*Int_cc1(E0))
            - Phi26_1_4body(E0);
}

double Bsgamma::Phi27_1(double E0, double z)
{
    double d = delta(E0);
    double d2 = d*d;
    double Pi2 = M_PI*M_PI;
    double st0 = sqrt(1. - 4.*z);
    double std = sqrt( (1. - d - 4.*z) * (1. - d) );
    double L0 = log( ( 1. + st0 ) / ( 2.*sqrt(z) ) );
    double Ld = log( ( sqrt(1. - d) + sqrt(1. - d - 4.*z) ) / ( 2.*sqrt(z) ) );
    
    if (d == 1) {
        return -2./27. + (2.*Pi2 - 7.)/9. * z + 4.*(3. - 2.*Pi2)/9. * z * z
                + 4./3. * z * (1. - 2.*z) * st0 * L0
                - 8./9. * z * (6.*z*z - 4.*z + 1.) * L0*L0 + 4./3. * Pi2 * z * z *z;
    } else return -2./27. * d * (3. - 3.*d + d2) + (2.*Pi2 - 7.)/9. * z * d * (2. - d)
            + 4.*(3. - 2.*Pi2)/9. * z * z * d 
            + 4./3. * z * (1. - 2.*z) * ( st0 * L0 - std * Ld ) 
            + 4./3. * z * d * std * Ld 
            - 8./9. * z * (6.*z*z - 4.*z + 1.) * ( L0*L0 - Ld*Ld ) 
            - 8./9. * z * d * (2. - d - 4.*z) * Ld * Ld;
}

double Bsgamma::Phi27_1_im(double E0, double z)
{
    if (z >= 1./4.) 
        throw std::runtime_error("Bsgamma::Phi27_1_im(): z can not be greater than 1/4");
    
    double d = delta(E0);
    double z2 = z*z;
    double st0 = sqrt(1. - 4.*z);
    double std = sqrt( (1. - d - 4.*z) * (1. - d) );
    double L0 = log( ( 1. + st0 ) / ( 2.*sqrt(z) ) );
    double Ld = log( ( sqrt(1. - d) + sqrt(1. - d - 4.*z) ) / ( 2.*sqrt(z) ) );
    
    if (z < (1. - d)/4.)
        return 8./9. * M_PI * z * ( (1. - 4. * z + 6. * z2)* (L0-Ld) - 3./4. * (1. - 2. * z) * (st0-std) 
                + d * (2. - d - 4. * z) * Ld - 3./4. * d * std );
    else
        return 8./9. * M_PI * z * ( (1. - 4. * z + 6. * z2) * L0 - 3./4. * (1. - 2. * z) * st0 );
}

double Bsgamma::Phi28_1(double E0, double z)
{
    return -Phi27_1(E0, z)/3.;
}

double Bsgamma::Phi33_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return 2./27. * (Int_b1(E0) + 8.*Int_b2(E0) - 4.*Int_b4(E0)
            + 33.*Int_bb1(E0) - 20.*Int_bb2(E0) + 4.*Int_bb4(E0))
            + 1./18. * d * ( 1./2. - 1./2.*d2 + 1./3.*d3 - 1./15.*d4 );
}

double Bsgamma::Phi34_1(double E0)
{
    return -1./3.*Phi33_1(E0);
}

double Bsgamma::Phi35_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return 32./27. * (2.*Int_b1(E0) + 4.*Int_b2(E0) - 2.*Int_b4(E0)
            + 18.*Int_bb1(E0) - 13.*Int_bb2(E0) + 2.*Int_bb4(E0))
            + 4./9. * d * ( 4./3. - d2 + 1./2.*d3 - 1./15.*d4 );
}

double Bsgamma::Phi36_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return 8./81. * (5.*Int_b1(E0) + Int_b2(E0) + 4.*Int_b4(E0)
            - 18.*Int_bb1(E0) + 8.*Int_bb2(E0) - 4.*Int_bb4(E0))
            + 6. * (Phi23_1(E0) + Phi23_1_4body(E0))
            - 2./27. * d * ( 4./3. - d2 + 1./2.*d3 - 1./15.*d4 );
}

double Bsgamma::Phi37_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    
    return -4./3. * Int_b3(E0) + 1./9. * d * (1. - d + 1./3.*d2) + 1./4.*ff7_sMP(E0);
}

double Bsgamma::Phi38_1(double E0)
{
    return -1./3. * ( Phi37_1(E0) - 1./4.*ff7_sMP(E0) ) + 1./4.*ff8_sMP(E0);
}

double Bsgamma::Phi44_1(double E0)
{
    return 1./36. * Phi33_1(E0);
}

double Bsgamma::Phi45_1(double E0)
{
    return -1./6. * Phi35_1(E0);
}

double Bsgamma::Phi46_1(double E0)
{
    return -1./6. * Phi36_1(E0);
}

double Bsgamma::Phi47_1(double E0)
{
    return -1./6. * ( Phi37_1(E0) - 1./4.*ff7_sMP(E0) ) 
            + 1./4. * (-1./6. * ff7_sMP(E0) + ff7_dMP(E0));
}

double Bsgamma::Phi48_1(double E0)
{
    return -1./3. * (Phi47_1(E0) - 1./4. * (-1./6. * ff7_sMP(E0) + ff7_dMP(E0)) ) 
            + 1./4. * (-1./6. * ff8_sMP(E0) + ff8_dMP(E0));
}

double Bsgamma::Phi55_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return 128./27. * (4.*Int_b1(E0) + 2.*Int_b2(E0) - Int_b4(E0)
            + 12.*Int_bb1(E0) - 8.*Int_bb2(E0) + Int_bb4(E0))
            + 8./9. * d * ( 11./3. - 2.*d2 + 2./3.*d3 - 1./15.*d4 );
}

double Bsgamma::Phi56_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return 32./81. * (20.*Int_b1(E0) + Int_b2(E0) + 4.*Int_b4(E0)
            + 24.*Int_bb1(E0) + 14.*Int_bb2(E0) - 4.*Int_bb4(E0))
            + 6. * (Phi25_1(E0) + Phi25_1_4body(E0))
            - 8./27. * d * ( 11./3. - 2.*d2 + 2./3.*d3 - 1./15.*d4 );
}

double Bsgamma::Phi57_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    
    return 16./9. * d * ( 1. - d + 1./3.*d2) + 4. * ff7_sMP(E0);
}

double Bsgamma::Phi58_1(double E0)
{
    return -1./3. * (Phi57_1(E0) - 4. *  ff7_sMP(E0)) 
            + 4. * ff8_sMP(E0);
}

double Bsgamma::Phi66_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    double d4=d3*d;
    
    return -8./243. * (56.*Int_b1(E0) + 10.*Int_b2(E0) + 4.*Int_b4(E0)
            + 15.*Int_bb1(E0) - 4.*Int_bb2(E0) - 4.*Int_bb4(E0))
            + 32./27. * (4.*Int_c1(E0) + Int_c2(E0) - Int_bc1(E0) 
            - 2.*Int_bc2(E0) + 9.*Int_cc1(E0))
            + 2./81. * d * ( 11./3. - 2.*d2 + 2./3.*d3 - 1./15.*d4 );
}

double Bsgamma::Phi67_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    
    return 8./3. * (Int_b3(E0) - 2.*Int_c3(E0)) - 8./27. * d * ( 1. - d + 1./3.*d2) 
            + 1./4. * (-8./3. * ff7_sMP(E0) + 10. * ff7_dMP(E0));
}

double Bsgamma::Phi68_1(double E0)
{
    return -1./3. * (Phi67_1(E0) - 1./4. * (-8./3. * ff7_sMP(E0) + 10. * ff7_dMP(E0)) ) 
            + 1./4. * (-8./3. * ff8_sMP(E0) + 10. * ff8_dMP(E0));
}

double Bsgamma::Phi77_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    return -2./3.*pow(log(d),2.) - 7./3.*log(d) - 31./9. + 10./3.*d + d2/3. - 2./9.*d3 + d*(d - 4.)*log(d)/3.;
}

double Bsgamma::Phi78_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    double Li2 = gsl_sf_dilog(1. - d);
    
    double pi2=M_PI*M_PI;
    
    return 8./9.*( Li2 -  pi2/6. - d*log(d) + 9./4.*d - d2/4. + d3/12.);
}

double Bsgamma::Phi88_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    double d3=d2*d;
    
    double Li2 = gsl_sf_dilog(1. - d);
    
    double pi2=M_PI*M_PI;
    
    return 1./27.*( -2.*log(Mb_kin/Ms)*( d2 + 2.*d + 4.*log(1. - d) ) + 4.*Li2 - 2./3.*pi2 - d*(2. + d)*log(d) 
            + 8.*log(1. - d) -  2./3.*d3 + 3.*d2 +7*d);
}

double Bsgamma::Kij_1(int i, int j, double E0, double mu)
{
    if (i > 8 || j>8 || i<1 || j<1) throw std::runtime_error("Bsgamma::Kij_1(): indexes (i,j) must be included in (1,..,8)");
    
    int temp;
    
    if (i > j) {temp=i; i=j; j=temp;}
    
    double gamma_i7[8] = {-208./243., 416./81., -176./81., -152./243., -6272./81., 4624./243., 32./3., -32./9.};
    double K_ij[8][8];
    double Lb = log(mu/Mb_kin);
    
    K_ij[0][0] = 4.*Phi11_1(E0);
    K_ij[0][1] = 2.*Phi12_1(E0);
    K_ij[0][2] = 2.*Phi13_1(E0);
    K_ij[0][3] = 2.*Phi14_1(E0);
    K_ij[0][4] = 2.*Phi15_1(E0);
    K_ij[0][5] = 2.*Phi16_1(E0);
    K_ij[0][6] = r1(1,zeta()).real() - gamma_i7[0]*Lb + 2.*Phi17_1(E0, zeta());
    K_ij[0][7] = 2.*Phi18_1(E0, zeta());
    
    K_ij[1][1] = 4.*Phi22_1(E0);
    K_ij[1][2] = 2.*Phi23_1(E0);
    K_ij[1][3] = 2.*Phi24_1(E0);
    K_ij[1][4] = 2.*Phi25_1(E0);
    K_ij[1][5] = 2.*Phi26_1(E0);
    K_ij[1][6] = r1(2,zeta()).real() - gamma_i7[1]*Lb + 2.*Phi27_1(E0, zeta());
    K_ij[1][7] = 2.*Phi28_1(E0, zeta());
    
    K_ij[2][2] = 4.*Phi33_1(E0);
    K_ij[2][3] = 2.*Phi34_1(E0);
    K_ij[2][4] = 2.*Phi35_1(E0);
    K_ij[2][5] = 2.*Phi36_1(E0);
    K_ij[2][6] = r1(3,zeta()).real() - gamma_i7[2]*Lb + 2.*Phi37_1(E0);
    K_ij[2][7] = 2.*Phi38_1(E0);
    
    K_ij[3][3] = 4.*Phi44_1(E0);
    K_ij[3][4] = 2.*Phi45_1(E0);
    K_ij[3][5] = 2.*Phi46_1(E0);
    K_ij[3][6] = r1(4,zeta()).real() - gamma_i7[3]*Lb + 2.*Phi47_1(E0);
    K_ij[3][7] = 2.*Phi48_1(E0);
    
    K_ij[4][4] = 4.*Phi55_1(E0);
    K_ij[4][5] = 2.*Phi56_1(E0);
    K_ij[4][6] = r1(5,zeta()).real() - gamma_i7[4]*Lb + 2.*Phi57_1(E0);
    K_ij[4][7] = 2.*Phi58_1(E0);
    
    K_ij[5][5] = 4.*Phi66_1(E0);
    K_ij[5][6] = r1(6,zeta()).real() - gamma_i7[5]*Lb + 2.*Phi67_1(E0);
    K_ij[5][7] = 2.*Phi68_1(E0);
    
    K_ij[6][6] = -182./9. + 8./9.*M_PI*M_PI - gamma_i7[6]*2.*Lb + 4.*Phi77_1(E0);
    K_ij[6][7] =  r1(8,zeta()).real() - gamma_i7[7]*Lb + 2.*Phi78_1(E0);
    
    K_ij[7][7] = 4.*Phi88_1(E0);
    
    return K_ij[i-1][j-1];
}

void Bsgamma::computeCoeff(double mu)
{
    allcoeff = SM.getMyFlavour()->ComputeCoeffsgamma(mu);
    allcoeffprime = SM.getMyFlavour()->ComputeCoeffprimesgamma(mu);
    
    C1_0 = (*(allcoeff[LO]))(0);
    C2_0 = (*(allcoeff[LO]))(1);
    C3_0 = (*(allcoeff[LO]))(2);
    C4_0 = (*(allcoeff[LO]))(3);
    C5_0 = (*(allcoeff[LO]))(4);
    C6_0 = (*(allcoeff[LO]))(5);
    C7_0 = (*(allcoeff[LO]))(6);
    C8_0 = (*(allcoeff[LO]))(7);
    
    C1_1 = (*(allcoeff[NLO]))(0)/Alstilde;
    C2_1 = (*(allcoeff[NLO]))(1)/Alstilde;
    C3_1 = (*(allcoeff[NLO]))(2)/Alstilde;
    C4_1 = (*(allcoeff[NLO]))(3)/Alstilde;
    C5_1 = (*(allcoeff[NLO]))(4)/Alstilde;
    C6_1 = (*(allcoeff[NLO]))(5)/Alstilde;
    C7_1 = (*(allcoeff[NLO]))(6)/Alstilde;
    C8_1 = (*(allcoeff[NLO]))(7)/Alstilde;
    
    C7p_0 = (*(allcoeffprime[LO]))(6);
    C7p_1 = (*(allcoeffprime[NLO]))(6)/Alstilde;

}

double Bsgamma::P0(double E0)
{
    return C7_0.abs2() + C7p_0.abs2() + P0_4body(E0,Mb_kin*Mb_kin/Ms/Ms);
}

double Bsgamma::P11()
{
    return 2.*((C7_0*C7_1).real() + (C7p_0*C7p_1).real()); /*CHECK SIGN*/
}

double Bsgamma::P21(double E0, double mu)
{
    int i,j;
    gslpp::complex C0[8]={C1_0,C2_0,C3_0,C4_0,C5_0,C6_0,C7_0,C8_0};
    gslpp::complex C0p[8]={C7p_0}; /*IMPLEMENT OTHER WC*/
    double p21=0.;
    
    for(i=0;i<8;i++)
    {
        for(j=0;j<8;j++)
        {
            p21 += (C0[i]*C0[j]).real() * Kij_1(i+1,j+1,E0,mu);
        }
    }
    
    for(i=6;i<7;i++) /*CHECK ALGORITHM*/
    {
        for(j=6;j<7;j++)
        {
            p21 += (C0p[i]*C0p[j]).real() * Kij_1(i+1,j+1,E0,mu);
        }
    }
    
    return p21;
}

double Bsgamma::Vub_NLO_2body(bool CPodd)
{
    double z = zeta();
    
    return 4. * Alstilde * (C7_0 * ( C2_0 - C1_0/6. )).real() *
            (CKMu.real()*( a(z) + b(z) ).real() - CPodd * CKMu.imag()*( a(z) + b(z) ).imag());
}

double Bsgamma::Vub_NLO_3body(double E0,bool CPodd)
{
    double d = delta(E0);
    
    return 64./27. * Alstilde * ( C2_0 - C1_0/6. ).abs2() *
            ( CKMu.real() * ( 2. * Int_cc1(E0) - Int_c1(E0) )
            + CKMu.abs2() *  ( Int_cc1(E0) - Int_c1(E0) + 1./8. * d * ( 1. - 1./3. * d*d ) )
            - CPodd * CKMu.imag() * Int_c1_im(E0) )
            + 4. * Alstilde * (( C7_0 - C8_0/3. ) * ( C2_0 - C1_0/6. )).real() *
            ( CKMu.real() * ( Phi27_1(E0,zeta()) + 2./9. * d * ( 1. - d + 1./3. * d*d ) )
            - CPodd * CKMu.imag() * Phi27_1_im(E0,zeta()) );
}

double Bsgamma::Vub_NLO_4body(double E0, bool CPodd)
{
    if (FOUR_BODY) {
        double d = delta(E0);
        double d2 = d*d;
        double d3 = d2*d;
        double Ld = log(d);
        double Lumd = log(1. - d);
        double Lq = log(Ms/Mb_kin);

        double uphib427 = ( 2. * d * (-63. + 30. * d + 35. * d2 - 2. * d3 
                    + 3. * d * (-18. - 7. * d + d2) * Ld) ) / ( 243. * (d - 1.) );
        double uphib428 = ( 108. * (d - 1.) * (d - 1.) * Lumd*Lumd 
                    - 12. * Lumd * (- 25. - 18. * Lq - 18. * d * (5. + 4. * Lq) 
                    + 9. * d2 * (5. + 4. * Lq) + (9. + 36. * d - 18. * d2) * Ld) 
                    + d * (24. * (17. + 9. * Lq) + 27. * d * (43. + 26. * Lq) 
                    - d2 * (127. + 72. * Lq) + 9. * (-12. - 39. * d + 4. * d2) * Ld) 
                    + 108. * (-1. - 4. * d + 2. * d2) * gsl_sf_dilog(d) ) / 729.;

        return 4. * Alstilde * ( ( C2_0 - C1_0/6. ).abs2() * 
                ( CKMu.real() * 0.005025213076791178 + CPodd * CKMu.imag() * 0.013978889449487913)
                + ( C2_0 - C1_0/6. ).real() * CKMu.real() * (C7_0.real() * uphib427 + C8_0.real() * uphib428) );
    }
    
    else return 0.;
}

double Bsgamma::Vub_NLO(double E0, bool CPodd)
{
    return Vub_NLO_2body(CPodd) + Vub_NLO_3body(E0,CPodd) + Vub_NLO_4body(E0,CPodd);
}

double Bsgamma::P(double E0, double mu_b, double mu_c, orders order, bool CPodd)
{
    switch(order) {
        case NLO:
            return P0(E0) + Alstilde * (P11() + P21(E0,mu_b)) + Vub_NLO(E0, CPodd);
            break;
        case LO:
            return P0(E0);
            break;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("Bsgamma::P(): order " + out.str() + " not implemented");
    }
}

double Bsgamma::N_27()
{
    double mcnorm = 1.131; // value fixed according to arXiv:1003.5012, in order to employ the remaining corrections given in that work
    double lambda2 = mu_G2/3.;
    
    return -1./18. * (C7_0 * ( 2.*C2_0 - C1_0/3. )).real() * lambda2/mcnorm/mcnorm;
}

double Bsgamma::N_77(double E0, double mu)
{
    double z = 1. - delta(E0);
    double z2 = z*z;
    double z3 = z2*z;
    double z4 = z3*z;
    double umz2 = (1.-z)*(1.-z);
    double Lz = log(1. - z);
    double Lz2 = Lz*Lz;
    double Lb = 2. * log(mu/Mb_kin);
    
    double corrLambda2_rad;
    double corrLambda2_sem;
    double corrLambda2_mix;
    double corrLambda2;
    double corrLambda3;
    
    double alsb = SM.Alstilde5(Mb_kin);
    double Lambda_pert = 64./9. * alsb * mu_kin * 
                (1. + 4. * alsb * (9./2. * (log(Mb_kin/2./mu_kin) + 8./3.)
                - 3. * (M_PI*M_PI/6. - 13./12.)) );
    double mu_pi2_pert = 3./4. * mu_kin * Lambda_pert - 48. * alsb*alsb * mu_kin*mu_kin;
    double rho_D3_pert = 1./2. * mu_kin*mu_kin * Lambda_pert - 128./3. * alsb*alsb * mu_kin*mu_kin*mu_kin;
    
    double lambda1 = -mu_pi2 + mu_pi2_pert;
    double lambda2 = mu_G2/3.;
    double rho1 = rho_D3 - rho_D3_pert;
    
    double f1EGN = 16./9. * ( 4. - M_PI*M_PI) - 8./3. * Lz2 - 
             ( 4. * z * ( 30. - 63. * z + 31. * z2 + 5. * z3))/(9. * umz2) -
             ( 4. * (30. - 72. * z + 51. * z2 - 2. * z3 - 3. * z4))/(9. * umz2) * Lz;
    double f2EGN = -2./9. * ( 87. + 32. * M_PI*M_PI) - 32./3. * Lz2 + 
             2. * ( 162. - 244. * z + 113. * z2 - 7. * z3)/(3. * (1. - z)) * Lz +
             2. * z * ( 54 - 49. * z + 15. * z2)/(1. - z);
    
    corrLambda2_rad = lambda1 * ( f1EGN/8.  - 4./3. * (Lb + 1.) ) 
            + lambda2 * (f2EGN/8.  + 12. * (Lb + 1.) );
    corrLambda2_sem = -3. * 4.98 * lambda2 + (25. - 4. * M_PI*M_PI)/12.*lambda1;
    corrLambda2_mix = 1./8. * (9. * lambda2 - lambda1) * Kij_1(7,7,E0,mu);
    
    corrLambda2 = corrLambda2_rad - corrLambda2_sem + corrLambda2_mix;
    
    corrLambda3 = (-88./6. + 16.*log(2.))* rho1 /Mb_kin/Mb_kin/Mb_kin;
    
    return (C7_0.abs2() + C7p_0.abs2()) * (4. * Alstilde / Mb_kin / Mb_kin * corrLambda2 + corrLambda3);
}

double Bsgamma::N(double E0, double mu)
{
    return N_27() + N_77(E0,mu) + BLNPcorr * P0(E0);
}

double Bsgamma::C_sem()
{
    double z=zeta();
    return (1. - 8. * z + 8. * z*z*z - z*z*z*z - 12. * z*z * log(z)) * ( 0.903 
            - 0.588 * (SM.Alstilde5(4.6)*4*M_PI - 0.22) + 0.0650 * (Mb_kin - 4.55) 
            - 0.1080 * (Mc - 1.05) - 0.0122  * mu_G2 - 0.199 * rho_D3 + 0.004 * rho_LS3);
}

void Bsgamma::updateParameters()
{
    mu_kin=SM.getGambino_mukin();
    BRsl=SM.getGambino_BRsem()/100.;
    Mb_kin=SM.getGambino_Mbkin();
    Mc=SM.getGambino_Mcatmuc();
    mu_pi2=SM.getGambino_mupi2();
    rho_D3=SM.getGambino_rhoD3();
    mu_G2=SM.getGambino_muG2();
    rho_LS3=SM.getGambino_rhoLS3();
    C=C_sem();
    
    ale=SM.getAle();
    E0=SM.getbsgamma_E0();
    mu_b=SM.getMub();
    mu_c=SM.getMuc();
    alsUps=8./M_PI * mu_kin/Mb_kin * ( 1. + 3./8. * mu_kin/Mb_kin );
    Alstilde = SM.Alstilde5(mu_b);
    Ms=SM.getQuarks(QCD::STRANGE).getMass();
    lambda_t=SM.computelamt_s();
    V_cb=SM.getCKM().getVcb();
    CKMu=SM.computelamu_s().conjugate() / SM.computelamt_s().conjugate();
    
    BLNPcorr=SM.getBLNPcorr();
    
    checkCache();
    
    if (Intb_updated == 0) {
        Intb1Cached = 0;
        Intb2Cached = 0;
        Intb3Cached = 0;
        Intb4Cached = 0;
        Intbb1Cached = 0;
        Intbb2Cached = 0;
        Intbb4Cached = 0;
    }
    if (Intbc_updated == 0) {
        Intbc1Cached = 0;
        Intbc2Cached = 0;
        Intc1Cached = 0;
        Intc1imCached = 0;
        Intc2Cached = 0;
        Intc3Cached = 0;
        IntccCached = 0;
        Intcc1Cached = 0;
        Intcc1p1Cached = 0;
    }
    
    computeCoeff(mu_b);
    
    overall = BRsl * (lambda_t/V_cb).abs2() * 6. * ale / (M_PI * C);
}

double Bsgamma::computeThValue()
{
    updateParameters();
    
    if (obs == 1) 
        return overall *  ( P(E0, mu_b, mu_c, NLO, false) + N(E0,mu_b) );
    if (obs == 2) 
        return overall *  ( P(E0, mu_b, mu_c, NLO, true) + N(E0,mu_b) );
    
    throw std::runtime_error("Bsgamma::computeThValue(): Observable type not defined. Can be only 1 or 2");
}