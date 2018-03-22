/* 
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

/* 
 * phi1.4body partially hardcoded, currently switched off by setting FOUR_BODY to false
 * EW matching missing
 */

#include "StandardModel.h"
#include "bsgamma.h"
#include "std_make_vector.h"
#include "gslpp_function_adapter.h"
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_clausen.h>
#include <boost/bind.hpp>

Bsgamma::Bsgamma(const StandardModel& SM_i, QCD::quark quark_i, int obsFlag)
: ThObservable(SM_i),
Intbc_cache(2, 0.)
{    
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR) or 2 (ACP)");
    
    quark = quark_i;
    
    SUM = false;
    EWflag = true;
    FOUR_BODY = false;
    WET_NP_btos = false;
    SMEFT_NP_btos = false;
    
    setParametersForObservable(make_vector<std::string>() << "Gambino_mukin" << "Gambino_BRsem" << "Gambino_Mbkin" << "Gambino_Mcatmuc" << "Gambino_mupi2" 
                                                          << "Gambino_rhoD3" << "Gambino_muG2" << "Gambino_rhoLS3" << "BLNPcorr" << "mu_b_bsgamma" << "mu_c_bsgamma");
    
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
    IntPhi772rCached = 0;
    
    w_INT = gsl_integration_cquad_workspace_alloc (100);
}

Bsgamma::Bsgamma(const StandardModel& SM_i, int obsFlag)
: ThObservable(SM_i),
Intbc_cache(2, 0.)
{
    if (obsFlag > 0 and obsFlag < 3) obs = obsFlag;
    else throw std::runtime_error("obsFlag in bsgamma can only be 1 (BR) or 2 (ACP)");
    
    SUM = true;
    EWflag = true;
    FOUR_BODY = false;
    
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
    IntPhi772rCached = 0;
    
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
    gslpp::complex la_u =-CKMu;
    
    double A1sq =C1_0.abs2()*CKMusq;
    double A2sq =C2_0.abs2()*CKMusq;
    
    double C13re = (C1_0.real()*C3_0.real() + C1_0.imag()*C3_0.imag());
    double C14re = (C1_0.real()*C4_0.real() + C1_0.imag()*C4_0.imag());
    double C15re = (C1_0.real()*C5_0.real() + C1_0.imag()*C5_0.imag());
    double C16re = (C1_0.real()*C6_0.real() + C1_0.imag()*C6_0.imag());
    
    double C13im = (C1_0.real()*C3_0.imag() + C1_0.imag()*C3_0.real());
    double C14im = (C1_0.real()*C4_0.imag() + C1_0.imag()*C4_0.real());
    double C15im = (C1_0.real()*C5_0.imag() + C1_0.imag()*C5_0.real());
    double C16im = (C1_0.real()*C6_0.imag() + C1_0.imag()*C6_0.real());
    
    double C23re = (C2_0.real()*C3_0.real() + C2_0.imag()*C3_0.imag());
    double C24re = (C2_0.real()*C4_0.real() + C2_0.imag()*C4_0.imag());
    double C25re = (C2_0.real()*C5_0.real() + C2_0.imag()*C5_0.imag());
    double C26re = (C2_0.real()*C6_0.real() + C2_0.imag()*C6_0.imag());
    
    double C23im = (C2_0.real()*C3_0.imag() + C2_0.imag()*C3_0.real());
    double C24im = (C2_0.real()*C4_0.imag() + C2_0.imag()*C4_0.real());
    double C25im = (C2_0.real()*C5_0.imag() + C2_0.imag()*C5_0.real());
    double C26im = (C2_0.real()*C6_0.imag() + C2_0.imag()*C6_0.real());
    
    double C13 = (C13re*la_u.real() - C13im*la_u.imag());
    double C14 = (C14re*la_u.real() - C14im*la_u.imag());
    double C15 = (C15re*la_u.real() - C15im*la_u.imag());
    double C16 = (C16re*la_u.real() - C16im*la_u.imag());
    double C23 = (C23re*la_u.real() - C23im*la_u.imag());
    double C24 = (C24re*la_u.real() - C24im*la_u.imag());
    double C25 = (C25re*la_u.real() - C25im*la_u.imag());
    double C26 = (C26re*la_u.real() - C26im*la_u.imag());
    double C33 = C3_0.abs2();
    double C34 = (C3_0.real()*C4_0.real() + C3_0.imag()*C4_0.imag());
    double C35 = (C3_0.real()*C5_0.real() + C3_0.imag()*C5_0.imag());
    double C36 = (C3_0.real()*C6_0.real() + C3_0.imag()*C6_0.imag());
    double C44 = C4_0.abs2();
    double C45 = (C4_0.real()*C5_0.real() + C4_0.imag()*C5_0.imag());
    double C46 = (C4_0.real()*C6_0.real() + C4_0.imag()*C6_0.imag());
    double C55 = C5_0.abs2();
    double C56 = (C5_0.real()*C6_0.real() + C5_0.imag()*C6_0.imag());
    double C66 = C6_0.abs2();
    
    return (C33 + 20. * C35 + 2./9. * C44 + 40./9. * C46 + 136. * C55 + 272./9. * C66) * T1(E0,t) +
            
            (2./9. * A1sq + A2sq 
            + 8./9. * C13 - 4./27. * C14 + 128./9. * C15 - 64./27. * C16
            + 2./3. * C23 + 8./9. * C24 + 32./3. * C25 + 128./9. * C26) * T2(E0,t) +
            
            (C33 + 8./3. * C34 + 32. * C35 + 128./3. * C36 - 2./9. * C44 + 128./3. * C45 
            - 64./9. * C46 + 256. * C55 + 2048./3 * C56 - 512./9. * C66) * T3(E0,t);
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

gslpp::complex Bsgamma::r1_ew(int i, double z)
{
    double Xb = -0.16844083981858157;
    double PI2 = M_PI*M_PI;
    gslpp::complex iPI = gslpp::complex::i()*M_PI;
    
    switch(i){
        case 1:
            return 3332./2187. - 4.*(a(z) + b(z))/9. + 160./729.*iPI;
        case 2:
            return 833./729. - (a(z) + b(z))/3. + 40./243.*iPI;
        case 3:
            return 748./729. + 2.*M_PI/9./sqrt(3.) - 2./81.*PI2 + 8./27.*Xb 
                    - a(1.)/12. + 7.*b(1.)/6. - 2.*b(z) + 26./243.*iPI;
        case 4:
            return 2680./2187. + 8.*M_PI/27./sqrt(3.) - 8./243.*PI2 + 32./81.*Xb 
                    - a(1.)/9. + 2.*b(1.)/9. + 56./729.*iPI;
        case 5:
            return 78301./729. + 8.*M_PI/9./sqrt(3.) - 40./81.*PI2 + 32./27.*Xb 
                    - 13.*a(1.)/3. + 38.*b(1.)/3. - 12.*a(z) - 20.*b(z) + 3908./243.*iPI;
        case 6:
            return 62440./2187. + 32.*M_PI/27./sqrt(3.) - 160./243.*PI2 + 128./81.*Xb 
                    - 16.*a(1.)/9. + 32.*b(1.)/9. + 896./729.*iPI;
        case 7:
            return -25./27. - 2./9.*iPI;
        default:
            std::stringstream out;
            out << i;
            throw std::runtime_error("Bsgamma::r1_ew(): index " + out.str() + " not implemented");
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

gslpp::complex Bsgamma::Int_bc1(double E0)
{
    if (Intbc1Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt = avaINT;
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_Kb_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt += gslpp::complex::i() * avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt2 = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_Kb_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt2 += gslpp::complex::i() * avaINT;

        CacheIntbc1 = delta(E0)*mt + mt2;
        Intbc1Cached = 1;
    }

    return CacheIntbc1;
}

gslpp::complex Bsgamma::Int_bc2(double E0)
{
    if (Intbc2Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt = avaINT;
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_Kb_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt += gslpp::complex::i() * avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_Kb_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt2 = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_Kb_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt2 += gslpp::complex::i() * avaINT;

        CacheIntbc2 = delta(E0)*mt + mt2;
        Intbc2Cached = 1;
    }

    return CacheIntbc2;
}

gslpp::complex Bsgamma::Int_c1(double E0)
{
    if (Intc1Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt = avaINT;
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt += gslpp::complex::i() * avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt2 = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt2 += gslpp::complex::i() * avaINT;

        CacheIntc1 = delta(E0)*mt + mt2;
        Intc1Cached = 1;
    }

    return CacheIntc1;
}

gslpp::complex Bsgamma::Int_c2(double E0)
{
    if (Intc2Cached == 0) {
        double t1 = (1. - delta(E0));

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt += gslpp::complex::i() * avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt2 = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_t_1mt2, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt2 += gslpp::complex::i() * avaINT;

        CacheIntc2 = delta(E0)*mt + mt2;
        Intc2Cached = 1;
    }

    return CacheIntc2;
}

gslpp::complex Bsgamma::Int_c3(double E0)
{
    if (Intc3Cached == 0) {
        double t1 = (1. - delta(E0));
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex t = avaINT;
    
        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_t, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        t += gslpp::complex::i() * avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_re_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        gslpp::complex mt = avaINT;

        INT = convertToGslFunction(boost::bind(&Bsgamma::getKc_im_t_1mt, &(*this), _1));
        if (gsl_integration_cquad(&INT, t1, 1., 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        mt += gslpp::complex::i() * avaINT;
        
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

gslpp::complex Bsgamma::Phi13_1(double E0)
{
    return -Phi23_1(E0)/6.;
}

gslpp::complex Bsgamma::Phi14_1(double E0)
{
    return -Phi24_1(E0)/6.;
}

gslpp::complex Bsgamma::Phi15_1(double E0)
{
    return -Phi25_1(E0)/6.;
}

gslpp::complex Bsgamma::Phi16_1(double E0)
{
    return -Phi26_1(E0)/6.;
}

gslpp::complex Bsgamma::Phi17_1(double E0, double z)
{
    return -Phi27_1(E0,z)/6.;
}

gslpp::complex Bsgamma::Phi18_1(double E0, double z)
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

gslpp::complex Bsgamma::Phi23_1(double E0)
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

gslpp::complex Bsgamma::Phi24_1(double E0)
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

gslpp::complex Bsgamma::Phi25_1(double E0)
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

gslpp::complex Bsgamma::Phi26_1(double E0)
{
    return 16./81. * (4.*Int_c1(E0) + Int_c2(E0) - 10.*Int_bc1(E0) - 2.*Int_bc2(E0) + 36.*Int_cc1(E0))
            - Phi26_1_4body(E0);
}

gslpp::complex Bsgamma::Phi27_1(double E0, double z)
{
    double d = delta(E0);
    double d2 = d*d;
    double z2 = z*z;
    double Pi2 = M_PI*M_PI;
    double st0 = sqrt(1. - 4.*z);
    double std = sqrt( (1. - d - 4.*z) * (1. - d) );
    double L0 = log( ( 1. + st0 ) / ( 2.*sqrt(z) ) );
    double Ld = log( ( sqrt(1. - d) + sqrt(1. - d - 4.*z) ) / ( 2.*sqrt(z) ) );
    
    gslpp::complex res;
    
    if (d == 1) {
        res = -2./27. + (2.*Pi2 - 7.)/9. * z + 4.*(3. - 2.*Pi2)/9. * z * z
                + 4./3. * z * (1. - 2.*z) * st0 * L0
                - 8./9. * z * (6.*z*z - 4.*z + 1.) * L0*L0 + 4./3. * Pi2 * z * z *z;
    } else res = -2./27. * d * (3. - 3.*d + d2) + (2.*Pi2 - 7.)/9. * z * d * (2. - d)
            + 4.*(3. - 2.*Pi2)/9. * z * z * d 
            + 4./3. * z * (1. - 2.*z) * ( st0 * L0 - std * Ld ) 
            + 4./3. * z * d * std * Ld 
            - 8./9. * z * (6.*z*z - 4.*z + 1.) * ( L0*L0 - Ld*Ld ) 
            - 8./9. * z * d * (2. - d - 4.*z) * Ld * Ld;
    
    if (z < (1. - d)/4.)
        res += gslpp::complex::i() * 8./9. * M_PI * z * ( (1. - 4. * z + 6. * z2)* (L0-Ld) - 3./4. * (1. - 2. * z) * (st0-std) 
                + d * (2. - d - 4. * z) * Ld - 3./4. * d * std );
    else
        res += gslpp::complex::i() * 8./9. * M_PI * z * ( (1. - 4. * z + 6. * z2) * L0 - 3./4. * (1. - 2. * z) * st0 );
    
    return res;
}

gslpp::complex Bsgamma::Phi28_1(double E0, double z)
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

gslpp::complex Bsgamma::Phi36_1(double E0)
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

gslpp::complex Bsgamma::Phi46_1(double E0)
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

gslpp::complex Bsgamma::Phi56_1(double E0)
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

gslpp::complex Bsgamma::Phi66_1(double E0)
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

gslpp::complex Bsgamma::Phi67_1(double E0)
{
    double d=delta(E0);
    double d2=d*d;
    
    return 8./3. * (Int_b3(E0) - 2.*Int_c3(E0)) - 8./27. * d * ( 1. - d + 1./3.*d2) 
            + 1./4. * (-8./3. * ff7_sMP(E0) + 10. * ff7_dMP(E0));
}

gslpp::complex Bsgamma::Phi68_1(double E0)
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

gslpp::complex Bsgamma::Kij_1(int i, int j, double E0, double mu)
{
    if (i > 8 || j>8 || i<1 || j<1) throw std::runtime_error("Bsgamma::Kij_1(): indices (i,j) must be included in (1,..,8)");
    
    double gamma_i7[8] = {-208./243., 416./81., -176./81., -152./243., -6272./81., 4624./243., 32./3., -32./9.};
    gslpp::complex K_ij[8][8] = {{0.}};
    double Lb = log(mu/Mb_kin);
    
    K_ij[0][0] = 4.*Phi11_1(E0);
    K_ij[0][1] = 2.*Phi12_1(E0);
    K_ij[0][2] = 2.*Phi13_1(E0);
    K_ij[0][3] = 2.*Phi14_1(E0);
    K_ij[0][4] = 2.*Phi15_1(E0);
    K_ij[0][5] = 2.*Phi16_1(E0);
    K_ij[0][6] = r1(1,zeta()) - gamma_i7[0]*Lb + 2.*Phi17_1(E0, zeta());
    K_ij[0][7] = 2.*Phi18_1(E0, zeta());
    
    K_ij[1][1] = 4.*Phi22_1(E0);
    K_ij[1][2] = 2.*Phi23_1(E0);
    K_ij[1][3] = 2.*Phi24_1(E0);
    K_ij[1][4] = 2.*Phi25_1(E0);
    K_ij[1][5] = 2.*Phi26_1(E0);
    K_ij[1][6] = r1(2,zeta()) - gamma_i7[1]*Lb + 2.*Phi27_1(E0, zeta());
    K_ij[1][7] = 2.*Phi28_1(E0, zeta());
    
    K_ij[2][2] = 4.*Phi33_1(E0);
    K_ij[2][3] = 2.*Phi34_1(E0);
    K_ij[2][4] = 2.*Phi35_1(E0);
    K_ij[2][5] = 2.*Phi36_1(E0);
    K_ij[2][6] = r1(3,zeta()) - gamma_i7[2]*Lb + 2.*Phi37_1(E0);
    K_ij[2][7] = 2.*Phi38_1(E0);
    
    K_ij[3][3] = 4.*Phi44_1(E0);
    K_ij[3][4] = 2.*Phi45_1(E0);
    K_ij[3][5] = 2.*Phi46_1(E0);
    K_ij[3][6] = r1(4,zeta()) - gamma_i7[3]*Lb + 2.*Phi47_1(E0);
    K_ij[3][7] = 2.*Phi48_1(E0);
    
    K_ij[4][4] = 4.*Phi55_1(E0);
    K_ij[4][5] = 2.*Phi56_1(E0);
    K_ij[4][6] = r1(5,zeta()) - gamma_i7[4]*Lb + 2.*Phi57_1(E0);
    K_ij[4][7] = 2.*Phi58_1(E0);
    
    K_ij[5][5] = 4.*Phi66_1(E0);
    K_ij[5][6] = r1(6,zeta()) - gamma_i7[5]*Lb + 2.*Phi67_1(E0);
    K_ij[5][7] = 2.*Phi68_1(E0);
    
    K_ij[6][6] = -182./9. + 8./9.*M_PI*M_PI - gamma_i7[6]*2.*Lb + 4.*Phi77_1(E0);
    K_ij[6][7] =  r1(8,zeta()) - gamma_i7[7]*Lb + 2.*Phi78_1(E0);
    
    K_ij[7][7] = 4.*Phi88_1(E0);
    
    if (j >= i ) return K_ij[i-1][j-1];
    else return K_ij[j-1][i-1].conjugate();
}

double Bsgamma::Rer22(double z)
{
    double L = log(z);
    double L2 = L*L;
    double L3 = L2*L;
    double L4 = L3*L;
    double z32 = sqrt(z)*z;
    double z2 = z*z;
    double z52 = z32*z;
    double z3 = z2*z;
    double z72 = z52*z;
    double z4 = z3*z;
    double Pi2 = M_PI*M_PI;
    double zeta3 = gsl_sf_zeta_int(3);
    
    return 67454./6561. - 124./729. * Pi2 
            - 4./1215. * (11280. - 1520. * Pi2 - 171. * Pi2*Pi2 - 5760. * zeta3 
            + 6840. * L - 1440. * Pi2*L - 2520. * zeta3*L 
            + 120. * L2 + 100. * L3 - 30. * L4) * z 
            - 64./243. * Pi2*( 43. - 12. * log(2.) - 3. * L) * z32 
            - 2./1215. * (11475. - 380. * Pi2 + 96. * Pi2*Pi2 
            + 7200. * zeta3 - 1110. * L - 1560. * Pi2*L + 1440. * zeta3*L 
            + 990. * L2 + 260. * L3 - 60. * L4) * z2 
            + 2240./243. * Pi2 * z52
            - 2./2187. * (62471. - 2424. * Pi2 - 33264. * zeta3 - 19494. * L 
            - 504. * Pi2*L - 5184. * L2 + 2160. * L3) * z3 
            - 2464./6075. * Pi2 * z72
            + ( - 15103841./546750. + 7912./3645. * Pi2 + 2368./81. * zeta3 
            + 147038./6075. * L + 352./243. * Pi2*L + 88./243. * L2 
            - 512./243. * L3 ) * z4; 
}

double Bsgamma::Phi22_2beta0(double E0, double mu)
{
    double Lb = 2*log(mu/Mb_kin);
    double d = delta(E0);
    double d2 = d*d;
    double mcmb = Mc/Mb_kin;
    double mcmb2 = mcmb*mcmb;
    double mcmb3 = mcmb2*mcmb;
    
    return SM.Beta0(5) * (Phi22_1(E0)*Lb
            + 0.013698269459646965 + 0.3356948452887703 * d 
            - 0.086677232161681  * d2 
            + ( 0.3575455009710223 + 1.8248223618702617 * d 
            - 0.374324331239819 * d2 ) * mcmb 
            + (-2.3059130759599302 - 5.799640881350228 * d 
            - 6.226247001127346 * d2 ) * mcmb2 
            + ( 3.4485885608332834 - 0.5479757965141787 * d 
            + 17.272487170738795 * d2 ) * mcmb3);
}

double Bsgamma::Phi28_2beta0(double E0, double mu)
{
    double Lb = 2*log(mu/Mb_kin);
    double d = delta(E0);
    double d2 = d*d;
    double mcmb = Mc/Mb_kin;
    double mcmb2 = mcmb*mcmb;
    double mcmb3 = mcmb2*mcmb;
    double mcmb4 = mcmb3*mcmb;
    double mcmb5 = mcmb4*mcmb;
    
    return SM.Beta0(5) * (Phi28_1(E0, zeta()).real()*Lb
            + 0.026054745293391798 + 0.1678721564514209 * d 
            - 0.19700988587274693 * d2 
            + ( -0.03801105485376407 + 0.601712887338462 * d 
            - 0.7557529126506585 * d2 ) * mcmb 
            + ( 2.7551159092192132 - 10.034450524236696 * d 
            + 11.271837772655209 * d2 ) * mcmb2 
            + ( -27.045289848315868 + 68.46851531490181 * d 
            - 72.50921751760909 * d2 ) * mcmb3 
            + ( 85.86574743951778 - 289.3441408351491 * d 
            + 297.6777008484198 * d2 ) * mcmb4 
            + ( -91.5260435658921 + 399.81982774456964 * d 
            - 399.85440571662446 * d2 ) * mcmb5);
}

double Bsgamma::Phi77_2beta0(double E0, double mu)
{
    double Lb = 2*log(mu/Mb_kin);
    double d = delta(E0);
    double d2 = d*d;
    double d3 = d2*d;
    double Ld = log(d);
    double zeta3 = gsl_sf_zeta_int(3);
    double Li2 = gsl_sf_dilog(1. - d);
    Polylogarithms Poly;
    double Li3 = Poly.Li3(d);
    
    return SM.Beta0(5) * (Phi77_1(E0)*Lb
            + ( -3. + 4./3. * d - 1./3. * d2 - 4./3. * Ld ) * Li2 
            + ( 13./18. + 2. * d - 1./2. * d2 - 4./3. * log(1. - d) 
            + 2./3. * Ld ) * Ld*Ld 
            - 8./3. * (Li3 - zeta3) 
            + ( 4./9. * M_PI*M_PI - 85./18. - 47./9. * d 
            + 19./18. * d2 + 2./9. * d3 ) * Ld 
            - 49./6. + 80./9. * d + 1./18. * d2 - 7./9. * d3);
}

double Bsgamma::Phi88_2beta0(double E0, double mu)
{
    double Lb = 2*log(mu/Mb_kin);
    double d = delta(E0);
    double d2 = d*d;
    double d3 = d2*d;
    double Ld = log(d);
    double L1d = log(1. - d);
    double Li2 = gsl_sf_dilog(1. - d);
    Polylogarithms Poly;
    double Li3 = Poly.Li3(d);
    
    return SM.Beta0(5) * (Phi88_1(E0)*Lb
            + 4./27. * ( - 2. * ( Li2 - 1./6. * M_PI*M_PI + 3. * L1d 
            - 1./4. * d * (2. + d) * Ld + 8./3. * d + 5./6. * d2 
            - 1./18. * d3 ) * log(Mb_kin/Ms) - 2. * Li3 
            + ( 5. - 2. * Ld ) * ( Li2 - 1./6. * M_PI*M_PI ) 
            + ( 1./2. * d + 1./4. * d2 - L1d ) * Ld * Ld
            - 1./12. * M_PI*M_PI * d * (2. + d) 
            + ( 151./18. - 1./3. * M_PI*M_PI ) * L1d 
            + ( - 53./12. * d - 19./12. * d2 + 2./9. * d3 ) * Ld 
            + 787./72. * d + 227./72. * d2 - 41./72. * d3 ));
}

double Bsgamma::dY1(double E0)
{
    double z0 = 1. - delta(E0);
    double Li2 = gsl_sf_dilog(z0);
    
    return + 2./9. * z0*(z0*z0 + 24.)- 8./3. * (z0 - 1.)*log(1. - z0) - 8./3. * Li2;
}

double Bsgamma::Y1(double E0, double mu)
{
    double Lb = log(mu/Mb_kin);
    
    return 4./9.*(29. - 2.*M_PI*M_PI) + 16./3.*Lb - dY1(E0);
}

double Bsgamma::Y2CF(double E0, double mu)
{
    double Lb = log(mu/Mb_kin);
    double z0 = 1. - delta(E0);
    double z02 = z0*z0;
    double z03 = z02*z0;
    double z04 = z03*z0;
    double z05 = z04*z0;
    double z06 = z05*z0;
    double z07 = z06*z0;
    double z08 = z07*z0;
    double z09 = z08*z0;
    double z010 = z09*z0;
    double z011 = z010*z0;
    double Lz = log(1. - z0);
    double Li2 = gsl_sf_dilog(z0);
    
    return -21.9087 - 112.464 * Lb - 42.6667 * Lb*Lb - 77.7675 * z0 + 10.6667 * Lb*z0 
            + 68.5848 * z02 - 5.33333 * Lb * z02 - 4.42133 * z03 + 6.22222 * Lb * z03 
            - 4.0317 * z04 + 6.64376 * z05 - 11.647 * z06 + 15.8697 * z07 
            - 14.8006 * z08 + 8.85514 * z09 - 2.9929 * z010 + 0.433254 * z011 
            + (-77.7675 + 85.8251 * z0 - 28.6855 * z02 
            + Lb * (-21.3333 - 21.3333 * z0 + 5.33333 * z02)) * Lz 
            + (-12.2881 - 10.6667 * Lb + 6.12213 * z0 + 0.27227 * z02) * Lz*Lz 
            + (-2.88573 + 5.77146 * z0 - 2.88573 * z02) * Lz*Lz*Lz - 32. * Lb*Li2;
}

double Bsgamma::Y2CA(double E0, double mu)
{
    double Lb = log(mu/Mb_kin);
    double z0 = 1. - delta(E0);
    double z02 = z0*z0;
    double z03 = z02*z0;
    double z04 = z03*z0;
    double z05 = z04*z0;
    double z06 = z05*z0;
    double z07 = z06*z0;
    double z08 = z07*z0;
    double z09 = z08*z0;
    double z010 = z09*z0;
    double z011 = z010*z0;
    double Lz = log(1. - z0);
    double Li2 = gsl_sf_dilog(z0);
    
    return 22.8959 + 76.5729 * Lb + 30.2222 * Lb*Lb + 2.94616 * z0 - 60.4444 * Lb*z0 
            - 13.2522 * z02 - 6.96907 * z03 - 2.51852 * Lb*z03 + 0.117907 * z04 
            - 2.02988 * z05 + 2.90402 * z06 - 3.53904 * z07 + 2.55728 * z08 
            - 0.941549 * z09 + 0.0173599 * z010 + 0.0598012 * z011 
            + (2.94616 - 30.2222 * Lb- 12.3947 * z0 + 30.2222 * Lb*z0 + 9.44855 * z02) * Lz 
            - 6.61587 * (1. - z0) * (1. - z0) * Lz*Lz + 30.2222 * Lb*Li2
            + (0.69995 - 1.3999 * z0 + 0.69995 * z02) * Lz*Lz*Lz;
}

double Bsgamma::Y2NL(double E0, double mu)
{
    double z0 = 1. - delta(E0);
    double Lz = log(1. - z0);
    double Lb = log(mu/Mb_kin);
    double zeta3 = gsl_sf_zeta_int(3);
    double Li2 = gsl_sf_dilog(z0);
    Polylogarithms Poly;
    double Li3 = Poly.Li3(z0);
    double Li3min = Poly.Li3(1. - z0);
    
    return -16./81. * (328. - 13.*M_PI*M_PI) - 64./27. * (18. - M_PI*M_PI)*Lb
            -64./9. * Lb*Lb + 64./3. * zeta3
            +4./27. * z0*(7.*z0*z0 - 17.*z0 + 238.) + 8./3. * Lb*dY1(E0)
            -8./27. * (z0*z0*z0 - 6.*z0*z0 + 80.*z0 - 75. + 6.*M_PI*M_PI)*Lz
            +16./3. * (z0 - 1.)*Lz*Lz + 16./3. * log(z0)*Lz*Lz
            +32./27. * (3.*z0 - 8.)*Li2 + 32./3. * Lz*Li2
            -32./9. * Li3 + 32./3. * Li3min - 32./3. * zeta3;
}

double Bsgamma::Y2NV_PHI1(double rho)
{
    double y = (1. - sqrt(1. - 4.*rho)) / (1. + sqrt(1. - 4.*rho));
    
    if (rho < 1./4.) 
        return log(y)*log(y) - M_PI*M_PI;
    else 
        return - acos( 1. - 1./(2. * rho) ) * acos( 1. - 1./(2. * rho) );
}


double Bsgamma::Y2NV_PHI2(double rho)
{
    double y = (1. - sqrt(1. - 4.*rho)) / (1. + sqrt(1. - 4.*rho));
    
    if (rho < 1./4.) 
        return log(y) * sqrt(1. - 4.*rho);
    else 
        return - acos( 1. - 1./(2. * rho) ) * sqrt(4.*rho - 1.);
}


double Bsgamma::Y2NV_PHI3(double rho)
{
    double y = (1. - sqrt(1. - 4.*rho)) / (1. + sqrt(1. - 4.*rho));
    
    if (rho < 1./4.) 
        return ( gsl_sf_dilog(-y) + log(y)*log(y)/4. + M_PI*M_PI/12. ) * sqrt(1. - 4.*rho);
    else 
        return - gsl_sf_clausen(2. * asin( 1./(2. * sqrt(rho)) )) * sqrt(4.*rho - 1.);
}


double Bsgamma::Y2NV_PHI4(double rho)
{
    double y = (1. - sqrt(1. - 4.*rho)) / (1. + sqrt(1. - 4.*rho));
    Polylogarithms Poly;
    ClausenFunctions Clausen;
    
    if (rho < 1./4.) 
        return Poly.Li3(- y) + log(y)*log(y)*log(y)/12. + log(y)*M_PI*M_PI/12.;
    else 
        return Clausen.Cl3(2. * asin( 1./(2. * sqrt(rho)) ));
}

double Bsgamma::Y2NV(double E0, double mu)
{
    double Lb = log(mu/Mb_kin);
    double rho = zeta();
    double rho2 = rho*rho;
    double rho32 = rho*sqrt(rho);
    double Lr = log(rho);
    double Li2 = gsl_sf_dilog(1. - rho);
    double Li2sqrt = gsl_sf_dilog(1. - sqrt(rho));
    Polylogarithms Poly;
    double Li3 = Poly.Li3(1. - rho);
    double Li3ov = Poly.Li3(1. - 1./rho);
    
    return -16./81. * (157. - 279.*rho - M_PI*M_PI*(5. + 9.*rho2 - 42.*rho32))
            -64./27. * (18. - M_PI*M_PI)*Lb - 64./9. * Lb*Lb
            +16./27. * (22. - M_PI*M_PI + 10.*rho)*Lr + 16./27. * (8. + 9.*rho2)*Lr*Lr
            -16./27. * Lr*Lr*Lr - 8./9. * (1. - 6.*rho2)*Y2NV_PHI1(rho)
            -8./27. * (19. - 46.*rho)*Y2NV_PHI2(rho) - 32./27. * (13. + 14.*rho)*Y2NV_PHI3(rho)
            -64./9. * Y2NV_PHI4(rho) - 32./9. * Lr*Li2 
            +32./27. * (5. + 9.*rho2 + 14.*rho32)*Li2 - 1792./27. * rho32*Li2sqrt
            +64./9. * Li3 + 64./9. * Li3ov + 4./3.*(2.*Lb - Lr)*dY1(E0);
}

double Bsgamma::Y2NH(double E0, double mu)
{
    double Lb = log(mu/Mb_kin);
    double zeta3 = gsl_sf_zeta_int(3);
    double Cl2 = gsl_sf_clausen(M_PI/3.);
    
    return 8./81. * (244. - 27.*sqrt(3.)*M_PI - 61.*M_PI*M_PI) 
            - 64./27.*(18. - M_PI*M_PI)*Lb - 64./9.*Lb*Lb - 64./27.*zeta3
            + 32.*sqrt(3.)*Cl2 + 8./3.*Lb*dY1(E0);
}

double Bsgamma::Y2(double E0, double mu)
{
    double CF = 4./3.;
    double CA = 3.;
    double TR = 1./2.;
    double NL = 3.;
    double NV = 1.;
    double NH = 1.;
    
    return CF*Y2CF(E0,mu) + CA*Y2CA(E0,mu) 
            + TR * ( NL*Y2NL(E0,mu) + NV*Y2NV(E0,mu) + NH*Y2NH(E0,mu) );
}

double Bsgamma::f_NLO_1(double z)
{
    return r1(2,z).real() + 2.*Phi27_1(0.,z).real();
}

double Bsgamma::zdz_f_NLO(double z, double E0)
{
    double d = delta(E0);
    double sqrt1d = sqrt(1. - d);
    double sqrt4z = sqrt(1. - 4. * z);
    double sqrt1d4z = sqrt(1. - d - 4. * z);
    double sqrtz = sqrt(z);
    double sqrt1ovz = sqrt(1./z);
    double sqrt4m1ovz = sqrt(-4. + 1./z);
    double SumSqrt = sqrt1d + sqrt1d4z;
    double ProdSqrt = sqrt1d * sqrt1d4z;
    double ProdSqrtz = sqrtz * sqrt1d4z;
    double LogSumSqrt = log(SumSqrt/(2. * sqrtz));
    double LogSqrt4z = log((1. + sqrt4z)/(2. * sqrtz));
    double LogSqrtov = log((sqrt4m1ovz + sqrt1ovz)/2.);
    
    double z2=z*z;
    double z3=z2*z;
    double z4=z3*z;
    double z5=z4*z;
    double Lz = log(z);
    
    double Pi2 = M_PI*M_PI;
    double zeta3 = gsl_sf_zeta_int(3);
    
    double zdz_f_NLO_E0;
    
    if (E0 == 0. ){
        zdz_f_NLO_E0 = 2./27. * (3. * (-7. + 2. * M_PI*M_PI) 
                + 2. * (36. - 24. * M_PI*M_PI) * z 
                + 108. * M_PI*M_PI * z2 
                - ( 36. * (-pow(1./z,3./2.)/2. - 1./(2. * sqrt4m1ovz * z2)) 
                * sqrt4m1ovz * (-1. + 2. * z))/((sqrt4m1ovz + sqrt1ovz) * pow(1./z,3./2.)) 
                - (72. * sqrt4m1ovz * LogSqrtov)/pow(1./z,3./2.) 
                - (54. * sqrt4m1ovz * (-1. + 2. * z) * LogSqrtov)/sqrt1ovz 
                + (18. * sqrt1ovz * (-1. + 2. * z) * LogSqrtov)/sqrt4m1ovz 
                - (48. * (-pow(1./z,3./2.)/2. - 1./(2. * sqrt4m1ovz * z2))
                * z * (1. - 4. * z + 6. * z2) * LogSqrtov)/(sqrt4m1ovz + sqrt1ovz) 
                - 24. * z * (-4. + 12. * z) * LogSqrtov * LogSqrtov 
                - 24. * (1. - 4. * z + 6. * z2) * LogSqrtov * LogSqrtov);
    } else zdz_f_NLO_E0 = 2. * ((2. - d) * d * (-7. + 2. * Pi2) / 9. 
            + 8./9. * d * (3. - 2. * Pi2) * z 
            + (8. * d * (-SumSqrt/( 4. * pow(z,3./2.) ) - 1./ProdSqrtz) 
            * ProdSqrt * pow(z,3./2.))/(3. * SumSqrt) 
            + 4./3. * d * ProdSqrt * LogSumSqrt 
            - (8. * (1. - d) * d * z * LogSumSqrt)/(3. * ProdSqrt) 
            - (32. * d * (-SumSqrt/( 4. * pow(z,3./2.) ) - 1./ProdSqrtz)
            * (2. - d - 4. * z) * pow(z,3./2.) * LogSumSqrt)/(9. * SumSqrt) 
            - 8./9. * d * (2. - d - 4. * z) * LogSumSqrt * LogSumSqrt 
            + 32./9. * d * z * LogSumSqrt * LogSumSqrt 
            + 4./3. * (1. - 2. * z) * z * 
            ((2. * (-( (1. + sqrt4z)/(4. * pow(z,3./2.)) ) 
            - 1./(sqrt4z * sqrtz)) * sqrt4z * sqrtz)/(1. + sqrt4z) 
            - (2. * ( -(SumSqrt/(4. * pow(z,3./2.)) ) - 1./ProdSqrtz)
            * ProdSqrt * sqrtz)/(SumSqrt) - 2. * LogSqrt4z/sqrt4z 
            + 2. * (1. - d) * LogSumSqrt/ProdSqrt) 
            + 4./3. * (1. - 2. * z) * (sqrt4z * LogSqrt4z - ProdSqrt * LogSumSqrt) 
            - 8./3. * z * (sqrt4z * LogSqrt4z - ProdSqrt * LogSumSqrt) 
            - 8./9. * z * (1 - 4. * z + 6. * z2) * 
            (( 4. * (-( (1. + sqrt4z)/(4. * pow(z,3./2.)) ) 
            - 1./(sqrt4z * sqrtz)) * sqrtz * LogSqrt4z)/(1. + sqrt4z) 
            - (4. * (-SumSqrt/(4. * pow(z,3./2.)) - 1./ProdSqrtz) * sqrtz * LogSumSqrt)/SumSqrt) 
            - 8./9. * (LogSqrt4z*LogSqrt4z - LogSumSqrt*LogSumSqrt) 
            * ( z * (-4. + 12. * z)  + (1. - 4. * z + 6. * z2) )) ;
    
    return z * (zdz_f_NLO_E0
            
            - 16./9. * (-4. + Pi2/6. - Pi2 * sqrtz 
            - Lz - z2 * (19./18. - 4. * Lz) + z * (-2. - Lz) 
            + z3 * (137./30. + 4. * Lz) 
            + z4 * (887./84. + 10. * Lz) 
            + z5 * (16597./540. + 28. * Lz) 
            - 3. * z2 * (25./12. + Pi2/9. + 19. * Lz/18. - 2. * Lz * Lz) 
            + 2. * z * (1./2. + Pi2 - 2. * Lz - Lz * Lz/2) 
            + 4. * z3 * (-1376./225. + 2. * Pi2/3 + 137. * Lz/30. + 2. * Lz * Lz) 
            + 5. * z4 * (-131317./11760. + 5. * Pi2/3. + 887. * Lz/84. + 5. * Lz * Lz) 
            + 6. * z5 * (-2807617./97200. + 14. * Pi2/3. + 16597. * Lz/540. + 14. * Lz * Lz)) 
   
            + 32./9. * (5./2. - Pi2/3. + (5./2. - 3. * Pi2/4.) * Lz 
            + Lz * Lz/4 + Lz * Lz * Lz/12. 
            + z5 * (-3303./800. - 63. * Lz/10.) 
            + z4 * (-185./144. - 35. * Lz/12.) 
            + z3 * (-1./72. - 5. * Lz/3.) 
            + z2 * (2. - 3. * Lz/2.) 
            + 6. * z5 * (67801./8000. - 21. * Pi2/20. 
            - 3303. * Lz/800. - 63. * Lz * Lz/20.) 
            + 5. * z4 * (35101./8640. - 35. * Pi2/72. 
            - 185. * Lz/144. - 35. * Lz * Lz/24.) 
            + 4. * z3 * (457./216. - 5. * Pi2/18. - Lz/72. - 5. * Lz * Lz/6.) 
            + 3. * z2 * (-7./6. - Pi2/4. + 2. * Lz - 3. * Lz * Lz/4.) 
            + z * (-Pi2/2. - Lz/2. + Lz * Lz/4.) 
            + 5./2. - 3. * Pi2/4. + Lz/2. + Lz * Lz/4. 
            + 2. * z * (7./4. + 2. * Pi2/3. - Pi2 * Lz/2. - Lz * Lz/4. 
            + Lz * Lz * Lz/12.) - 3. * zeta3));
}

double Bsgamma::mddel_f_NLO(double z, double E0)
{
    double d = delta(E0);
    double d2 = d*d;
    double sqrt1d = sqrt(1. - d);
    double sqrt1d4z = sqrt(1. - d - 4. * z);
    double sqrtz = sqrt(z);
    double LogSqrt = log((sqrt1d + sqrt1d4z)/(2. * sqrtz));
    double SumSqrt = sqrt1d + sqrt1d4z;
    double ProdSqrt = sqrt1d * sqrt1d4z;
    
    return 2. * (1. - d) * ( -2./27. * d * (-3. + 2. * d) 
            - 2./27. * (3. - 3. * d + d2) 
            + 1./9. * (2. - d) * (-7. + 2. * M_PI * M_PI) * z 
            - 1./9. * d * (-7. + 2. * M_PI * M_PI) * z 
            + 4. * d * (-1./(2. * sqrt1d) - 1./(2. * sqrt1d4z))
            * ProdSqrt * z  / (3. * SumSqrt) 
            + 4./9. * (3. - 2. * M_PI * M_PI) * z * z 
            + 4./3. * ProdSqrt * z * LogSqrt 
            - 16. * d * (-1./(2. * sqrt1d) - 1./(2. * sqrt1d4z)) 
            * (2. - d - 4. * z) * z * LogSqrt / (9. * SumSqrt) 
            + 2. * d * z * (-2. + 2. * d + 4. * z) * LogSqrt / (3. * ProdSqrt) 
            + 16. * (-1./(2. * sqrt1d) - 1./(2. * sqrt1d4z)) * z
            * (1 - 4. * z + 6. * z * z) * LogSqrt / (9. * SumSqrt) 
            + 8./9. * d * z * LogSqrt * LogSqrt 
            - 8./9. * (2. - d - 4. * z) * z * LogSqrt * LogSqrt 
            + 4./3. * (1. - 2. * z) * z 
            * ( (1./(2. * sqrt1d) + 1./(2. * sqrt1d4z)) * ProdSqrt/SumSqrt 
            - (-2. + 2. * d + 4. * z) * LogSqrt/(2. * ProdSqrt)));
}

double Bsgamma::h27_2(double z, double E0)
{
    double d = delta(E0);
    double d2 = d*d;
    
    if (E0 == 0.){
        return ( 41./27. - 2./9. * M_PI*M_PI
            - 2.24 * sqrt(z) - 7.04 * z + 23.72 * pow(z,3./2.)
            + ( -9.86 * z + 31.28 * z * z ) * log(z));
    } else return - 0.1755402735503456 - 1.4553730660088837 * d 
            + 1.1192806367180177 * d2 
            + ( 0.7259818237183779 - 7.230418135384073 * d 
            + 5.977206932166958 * d2 ) * sqrt(z) 
            + ( 13.786205094458156 + 113.71026116073105 * d 
            - 100.3588074342665 * d2 ) * z 
            + ( -145.05588751363894 - 307.05884309429547 * d 
            + 388.54181686721904 * d2 ) * pow(z,3./2.) 
            + ( 475.2039505292043 + 312.9832308573048 * d 
            - 775.8088176670707 * d2 ) * z * z
            + ( -509.7299390734172 - 126.08888075477071 * d 
            + 646.2084041395774 * d2 ) * pow(z,5./2.);
}

double Bsgamma::f_q(double z, double E0)
{
    return Rer22(z) - 4./3. * h27_2(z,E0);
}

double Bsgamma::f_b(double z)
{
    return -1.836 + 2.608 * z + 0.8271 * z * z - 2.441 * z * log(z);
}

double Bsgamma::f_c(double z)
{
    return 9.099 + 13.20 * z - 19.68 * z * z + 25.71 * z * log(z);
}

double Bsgamma::F_1(double z)
{
    return - 23.74697061848885 + 35./12. * f_q(z,0.) 
            + (2129./936. - 9./52. * M_PI*M_PI) * f_NLO_1(z)
            - 0.8444138663102 * zdz_f_NLO(z,0.);
}

double Bsgamma::F_2(double z)
{
    return - 3.006537367876035 - 592./81. * f_q(z,0.) 
            - 10.344289655256379 * f_NLO_1(z)
            - 9.550817514525745 * zdz_f_NLO(z,0.);
}

double Bsgamma::delddel_Phi22_1(double E0)
{
    double d = delta(E0);
    
    return 4. * (1. - d)/d * 16./27. * Int_cc1_part1(E0);
}

double Bsgamma::zdz_Phi22_1(double E0)
{
    return 64./27. * ( Int_cc1(E0) - Int_cc(E0));
}

double Bsgamma::delddel_Phi28_1(double z, double E0)
{
    double d = delta(E0);
    double Sq = sqrt( (1. - d) * (1. - d - 4.*z) );
    double Log = log( ( sqrt(1. - d) + sqrt(1. - d - 4.*z) ) / 2. / sqrt(z) );
    double Log2 = Log*Log;
    
    return 4. / (27. * Sq * (1 - d - 4. * z)) * Sq * Sq *
            (-8. * (-1. + d) * Log * z * (-1. + d + 4. * z) + 
            Sq * (1. + d*d + (4. + 8. * Log2 - 2. * M_PI*M_PI) * z + 
            4. * (-4. * Log2 + M_PI*M_PI) * z * z - 
            2. * d * (1. + (2. + 4. * Log2 - M_PI*M_PI) * z)));
}

double Bsgamma::zdz_Phi28_1(double z, double E0)
{
    double d = delta(E0);
    double Sq = sqrt( (1. - d) * (1. - d - 4.*z) );
    double Log1 = log( ( sqrt(1. - d) + sqrt(1. - d - 4.*z) ) / 2. / sqrt(z) );
    double Log2 = log( ( 1. + sqrt(1. - 4.*z) ) / 2. / sqrt(z) );
    
    return 2./27. * z * (d*d * (-7. - 8. * (-1. + Log1) * Log1 + 2. * M_PI*M_PI) 
            - 20. * Log2 * sqrt(1. - 4. * z) + 72. * Log2 * sqrt(1. - 4. * z) * z 
            + (48. * Log1 * Sq * z * z) / (-1. + d + 4. * z) 
            + (8. * z * ( -(3. + 4. * Log1) * Sq*Sq + 2. * (3. + 8. * Log1) * Sq * z 
            - 24. * Log1 * z * z )) / (-1. + d - Sq + 4. * z) 
            - 2. * d * (-7. - 3. * Sq + Log1 * (8. + 6. * Sq) + M_PI*M_PI * (2. - 8. * z) 
            + 12. * z + 8. * Log1*Log1 * (-1. + 4. * z)) 
            + 2. * (-3. * (-1. + Sq + 2. * (1. + Sq) * z) 
            + Log1 * (4. + Sq * (6. - 52. * z) + 24. * z * z) 
            + 4. * ( Log2 * Log2 - Log1 * Log1) * (1. + 2. * z * (-4. + 9. * z))));
}

double Bsgamma::delddel_Phi88_1(double E0)
{
    double d = delta(E0);
    double Ld = log(d);
    
    return 4./27. * (1. - d) * (5. - 8./(1. - d) + 5. * d - 2. * d * d - 
            2. * (2. - 4./(1. - d) + 2. * d) * log(Mb_kin/Ms) + (4. * Ld)/(1. - d)
            - d * Ld - (2. + d) * Ld);
}

double Bsgamma::f(double r)
{
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r3*r;
    double Lr = log(r);
    double zeta3 = gsl_sf_zeta_int(3);
    Polylogarithms Poly;
    
    if (r==1.) 
        return 7126./81. - 356./27. * M_PI*M_PI - 16./3. * zeta3;
    else
        return - 16./3. * Poly.Li3(r2) 
                + 8. * r * (35./9. * r2 + 9.) * ( gsl_sf_dilog(r) - atanh(r) * Lr - 1./4. * M_PI*M_PI) 
                + 2. * (8./3. * Lr - 6. * r4 - 35./9. * r3 - 9. * r - 62./9. ) * gsl_sf_dilog(r2) 
                - 8. * (3. * r4 + 31./9.) * log(1. - r2) * Lr + 32./9. * Lr*Lr*Lr 
                + 8. * (3. * r4 + 25./9.) * Lr*Lr + 64./9. * r2 * Lr 
                + 4. * M_PI*M_PI * r4 + 172./9. * r2 + 5578./81.;
}

double Bsgamma::Delta(double r) 
{
    double r2 = r*r;
    double r3 = r2*r;
    double Lmr = log(1. - r);
    double Lpr = log(1. + r);
    double Lr = log(r);
    double Lr2 = Lr*Lr;
    
    if (r==1.) 
        return -3./8. + 1./8. * M_PI*M_PI;
    else
        return 1./4. * (1. - r) * (1. - r3) * ( gsl_sf_dilog(r) + Lr * Lmr - 1./2. * Lr2 - 1./3. * M_PI*M_PI ) 
                - 1./4. * (1. + r) * (1. + r3) * ( gsl_sf_dilog(-1./r) - Lr * Lpr + Lr2 ) 
                + 1./4. * Lr2 - 1./4. * r2 * Lr - 3./8. * r2 + 1./24. * M_PI*M_PI;
}

double Bsgamma::f_u(double r) 
{
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r3*r;
    double r5 = r4*r;
    double r6 = r5*r;
    double r7 = r6*r;
    double Lr = log(r);
    double Lr2 = Lr*Lr;
    double Pi2 = M_PI*M_PI;
    double zeta3 = gsl_sf_zeta(3.);
    
    if (r==1.) 
        return 6335./288. - 1./2. * M_PI*M_PI - 16. * zeta3;
    else
        return -5./6. * Pi2 * r + ( 14. + 16./9. * Pi2) * r2 
                + (64./9. * Lr + 128./9. * log(2.) - 95./54.) * Pi2 * r3 
                + (-16./3. * Lr2 + 365./9. * Lr + 4. * Pi2 * Lr 
                - 4375./54. - 25./9. * Pi2 + 32. * zeta3) * r4 
                - 224./45. * Pi2 * r5 + (-128./27. * Lr2 + 16./15. * Lr 
                + 15608./2025. + 128./81. * Pi2) * r6 - 16./7. * Pi2 * r7;
}

double Bsgamma::omega77(double z)
{
    double z2 = z*z;
    double z3 = z2*z;
    double z4 = z3*z;
    double z5 = z4*z;
    double z6 = z5*z;
    double z7 = z6*z;
    double z8 = z7*z;
    double omz = 1. - z;
    double omz2 = omz * omz;
    double omz3 = omz2 * omz;
    double Lz = log(z);
    double Lomz = log (1. - z);
    double Lomz2 = Lomz*Lomz;
    double Lomz3 = Lomz2*Lomz;
    double Ltmz = log (2. - z);
    double Ltmz2 = Ltmz*Ltmz;
    double Li2omz = gsl_sf_dilog(1. - z);
    double Li2zmo = gsl_sf_dilog(z - 1.);
    Polylogarithms Poly;
    double Pi2 = M_PI*M_PI;
    
    return 4./9. * (z3 - 4. * z2 + 4. * z + 1.)/omz * ( 2. * Poly.Li3( 1./(2.-z) ) 
            - Poly.Li3( z/(2.-z) ) + Poly.Li3( z/(z-2.) ) 
            + Ltmz * ( Lomz2 - 1./3. * Ltmz2 + 1./6. * Pi2 ) ) 
            + 4./9. * (z3 + 36. * z - 43.)/omz * Poly.Li3(z) 
            + 8./9. * (z3 - 2. * z2 + 19. * z - 22.)/omz * Poly.Li3(1.-z) 
            - 16./9. * omz2 * Poly.Li3(z-1.) 
            - 4./9. * (z3 + 35. * z - 44.)/omz * Li2omz * Lomz 
            - 4./9. * (z3 - 2. * z2 + 2. * z - 3.)/omz * Li2zmo * Lomz 
            - 4./27. * (23. * z6 - 106. * z5 + 145. * z4 + 3. * z3 
            - 180. * z2 + 147. * z - 36.)/(z * omz3) * (Li2omz + Lomz * Lz) 
            + 2./27. * (z8 - 6. * z7 + 9. * z6 + 27. * z5 - 140. * z4 + 219. * z3 
            - 124. * z2 + 28. * z - 6.)/(z * omz3) * (Li2zmo + Lomz * Ltmz) 
            - 8./9. * (z2 + 8. * z - 11.)/omz * Lomz2 * Lz 
            - 2./9. * (z4 - 3. * z3 - 5. * z2 + 15. * z + 8.)/(z * omz) * Lomz3 
            - (z6 - 4. * z5 - 46. * z4 + 101. * z3 - 461. * z2 + 1057. * z - 72.)/(27. * z * omz) * Lomz2 
            + 2./27. * (z3 - 2. * z2 + 4. * z - 5.)/omz * Pi2 * Lomz 
            + (2. * z5 - 29. * z4 - 113. * z3 + 153. * z2 - 827. * z - 162.)/(27. * z * omz) * Lomz 
            - (3. * z3 - 8. * z2 + 144. * z - 157.)/(9. * omz) * gsl_sf_zeta(3.) 
            + (z6 - 4. * z5 + 48. * z4 - 106. * z3 - 58. * z2 + 158. * z - 75.)/(81. * z * omz) * Pi2 
            + (2. * z4 - 92. * z3 + 88. * z2 - 713. * z - 18.)/(27. * omz);
}

double Bsgamma::Int_Phi77_2rem(double E0)
{
    if (IntPhi772rCached == 0) {
        double t1 = (1. - delta(E0));

        INT = convertToGslFunction(boost::bind(&Bsgamma::omega77, &(*this), _1));
        if (gsl_integration_cquad(&INT, 0., t1, 1.e-2, 1.e-1, w_INT, &avaINT, &errINT, NULL) != 0) return std::numeric_limits<double>::quiet_NaN();
        
        CacheIntPhi772r = avaINT;
        IntPhi772rCached = 1;
    }
    
    return CacheIntPhi772r;
}

double Bsgamma::Phi77_2rem(double E0)
{
    double xm = 8./9. * M_PI * alsUps;
    double d = delta(E0);
    double d2 = d*d;
    double d3 = d2*d;
    double d4 = d3*d;
    
    return Int_Phi77_2rem(E0) - xm/3./ d * 
            ((4. - 6. * d2 + 2. * d3) * log(d) + 7. - 13. * d + 3. * d2 + 5. * d3 - 2. * d4);
}

double Bsgamma::K77_2_z1(double E0, double mu)
{
    double K77_1 = Kij_1(7,7,E0,mu).real();
    double Pi2 = M_PI*M_PI;
    double xm = 8./9. * M_PI * alsUps;
    double Lb = 2.*log(mu_b/Mb_kin);
    
    return ( K77_1 - 4. * Phi77_1(E0) ) * K77_1 - 1178948./729. + 18593./729. * Pi2 
            - 628./405. * Pi2*Pi2 + 428./27. * Pi2 * log(2.) + 61294./81. * gsl_sf_zeta(3.) 
            - 880./9. * Lb * Lb + ( 440./27. * Pi2 - 14698./27. ) * Lb 
            + 64./3. * xm + 4. * (Phi77_2beta0(E0,mu) + Phi77_2rem(E0));
}

double Bsgamma::Kij_2(int i, int j, double E0, double mu_b, double mu_c)
{
    if (i<1 || i > 2)
        if (i < 7)
            throw std::runtime_error("Bsgamma::Kij_2(): index i must be included in (1,2,7,8)");
    if (j<1 || j > 2)
        if (j < 7)
            throw std::runtime_error("Bsgamma::Kij_2(): index j must be included in (1,2,7,8)");
    
    int temp;
    
    if (i > j) {temp=i; i=j; j=temp;}
    
    double K_ij[8][8] = {{0.}};
    double z = zeta();
    double d = delta(E0);
    double r = sqrt(z);
    double Lb = 2.*log(mu_b/Mb_kin);
    double Lb2 = Lb*Lb;
    double Lc = 2.*log(mu_c/Mc);
    double Lcb = log(Mc/Mb_kin);
    double xm = 8./9. * M_PI * alsUps;
    
    double A1 = 22.604961348474838;
    double A2 = 75.60281607240395;
    
    K_ij[1][1] = (r1(2,zeta()).real() - 208./81.*Lb) * (r1(2,zeta()).real() - 208./81.*Lb)
            + r1(2,zeta()).imag() * r1(2,zeta()).imag()
            + 4.*Phi22_2beta0(E0,mu_b) * SM.Beta0(3)/SM.Beta0(5) + 16./3. * Phi22_1(E0) * (Lcb - Lb)
            + xm * (delddel_Phi22_1(E0) - 2. * zdz_Phi22_1(E0));
    K_ij[1][6] = A2 + F_2(z) - 27./2. * f_q(z,E0) + f_b(z) + f_c(z)
            + 4./3. * Phi27_1(E0,z).real() * log(z) + (8. * Lc - 2. * xm) * zdz_f_NLO(z,E0)
            + xm * mddel_f_NLO(z,E0) + 416./81. * xm
            + (10./3. * Kij_1(2,7,E0,mu_b).real() - 2./3. * Kij_1(4,7,E0,mu_b).real()
            - 208./81. * Kij_1(7,7,E0,mu_b).real() - 35./27. * Kij_1(7,8,E0,mu_b).real()
            - 254./81.) * Lb - 5948./729. * Lb2;
    K_ij[1][7] = (r1(2,zeta()).real() - 208./81.*Lb) * (r1(8,zeta()).real() + 16./9.*Lb)
            + r1(2,zeta()).imag() * r1(8,zeta()).imag()
            + 2.*Phi28_2beta0(E0,mu_b) * SM.Beta0(3)/SM.Beta0(5) + 8./3. * Phi28_1(E0,z).real() * (Lcb - Lb)
            + xm * (delddel_Phi28_1(z,E0) - 2. * zdz_Phi28_1(z,E0));
    
    K_ij[0][0] = 1./36. * K_ij[1][1];
    K_ij[0][1] = -1./6. * K_ij[1][1];
    K_ij[0][6] = - 1./6. * K_ij[1][6] + A1 + F_1(z)
            + (- 3./2. * Kij_1(2,7,E0,mu_b).real() - 3./4. * Kij_1(7,8,E0,mu_b).real()
            + 94./81.) * Lb - 34./27. * Lb2;
    K_ij[0][7] = -1./6. * K_ij[1][7];
    
    K_ij[6][6] = K77_2_z1(E0,mu_b) + ( 1972./81. - 16./27. * M_PI*M_PI + 8./3. * Phi77_1(E0)) * log(zeta())
            + 2./3. * (f(r) - f(1.)) - 128./3. * (Delta(r) - Delta(1.)) 
            - 16. * (f_u(r) - f_u(1.));
    K_ij[6][7] = 2./3. * Y2(E0,mu_b) + (16./9.*M_PI*M_PI - 164./9. - 32./6. * Lb) * Y1(E0,mu_b) 
            - 32./81. * alsUps * M_PI * (3. + 7.*d - 3.*d*d + d*d*d - 4.*d*log(d) );
    
    K_ij[7][7] = (r1(8,zeta()).real() + 16./9.*Lb) * (r1(8,zeta()).real() + 16./9.*Lb)
            + r1(8,zeta()).imag() * r1(8,zeta()).imag()
            + 4.*Phi88_2beta0(E0,mu_b) * SM.Beta0(3)/SM.Beta0(5) + 16./3. * Phi88_1(E0) * (Lcb - Lb)
            + xm * (delddel_Phi88_1(E0));
    
    return K_ij[i-1][j-1];
}

void Bsgamma::computeCoeff(double mu)
{
    
    /*allcoeff = SM.getMyFlavour()->ComputeCoeffsgamma(160.);
    
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
    
    std::cout << "C_0(MuW): (" << C1_0.real() << "," << C2_0.real() << "," 
            << C3_0.real() << "," << C4_0.real() << "," << C5_0.real() << "," 
            << C6_0.real() << "," << C7_0.real() << "," << C8_0.real() << ")" << std::endl;
    std::cout << "C_1(MuW): (" << C1_1.real() << "," << C2_1.real() << "," 
            << C3_1.real() << "," << C4_1.real() << "," << C5_1.real() << "," 
            << C6_1.real() << "," << C7_1.real() << "," << C8_1.real() << ")" << std::endl << std::endl;*/
    
    allcoeff = SM.getFlavour().ComputeCoeffsgamma(mu);
    allcoeffprime = SM.getFlavour().ComputeCoeffprimesgamma(mu);
    
    C1_0 = (*(allcoeff[LO]))(0);
    C2_0 = (*(allcoeff[LO]))(1);
    C3_0 = (*(allcoeff[LO]))(2);
    C4_0 = (*(allcoeff[LO]))(3);
    C5_0 = (*(allcoeff[LO]))(4);
    C6_0 = (*(allcoeff[LO]))(5);
    C7_0 = (*(allcoeff[LO]))(6)  + C_7_NP;
    C8_0 = (*(allcoeff[LO]))(7);
    
    C1_1 = (*(allcoeff[NLO]))(0)/Alstilde;
    C2_1 = (*(allcoeff[NLO]))(1)/Alstilde;
    C3_1 = (*(allcoeff[NLO]))(2)/Alstilde;
    C4_1 = (*(allcoeff[NLO]))(3)/Alstilde;
    C5_1 = (*(allcoeff[NLO]))(4)/Alstilde;
    C6_1 = (*(allcoeff[NLO]))(5)/Alstilde;
    C7_1 = (*(allcoeff[NLO]))(6)/Alstilde;
    C8_1 = (*(allcoeff[NLO]))(7)/Alstilde;

    C7_2 = (*(allcoeff[NNLO]))(6)/Alstilde/Alstilde;

    C7p_0 = (*(allcoeffprime[LO]))(6) + Ms/Mb*((*(allcoeff[LO]))(6)) + C_7p_NP;
    C7p_1 = ((*(allcoeffprime[NLO]))(6) + Ms/Mb*((*(allcoeff[NLO]))(6)))/Alstilde; /*Implement the other WCs*/

    /*std::cout << "C_0(mu): (" << C1_0.real() << "," << C2_0.real() << "," 
            << C3_0.real() << "," << C4_0.real() << "," << C5_0.real() << "," 
            << C6_0.real() << "," << C7_0.real() << "," << C8_0.real() << ")" << std::endl;
    std::cout << "C_1(mu): (" << C1_1.real() << "," << C2_1.real() << "," 
            << C3_1.real() << "," << C4_1.real() << "," << C5_1.real() << "," 
            << C6_1.real() << "," << C7_1.real() << "," << C8_1.real() << ")" << std::endl << std::endl;
    std::cout << "C_2^7(mu): " << C7_2.real() << std::endl << std::endl;*/
    
    C7_1ew = 4.868;
    
}

double Bsgamma::P0(double E0)
{
    return C7_0.abs2() + C7p_0.abs2() + P0_4body(E0,Mb_kin*Mb_kin/Ms/Ms);
}

double Bsgamma::P11()
{
    return 2.*( C7_0.real()*C7_1.real() + C7_0.imag()*C7_1.imag()
            + C7p_0.real()*C7p_1.real() + C7p_0.imag()*C7p_1.imag() ); /*CHECK SIGN*/
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
            p21 += ( C0[i].real()*C0[j].real() + C0[i].imag()*C0[j].imag() ) * Kij_1(i+1,j+1,E0,mu).real();
        }
    }
    
    for(i=6;i<7;i++) /*CHECK ALGORITHM*/
    {
        for(j=6;j<7;j++)
        {
            p21 += (C0p[i].real()*C0p[j].real() + C0p[i].imag()*C0p[j].imag()) * Kij_1(i+1,j+1,E0,mu).real();
        }
    }
    
    return p21;
}

double Bsgamma::P21_CPodd(double E0, double mu)
{
    int i,j;
    gslpp::complex C0[8]={C1_0,C2_0,C3_0,C4_0,C5_0,C6_0,C7_0,C8_0};
    gslpp::complex C0p[8]={C7p_0}; /*IMPLEMENT OTHER WC*/
    double p21=0.;
    
    for(i=0;i<8;i++)
    {
        for(j=0;j<8;j++)
        {
            p21 += - ( C0[i].real()*C0[j].imag() - C0[i].imag()*C0[j].real() ) * Kij_1(i+1,j+1,E0,mu).imag();
        }
    }
    
    for(i=6;i<7;i++) /*CHECK ALGORITHM*/
    {
        for(j=6;j<7;j++)
        {
            p21 += - ( C0p[i].real()*C0p[j].imag() - C0p[i].imag()*C0p[j].real() ) * Kij_1(i+1,j+1,E0,mu).imag();
        }
    }
    
    return p21;
}

double Bsgamma::P12()
{
    
   return C7_1.abs2() + C7p_1.abs2() + 2.*(C7_0*C7_2).real(); /*CHECK SIGN*/
}

double Bsgamma::P22(double E0, double mu_b, double mu_c)
{
    
    int i,j, temp_i,temp_j;
    gslpp::complex C0[4]={C1_0,C2_0,C7_0,C8_0};
    double p22=0.;
    
    for(i=0;i<4;i++)
    {
        for(j=0;j<4;j++)
        {
            if (i > 1) {
                temp_i=i+4;
            } else temp_i=i;
            if (j > 1) {
                temp_j=j+4;
            } else temp_j=j;
            p22 += (C0[i]*C0[j]).real() * Kij_2(temp_i+1,temp_j+1,E0,mu_b,mu_c);
        }
    }
    
    return p22;
}

double Bsgamma::P32(double E0, double mu)
{
    
    int i,j;
    gslpp::complex C0[8]={C1_0,C2_0,C3_0,C4_0,C5_0,C6_0,C7_0,C8_0};
    gslpp::complex C1[8]={C1_1,C2_1,C3_1,C4_1,C5_1,C6_1,C7_1,C8_1};
    double p32=0.;
    
    for(i=0;i<8;i++)
    {
        for(j=0;j<8;j++)
        {
            p32 += 2.*(C0[i]*C1[j]).real() * Kij_1(i+1,j+1,E0,mu).real();
        }
    }

    return p32;
}

double Bsgamma::EW_NLO(double mu)
{
    
    if(EWflag) {
        double ew_nlo = 0.;
        double ga_eff_ew_7[7] = {-832./729., -208./243., -20./243., -176./729., -22712./243., -6272./729., 16./9.};
        double Lb = log(mu/Mb_kin);
        double Lz = 2. * log(Mz/mu);
        gslpp::complex C[7] = {C1_0, C2_0, C3_0, C4_0, C5_0, C6_0, C7_0};
        gslpp::complex r[7] = {0.};
        
        r[0] = r1_ew(1,zeta()) - ga_eff_ew_7[0] * Lb;
        r[1] = r1_ew(2,zeta()) - ga_eff_ew_7[1] * Lb;
        r[2] = r1_ew(3,zeta()) - ga_eff_ew_7[2] * Lb;
        r[3] = r1_ew(4,zeta()) - ga_eff_ew_7[3] * Lb;
        r[4] = r1_ew(5,zeta()) - ga_eff_ew_7[4] * Lb;
        r[5] = r1_ew(6,zeta()) - ga_eff_ew_7[5] * Lb;
        r[6] = r1_ew(7,zeta()) - ga_eff_ew_7[6] * Lb;
    
        for(int i=0;i<7;i++){
            ew_nlo += 2. * (C7_0.real()*C[i].real() + C7_0.imag()*C[i].imag()) * r[i].real()
                - 2. * (C7_0.real()*C[i].imag() - C7_0.imag()*C[i].real()) * r[i].imag();
        }
        
        ew_nlo += 2. * (C7_0.real() * (-2. * C7_0 * Lz + C7_1ew).real() 
                + C7_0.imag() * (-2. * C7_0 * Lz + C7_1ew).imag());
        
        return ew_nlo;
    }
    
    else return 0.;
}

double Bsgamma::Vub_NLO_2body()
{
    double z = zeta();
    
    return 4. * Alstilde * 
            ((C7_0.real()*( C2_0-C1_0/6. ).real() + C7_0.imag()*( C2_0-C1_0/6. ).imag()) * CKMu.real() +
            (C7_0.real()*( C2_0-C1_0/6. ).imag() - C7_0.imag()*( C2_0-C1_0/6. ).real()) * CKMu.imag())
            * ( a(z)+b(z) ).real();
}

double Bsgamma::Vub_NLO_2body_CPodd()
{
    double z = zeta();
    
    return - 4. * Alstilde * 
            ((C7_0.real()*( C2_0-C1_0/6. ).real() + C7_0.imag()*( C2_0-C1_0/6. ).imag()) * CKMu.imag() +
            (C7_0.real()*( C2_0-C1_0/6. ).imag() - C7_0.imag()*( C2_0-C1_0/6. ).real()) * CKMu.real() ) 
            * ( a(z)+b(z) ).imag();
}

double Bsgamma::Vub_NLO_3body_A(double E0)
{
    double d = delta(E0);
    
    return 64./27. * Alstilde * ( C2_0 - C1_0/6. ).abs2() *
            ( CKMu.real() * ( 2. * Int_cc1(E0) - Int_c1(E0).real() )
            + CKMusq *  ( Int_cc1(E0) - Int_c1(E0).real() + 1./8. * d * ( 1. - 1./3. * d*d ) ));
}

double Bsgamma::Vub_NLO_3body_A_CPodd(double E0)
{
    return - 64./27. * Alstilde * ( C2_0 - C1_0/6. ).abs2() * CKMu.imag() * Int_c1(E0).imag();
}

double Bsgamma::Vub_NLO_3body_B(double E0)
{
    double d = delta(E0);
    gslpp::complex wc1 = C7_0 - C8_0/3.;
    gslpp::complex wc2 = C2_0 - C1_0/6.;
    gslpp::complex me = Phi27_1(E0,zeta()) + 2./9. * d * ( 1. - d + 1./3. * d*d );
    
    return 4. * Alstilde * 
            (( wc1.real() * wc2.real() + wc1.imag() * wc2.imag()) * CKMu.real() +
            ( wc1.real() * wc2.imag() - wc1.imag() * wc2.real()) * CKMu.imag() )
            * me.real();
}

double Bsgamma::Vub_NLO_3body_B_CPodd(double E0)
{
    double d = delta(E0);
    gslpp::complex wc1 = C7_0 - C8_0/3.;
    gslpp::complex wc2 = C2_0 - C1_0/6.;
    gslpp::complex me = Phi27_1(E0,zeta()) + 2./9. * d * ( 1. - d + 1./3. * d*d );
    
    return - 4. * Alstilde * 
            (( wc1.real() * wc2.real() + wc1.imag() * wc2.imag()) * CKMu.imag() +
            ( wc1.real() * wc2.imag() - wc1.imag() * wc2.real()) * CKMu.real() )
            * me.imag();
}

double Bsgamma::Vub_NLO_4body(double E0)
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

        return 4. * Alstilde * ( 
                ( C2_0 - C1_0/6. ).abs2() * CKMu.real() * 0.005025213076791178
                + (( C2_0 - C1_0/6. ).real() * (C7_0.real() * uphib427 + C8_0.real() * uphib428)
                + ( C2_0 - C1_0/6. ).imag() * (C7_0.imag() * uphib427 + C8_0.imag() * uphib428) ) * CKMu.real() 
                + (( C2_0 - C1_0/6. ).real() * (C7_0.imag() * uphib427 + C8_0.imag() * uphib428)
                - ( C2_0 - C1_0/6. ).imag() * (C7_0.real() * uphib427 + C8_0.real() * uphib428) ) * CKMu.imag() );
    }
    
    else return 0.;
}

double Bsgamma::Vub_NLO_4body_CPodd(double E0)
{
    if (FOUR_BODY) {
        return 4. * Alstilde * ( C2_0 - C1_0/6. ).abs2() * CKMu.imag() * 0.013978889449487913;
    }
    
    else return 0.;
}

double Bsgamma::Vub_NLO(double E0)
{
    return Vub_NLO_2body() + Vub_NLO_3body_A(E0) + Vub_NLO_3body_B(E0) + Vub_NLO_4body(E0);
}

double Bsgamma::Vub_NLO_CPodd(double E0)
{
    return Vub_NLO_2body_CPodd() + Vub_NLO_3body_A_CPodd(E0) + Vub_NLO_3body_B_CPodd(E0) + Vub_NLO_4body_CPodd(E0);
}

double Bsgamma::Vub_NNLO(double E0)
{
    double r12 = (( C2_1 - C1_1/6. )/( C2_0 - C1_0/6. )).real();
    double r78 = (( C7_1 - C8_1/3. )/( C7_0 - C8_0/3. )).real();
    double r7 = (C7_1/C7_0).real();
    
    return Alstilde * ( (r12 + r7) * Vub_NLO_2body() 
            + 2. * r12 * Vub_NLO_3body_A(E0) + (r12 + r78) * Vub_NLO_3body_B(E0));
}

double Bsgamma::P(double E0, double mu_b, double mu_c, orders order)
{

    switch(order) {
        case NNLO:
            /*std::cout << "p0 w/ tree, VubLO: " << P0(E0) << std::endl;
            std::cout << "p11: " << P11() << std::endl;
            std::cout << "p21: " << P21(E0,mu_b) << std::endl;
            std::cout << "p12: " << P12() << std::endl;
            std::cout << "p22: " << P22(E0,mu_b,mu_c) << std::endl;
            std::cout << "p32: " << P32(E0,mu_b) << std::endl;
            std::cout << "Vub_NLO: " << Vub_NLO(E0) << std::endl;
            std::cout << "Vub_NNLO: " << Vub_NNLO(E0) << std::endl;
            std::cout << "EW_NLO: " << EW_NLO(mu_b) << std::endl;*/
            return P0(E0) 
                    + Alstilde * (P11() + P21(E0,mu_b)) + Vub_NLO(E0) + AleatMztilde * EW_NLO(mu_b)
                    + Alstilde * Alstilde * (P12() + P22(E0,mu_b,mu_c) + P32(E0,mu_b)) + Vub_NNLO(E0);
            break;
        case NLO:
            return P0(E0) 
                    + Alstilde * (P11() + P21(E0,mu_b)) + Vub_NLO(E0) + AleatMztilde * EW_NLO(mu_b);
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
    corrLambda2_mix = 1./8. * (9. * lambda2 - lambda1) * Kij_1(7,7,E0,mu).real();
    
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
    mu_kin=SM.getOptionalParameter("Gambino_mukin");
    BRsl=SM.getOptionalParameter("Gambino_BRsem")/100.;
    Mb_kin=SM.getOptionalParameter("Gambino_Mbkin");
    Mc=SM.getOptionalParameter("Gambino_Mcatmuc");
    mu_pi2=SM.getOptionalParameter("Gambino_mupi2");
    rho_D3=SM.getOptionalParameter("Gambino_rhoD3");
    mu_G2=SM.getOptionalParameter("Gambino_muG2");
    rho_LS3=SM.getOptionalParameter("Gambino_rhoLS3");
    mu_b=SM.getOptionalParameter("mu_b_bsgamma");
    mu_c=SM.getOptionalParameter("mu_c_bsgamma");
    C=C_sem();
    
    ale=SM.getAle();
    alsUps=8./M_PI * mu_kin/Mb_kin * ( 1. + 3./8. * mu_kin/Mb_kin );
    Alstilde = SM.Alstilde5(mu_b);
    AleatMztilde=SM.ale_OS(SM.getMz())/4./M_PI;
    Ms=SM.getQuarks(QCD::STRANGE).getMass();
    Mb=SM.getQuarks(QCD::BOTTOM).getMass();
    Mz=SM.getMz();
    V_ub=SM.getCKM().V_ub();
    V_cb=SM.getCKM().V_cb();
    V_tb=SM.getCKM().V_tb();
    
    if(WET_NP_btos){
        C_7_NP = SM.getOptionalParameter("C7_NP");
        C_7p_NP = SM.getOptionalParameter("C7p_NP");
    }
    else if(SMEFT_NP_btos){
            gslpp::complex SMEFT_factor = (M_PI/SM.getAle())*(SM.v()*1.e-3)*(SM.v()*1.e-3)/SM.computelamt_s();
            C_7_NP = SM.getOptionalParameter("CdB")-SM.getOptionalParameter("CdW");
            C_7_NP *= SMEFT_factor*SM.getAle()*8.*M_PI*SM.v()/Mb;
            C_7p_NP = SM.getOptionalParameter("CpdB")-SM.getOptionalParameter("CpdW");
            C_7p_NP *= SMEFT_factor*SM.getAle()*8.*M_PI*SM.v()/Mb;
    }
    else{
        C_7_NP = 0.;
        C_7p_NP = 0.;        
    }
    
    if (SUM) {
        CKMratio=(V_tb/V_cb).abs2()*(1. - V_tb.abs2());
        CKMu=-V_ub.abs2()/(1. - V_tb.abs2());
        CKMusq = (V_ub/V_tb).abs2() * (1. - V_ub.abs2())/(1. - V_tb.abs2());
        }
    else
        switch (quark) {
            case StandardModel::STRANGE:
                CKMratio=(SM.computelamt_s()/V_cb).abs2();
                CKMu=SM.computelamu_s().conjugate() / SM.computelamt_s().conjugate(); // -0.00802793 + 0.0180942*gslpp::complex::i(); //
                CKMusq = CKMu.abs2();
                break;
            case StandardModel::DOWN:
                CKMratio=(SM.computelamt_d()/V_cb).abs2();
                CKMu=SM.computelamu_d().conjugate() / SM.computelamt_d().conjugate(); //  0.00745398 - 0.40416*gslpp::complex::i(); //
                CKMusq = CKMu.abs2();
                break;
            default:
                std::stringstream out;
                out << quark;
                throw std::runtime_error("bqgamma: quark " + out.str() + " not implemented");
        }
    
    BLNPcorr=SM.getOptionalParameter("BLNPcorr");
    
    checkCache();
    
    if (Intb_updated == 0) {
        Intb1Cached = 0;
        Intb2Cached = 0;
        Intb3Cached = 0;
        Intb4Cached = 0;
        Intbb1Cached = 0;
        Intbb2Cached = 0;
        Intbb4Cached = 0;
        IntPhi772rCached = 0;
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
    
    overall = BRsl * CKMratio * 6. * ale / (M_PI * C);
    
}

double Bsgamma::computeThValue()
{
    double E0 = getBinMin();
    
    updateParameters();
    
    if (obs == 1)
        return overall *  ( P(E0, mu_b, mu_c, NNLO) + N(E0, mu_b) );
    if (obs == 2)
        return (Alstilde * P21_CPodd(E0, mu_b) + Vub_NLO_CPodd(E0) ) / (P(E0, mu_b, mu_c, NNLO) + N(E0, mu_b) );
    
    throw std::runtime_error("Bsgamma::computeThValue(): Observable type not defined. Can be only 1 or 2");
}
