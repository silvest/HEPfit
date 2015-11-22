/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "EWSMThreeLoopQCD.h"

EWSMThreeLoopQCD::EWSMThreeLoopQCD(const EWSMcache& cache_i)
: cache(cache_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWSMThreeLoopQCD::DeltaAlpha_l(const double s) const
{
    return (0.0);
}

double EWSMThreeLoopQCD::DeltaAlpha_t(const double s) const
{
    double xt = s / cache.getSM().getMtpole() / cache.getSM().getMtpole();
    double log_t, als;
    if (s == cache.getSM().getMz() * cache.getSM().getMz()) {
        log_t = 2.0 * cache.logMZtoMTOP();
        als = cache.getSM().getAlsMz();
    } else {
        double mu = sqrt(s);
        log_t = log(s / pow(cache.mf(cache.getSM().getQuarks(QCD::TOP), mu), 2.0));
        als = cache.Als(sqrt(s), FULLNNLO);
    }
    double tmp = ((28.220 + 9.702 * log_t)
            + xt * (6.924 + 1.594 * log_t))
            * pow(als / M_PI, 2.0);
    tmp *= -4.0 / 45.0 * cache.getSM().getAle() / M_PI*xt;
    return tmp;
}

double EWSMThreeLoopQCD::DeltaRho(const double Mw_i) const
{
    double Mw = Mw_i;
    return ( 3.0 * cache.Xt_alpha(Mw) * pow(cache.alsMt() / M_PI, 2.0) * deltaQCD_3(Mw));
}

double EWSMThreeLoopQCD::DeltaR_rem(const double Mw_i) const
{
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double cW2 = cache.getSM().cW2(Mw);

    /* Logarithm */
    double log_cW2 = cache.log_cW2(Mw);

    // O(alpha_s) correction to Delta r^{ud} of O(alpha alpha_s). 
    double DeltaR;
    DeltaR = -log_cW2;
    DeltaR *= (cW2 - sW2) / 4.0 / sW2 / sW2;
    DeltaR *= cache.getSM().getAle() * cache.getSM().getAlsMz() / M_PI / M_PI;
    DeltaR *= 1.4097 * cache.getSM().getAlsMz() / M_PI;
    return DeltaR;
}

gslpp::complex EWSMThreeLoopQCD::deltaRho_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    return ( gslpp::complex(0.0, 0.0, false));
}

gslpp::complex EWSMThreeLoopQCD::deltaKappa_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP"))
        return ( gslpp::complex(0.0, 0.0, false));
    else {
        double Mw = Mw_i;
        return ( -3.0 * cache.Xt_alpha(Mw) * cache.getSM().cW2(Mw) / cache.getSM().sW2(Mw)
                * pow(cache.alsMt() / M_PI, 2.0)
                *(deltaQCD_3(Mw) + deltaQCD_kappa3(Mw).real()));
    }
}


////////////////////////////////////////////////////////////////////////

double EWSMThreeLoopQCD::deltaQCD_3(const double Mw_i) const
{
    double dQCD3;
    double lZ = 2.0 * cache.logMZtoMTOP();
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double log2 = cache.getLog2();
    double zeta2 = cache.getZeta2();
    double zeta3 = cache.getZeta3();
    double zeta4 = cache.getZeta4();
    double S2 = cache.getS2(), D3 = cache.getD3(), B4 = cache.getB4();
    double MZtoMT = cache.getSM().getMz() / cache.getSM().getMtpole();
    double nf = 6.0;
    dQCD3 = 157.0 / 648.0 - 3313.0 / 162.0 * zeta2 - 308.0 / 27.0 * zeta3
            + 143.0 / 18.0 * zeta4 - 4.0 / 3.0 * zeta2 * log2
            + 441.0 / 8.0 * S2 - B4 / 9.0 - D3 / 18.0
            - (1.0 / 18.0 - 13.0 / 9.0 * zeta2 + 4.0 / 9.0 * zeta3) * nf
            + pow(MZtoMT, 2.0)
            *(-17.224 + 0.08829 * lZ + 0.4722 * lZ * lZ
            + (22.6367 + 1.2527 * lZ - 0.8519 * lZ * lZ) * sW2)
            + pow(MZtoMT, 4.0)
            *(-7.7781 - 0.07226 * lZ + 0.004938 * lZ * lZ
            + (21.497 + 0.05794 * lZ - 0.006584 * lZ * lZ) * sW2
            - 21.0799 * sW2 * sW2);
    return dQCD3;
}

gslpp::complex EWSMThreeLoopQCD::deltaQCD_kappa3(const double Mw_i) const
{
    gslpp::complex dQCDk3;
    double lZ = 2.0 * cache.logMZtoMTOP();
    double Mw = Mw_i;
    double sW2 = cache.getSM().sW2(Mw);
    double MZtoMT = cache.getSM().getMz() / cache.getSM().getMtpole();
    dQCDk3.real() = -deltaQCD_3(Mw)
            + pow(MZtoMT, 2.0)
            *((22.6367 + 1.2527 * lZ - 0.8519 * lZ * lZ) * sW2
            + (-11.3184 - 0.6263 * lZ + 0.4259 * lZ * lZ) * sW2)
            + pow(MZtoMT, 4.0)
            *((21.497 + 0.05794 * lZ - 0.006584 * lZ * lZ) * sW2
            + (-16.0186 - 0.02897 * lZ + 0.003292 * lZ * lZ) * sW2
            - 21.0799 * sW2 * sW2 + 10.54 * sW2 * sW2);
    dQCDk3.imag() = pow(MZtoMT, 2.0)
            *((-1.968 + 2.676 * lZ) * sW2
            + (2.6235 - 3.5682 * lZ) * sW2 * sW2)
            + pow(MZtoMT, 4.0)
            *((-0.09102 + 0.02069 * lZ) * sW2
            + (0.1214 - 0.02758 * lZ) * sW2 * sW2);
    return dQCDk3;
}







