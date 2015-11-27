/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "EWSMThreeLoopEW2QCD.h"

EWSMThreeLoopEW2QCD::EWSMThreeLoopEW2QCD(const EWSMcache& cache_i)
: cache(cache_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWSMThreeLoopEW2QCD::DeltaAlpha_l(const double s) const
{
    return (0.0);
}

double EWSMThreeLoopEW2QCD::DeltaAlpha_t(const double s) const
{
    return (0.0);
}

double EWSMThreeLoopEW2QCD::DeltaRho(const double Mw_i) const
{
    double Mw = Mw_i;
    double mh = cache.getSM().getMHl();
    double Mt = cache.getSM().getMtpole();
    double DeltaRho;
    if (mh == 0.0) {
        DeltaRho = 4.0 * (185.0 / 3.0 + 729.0 / 4.0 * cache.getS2()
                - 48.0 * cache.getZeta2() * cache.getLog2()
                - 151.0 / 6.0 * cache.getZeta2() + 29.0 * cache.getZeta3()
                - 24.0 * cache.getZeta4() + 12.0 * cache.getB4());
    } else if (mh > 0.0 && mh <= 2.5 * Mt) {
        double delta = mh / Mt - 1.0;
        DeltaRho = 157.295 + 112.00 * delta - 24.73 * delta * delta
                + 7.39 * pow(delta, 3.0) - 3.52 * pow(delta, 4.0)
                + 2.06 * pow(delta, 5.0);
    } else if (mh > 2.5 * Mt) {
        double Y = 4.0 * pow(Mt / mh, 2.0);
        double logY = 2.0 * (cache.getLog2() + cache.logMTOPtoMH());
        double logY2 = logY*logY;
        double logY3 = logY2*logY;
        DeltaRho = 79.73 - 47.77 * logY + 42.07 * logY2 + 9.00 * logY3
                + Y * (225.16 - 179.74 * logY + 70.22 * logY2 - 19.22 * logY3)
                + Y * Y * (-76.07 + 25.33 * logY - 9.17 * logY2 - 5.57 * logY3)
                + Y * Y * Y * (-10.10 - 24.69 * logY - 0.30 * logY2 - 5.46 * logY3)
                + Y * Y * Y * Y * (-4.52 - 32.85 * logY + 0.72 * logY2 - 5.25 * logY3)
                + Y * Y * Y * Y * Y * (-2.55 - 36.61 * logY + 1.06 * logY2 - 5.14 * logY3);
    } else {
        throw std::runtime_error("Higgs mass is out of range in EWSMThreeLoopEW2QCD::DeltaRho()");
    }
    DeltaRho *= pow(cache.Xt_alpha(Mw), 2.0) * cache.alsMt() / M_PI;
    return DeltaRho;
}

double EWSMThreeLoopEW2QCD::DeltaR_rem(const double Mw_i) const
{
    return (0.0);
}

gslpp::complex EWSMThreeLoopEW2QCD::deltaRho_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    return ( gslpp::complex(0.0, 0.0, false));
}

gslpp::complex EWSMThreeLoopEW2QCD::deltaKappa_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    return ( gslpp::complex(0.0, 0.0, false));
}




