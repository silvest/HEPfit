/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <stdexcept>
#include "EWSMThreeLoopEW.h"

EWSMThreeLoopEW::EWSMThreeLoopEW(const EWSMcache& cache_i)
: cache(cache_i)
{
}


////////////////////////////////////////////////////////////////////////

double EWSMThreeLoopEW::DeltaAlpha_l(const double s) const
{
    double log_l[3];
    if (s == cache.getSM().getMz() * cache.getSM().getMz()) {
        log_l[0] = 2.0 * cache.logMZtoME();
        log_l[1] = 2.0 * cache.logMZtoMMU();
        log_l[2] = 2.0 * cache.logMZtoMTAU();
    } else {
        log_l[0] = log(s / cache.mf(cache.getSM().getLeptons(StandardModel::ELECTRON)) / cache.mf(cache.getSM().getLeptons(StandardModel::ELECTRON)));
        log_l[1] = log(s / cache.mf(cache.getSM().getLeptons(StandardModel::MU)) / cache.mf(cache.getSM().getLeptons(StandardModel::MU)));
        log_l[2] = log(s / cache.mf(cache.getSM().getLeptons(StandardModel::TAU)) / cache.mf(cache.getSM().getLeptons(StandardModel::TAU)));
    }

    double threeLoop[3];
    for (int i = 0; i < 3; i++) {
        threeLoop[i] = -121.0 / 48.0 + (-5.0 + 8.0 * cache.getLog2()) * cache.getZeta2()
                - 99.0 / 16.0 * cache.getZeta3() + 10.0 * cache.getZeta5()
                + log_l[i] / 8.0;
        for (int j = 0; j < 3; j++) {
            if (i > j) { /* Pi^{(2)}_l */
                threeLoop[i] += -116.0 / 27.0 + 4.0 / 3.0 * cache.getZeta2()
                        + 38.0 / 9.0 * cache.getZeta3() + 14.0 / 9.0 * log_l[i]
                        + (5.0 / 18.0 - 4.0 / 3.0 * cache.getZeta3()) * log_l[j]
                        + log_l[i] * log_l[i] / 6.0
                        - log_l[i] * log_l[j] / 3.0;
                threeLoop[i] += M_PI * M_PI / 6.0;
            } else if (i == j) { /* Pi^{(2)}_F */
                threeLoop[i] += -307.0 / 216.0 - 8.0 / 3.0 * cache.getZeta2()
                        + 545.0 / 144.0 * cache.getZeta3()
                        + (11.0 / 6.0 - 4.0 / 3.0 * cache.getZeta3()) * log_l[i]
                        - log_l[i] * log_l[i] / 6.0;
                threeLoop[i] += M_PI * M_PI / 6.0;
            } else { /* Pi^{(2)}_h */
                threeLoop[i] += -37.0 / 6.0 + 38.0 / 9.0 * cache.getZeta3()
                        + (11.0 / 6.0 - 4.0 / 3.0 * cache.getZeta3()) * log_l[j]
                        - log_l[j] * log_l[j] / 6.0;
                threeLoop[i] += M_PI * M_PI / 6.0;
            }
        }
        threeLoop[i] /= -4.0;
    }

    return ( pow(cache.getSM().getAle() / M_PI, 3.0)
            *(threeLoop[0] + threeLoop[1] + threeLoop[2]));
}

double EWSMThreeLoopEW::DeltaAlpha_t(const double s) const
{
    return (0.0);
}

double EWSMThreeLoopEW::DeltaRho(const double Mw_i) const
{
    double Mw = Mw_i;
    double mh = cache.getSM().getMHl();
    double Mt = cache.getSM().getMtpole();
    double DeltaRho;
    if (mh == 0.0) {
        DeltaRho = 3.0 * (68.0 + 729.0 * cache.getS2() + 36.0 * cache.getD3()
                + 96.0 * cache.getZeta2() * cache.getLog2()
                + 6.0 * cache.getZeta2() - 612.0 * cache.getZeta3()
                + 324.0 * cache.getZeta4() - 72.0 * cache.getB4()
                + 3.0 * (-6572.0 / 15.0 - 4374.0 / 5.0 * cache.getS2()
                + 1472.0 / 15.0 * cache.getZeta2()
                + 440.0 * cache.getZeta3()));
    } else if (mh > 0.0 && mh <= 2.5 * Mt) {
        double delta = mh / Mt - 1.0;
        DeltaRho = 95.92 - 111.98 * delta + 8.099 * delta * delta
                + 9.36 * pow(delta, 3.0) + 7.27 * pow(delta, 4.0)
                - 15.60 * pow(delta, 5.0);
    } else if (mh > 2.5 * Mt) {
        double Y = 4.0 * pow(Mt / mh, 2.0);
        double logY = 2.0 * (cache.getLog2() + cache.logMTOPtoMH());
        double logY2 = logY*logY;
        double logY3 = logY2*logY;
        DeltaRho = 1.0 / Y * (-3.17 - 83.25 * logY)
                - 189.93 - 231.48 * logY - 142.06 * logY2 + 2.75 * logY3
                + Y * (-332.34 + 77.71 * logY - 68.67 * logY2 + 51.79 * logY3)
                + Y * Y * (227.55 - 510.55 * logY + 87.77 * logY2 + 6.41 * logY3)
                + Y * Y * Y * (-58.40 - 329.18 * logY + 20.42 * logY2 + 14.54 * logY3)
                + Y * Y * Y * Y * (-36.14 - 381.88 * logY + 18.63 * logY2 + 15.04 * logY3)
                + Y * Y * Y * Y * Y * (-39.08 - 416.36 * logY + 13.76 * logY2 + 17.19 * logY3);
    } else {
        throw std::runtime_error("Higgs mass is out of range in EWSMThreeLoopEW::DeltaRho()");

    }
    DeltaRho *= pow(cache.Xt_alpha(Mw), 3.0);
    return DeltaRho;
}

double EWSMThreeLoopEW::DeltaR_rem(const double Mw_i) const
{
    return (0.0);
}

gslpp::complex EWSMThreeLoopEW::deltaRho_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    return ( gslpp::complex(0.0, 0.0, false));
}

gslpp::complex EWSMThreeLoopEW::deltaKappa_rem_f(const Particle f, const double Mw_i) const
{
    if (f.is("TOP")) return ( gslpp::complex(0.0, 0.0, false));
    return ( gslpp::complex(0.0, 0.0, false));
}




