/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 * (Made by Danilo Ricciardella)
 */

#include "BClnuObservables.h"
#include "BClnu.h"
#include "StandardModel.h"

Q2moments_BClnu::Q2moments_BClnu(const StandardModel& SM_i, QCD::meson meson_i, unsigned int typ_i) 
: ThObservable(SM_i)
{
    meson1 = meson_i;
    typ = typ_i;

    setParametersForObservable(SM.getFlavour().getBClnu(meson1).initializeBClnuParameters());
}

double Q2moments_BClnu::computeThValue()
{
    double q2cut = getBinMin();
    
    if (typ == 1) return (SM.getFlavour().getBClnu(meson1).Q2moment_1(q2cut));
    else if (typ == 2) return (SM.getFlavour().getBClnu(meson1).Q2moment_2(q2cut));
    else if (typ == 3) return (SM.getFlavour().getBClnu(meson1).Q2moment_3(q2cut));
    else throw std::runtime_error("BClnuObservables::Q2moments_BClnu: incorrect type");
}

Elmoments_BClnu::Elmoments_BClnu(const StandardModel& SM_i, QCD::meson meson_i, unsigned int typ_i) 
: ThObservable(SM_i)
{
    meson1 = meson_i;
    typ = typ_i;

    setParametersForObservable(SM.getFlavour().getBClnu(meson1).initializeBClnuParameters());
}

double Elmoments_BClnu::computeThValue()
{
    double elcut = getBinMin();

    if (typ == 1) return (SM.getFlavour().getBClnu(meson1).Elmoment_1(elcut));
    else if (typ == 2) return (SM.getFlavour().getBClnu(meson1).Elmoment_2(elcut));
    else if (typ == 3) return (SM.getFlavour().getBClnu(meson1).Elmoment_3(elcut));
    else throw std::runtime_error("BClnuObservables::Elmoments_BClnu: incorrect type");
}

MXmoments_BClnu::MXmoments_BClnu(const StandardModel& SM_i, QCD::meson meson_i, unsigned int typ_i) 
: ThObservable(SM_i)
{
    meson1 = meson_i;
    typ = typ_i;

    setParametersForObservable(SM.getFlavour().getBClnu(meson1).initializeBClnuParameters());
}

double MXmoments_BClnu::computeThValue()
{
   double elcut = getBinMin();

    if (typ == 1) return (SM.getFlavour().getBClnu(meson1).Mxmoment_1(elcut));
    else if (typ == 2) return (SM.getFlavour().getBClnu(meson1).Mxmoment_2(elcut));
    else if (typ == 3) return (SM.getFlavour().getBClnu(meson1).Mxmoment_3(elcut));
    else throw std::runtime_error("BClnuObservables::Mxmoments_BClnu: incorrect type");
}

PartialAverageBR_BClnu::PartialAverageBR_BClnu(const StandardModel& SM_i, QCD::meson meson_i)
: ThObservable(SM_i)
{
    meson1 = meson_i;

    setParametersForObservable(SM.getFlavour().getBClnu(meson1).initializeBClnuParameters());
}

double PartialAverageBR_BClnu::computeThValue()
{
    double elcut = getBinMin();

    return ((SM.getFlavour().getBClnu(meson1).XDeltaBR(elcut) / SM.getFlavour().getBClnu(meson1).XDeltaBR(0.0)) * SM.getFlavour().getBClnu(meson1).getBR());
}

Vcb::Vcb(const StandardModel& SM_i, QCD::meson meson_i) 
: ThObservable(SM_i)
{
    meson1 = meson_i;

    setParametersForObservable(SM.getFlavour().getBClnu(meson1).initializeBClnuParameters());
}

double Vcb::computeThValue()
{
    double hbar = 6.582119569e-25; // Gev s
    double tauBzero = 1519e-15; // s 
    double tauBplus = 1638e-15 ; // s 

    return (sqrt(SM.getFlavour().getBClnu(meson1).getBR() / (SM.getFlavour().getBClnu(meson1).XGamma() * SM.getFlavour().getBClnu(meson1).getAmplsqfactor() * (tauBzero + tauBplus) / 2) * hbar));
}
