/*
 * Copyright (C) 2019 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "HiggsThObservables.h"
#include "NPbase.h"

muggH::muggH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggH called with a class whose parent is not NPbase");
}

double muggH::computeThValue()
{
    return myNPbase->muggH(sqrt_s);
}

muVBF::muVBF(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBF called with a class whose parent is not NPbase");

}

double muVBF::computeThValue()
{
    return myNPbase->muVBF(sqrt_s);
}

muVBFgamma::muVBFgamma(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFgamma called with a class whose parent is not NPbase");

}

double muVBFgamma::computeThValue()
{
    return myNPbase->muVBFgamma(sqrt_s);
}

mueeWBF::mueeWBF(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBF called with a class whose parent is not NPbase");

}

double mueeWBF::computeThValue()
{
    return myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep);
}


mueeHvv::mueeHvv(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvv called with a class whose parent is not NPbase");
}

double mueeHvv::computeThValue()
{
    return myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep);
}


mueeZBF::mueeZBF(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZBF called with a class whose parent is not NPbase");

}

double mueeZBF::computeThValue()
{
    return myNPbase->mueeZBF(sqrt_s, Pol_em, Pol_ep);
}


muepWBF::muepWBF(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBF called with a class whose parent is not NPbase");

}

double muepWBF::computeThValue()
{
    return myNPbase->muepWBF(sqrt_s);
}

muWH::muWH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWH called with a class whose parent is not NPbase");
}

double muWH::computeThValue()
{
    return myNPbase->muWH(sqrt_s);
}

muWHpT250::muWHpT250(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHpT250 called with a class whose parent is not NPbase");
}

double muWHpT250::computeThValue()
{
    return myNPbase->muWHpT250(sqrt_s);
}

muepZBF::muepZBF(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBF called with a class whose parent is not NPbase");

}

double muepZBF::computeThValue()
{
    return myNPbase->muepZBF(sqrt_s);
}

muZH::muZH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZH called with a class whose parent is not NPbase");
}

double muZH::computeThValue()
{
    return myNPbase->muZH(sqrt_s);
}

muZHpT250::muZHpT250(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHpT250 called with a class whose parent is not NPbase");
}

double muZHpT250::computeThValue()
{
    return myNPbase->muZHpT250(sqrt_s);
}


mueeZHGen::mueeZHGen(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHGen called with a class whose parent is not NPbase");
}

double mueeZHGen::computeThValue()
{
    return ( (myNPbase->mueeZHGen(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->delta2sH3(myNPbase->C1eeZH(sqrt_s))) );
}


mueeZH::mueeZH(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZH called with a class whose parent is not NPbase");
}

double mueeZH::computeThValue()
{
    return ( (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->delta2sH3(myNPbase->C1eeZH(sqrt_s))) );
}

mueeZllH::mueeZllH(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZllH called with a class whose parent is not NPbase");
}

double mueeZllH::computeThValue()
{
    return ( (myNPbase->mueeZllH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->delta2sH3(myNPbase->C1eeZH(sqrt_s))) );
}


mueeZqqH::mueeZqqH(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZqqH called with a class whose parent is not NPbase");
}

double mueeZqqH::computeThValue()
{
    return ( (myNPbase->mueeZqqH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->delta2sH3(myNPbase->C1eeZH(sqrt_s))) );
}


aPsk::aPsk(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("aPsk called with a class whose parent is not NPbase");
}

double aPsk::computeThValue()
{
    return myNPbase->aPskPol(sqrt_s, Pol_em, Pol_ep);
}


bPsk::bPsk(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("bPsk called with a class whose parent is not NPbase");
}

double bPsk::computeThValue()
{
    return myNPbase->bPskPol(sqrt_s, Pol_em, Pol_ep);
}


muVH::muVH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVH called with a class whose parent is not NPbase");
}

double muVH::computeThValue()
{
    return myNPbase->muVH(sqrt_s);
}

muVHpT250::muVHpT250(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHpT250 called with a class whose parent is not NPbase");
}

double muVHpT250::computeThValue()
{
    return myNPbase->muVHpT250(sqrt_s);
}

muVBFpVH::muVBFpVH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFpVH called with a class whose parent is not NPbase");
}

double muVBFpVH::computeThValue()
{
    return myNPbase->muVBFpVH(sqrt_s);
}

muttH::muttH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttH called with a class whose parent is not NPbase");
}

double muttH::computeThValue()
{
    return myNPbase->muttH(sqrt_s);
}

mutHq::mutHq(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mutHq called with a class whose parent is not NPbase");
}

double mutHq::computeThValue()
{
    return myNPbase->mutHq(sqrt_s);
}

muggHpttH::muggHpttH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpttH called with a class whose parent is not NPbase");
}

double muggHpttH::computeThValue()
{
    return myNPbase->muggHpttH(sqrt_s);
}

mueettH::mueettH(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueettH called with a class whose parent is not NPbase");

}

double mueettH::computeThValue()
{
    return myNPbase->mueettH(sqrt_s, Pol_em, Pol_ep);
}


mummH::mummH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummH called with a class whose parent is not NPbase");
}

double mummH::computeThValue()
{
    return myNPbase->mummH(sqrt_s);
}

mummHNWA::mummHNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHNWA called with a class whose parent is not NPbase");
}

double mummHNWA::computeThValue()
{
    return myNPbase->mummHNWA(sqrt_s);
}

// -----------------------------------------------------------------------------
// Decay widths
// -----------------------------------------------------------------------------

GammaHtoggRatio::GammaHtoggRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtoggRatio called with a class whose parent is not NPbase");
}

double GammaHtoggRatio::computeThValue()
{
    return myNPbase->GammaHggRatio();
}

GammaHtoWWRatio::GammaHtoWWRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtoWWRatio called with a class whose parent is not NPbase");
}

double GammaHtoWWRatio::computeThValue()
{
    return myNPbase->GammaHWWRatio();
}

GammaHtoZZRatio::GammaHtoZZRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtoZZRatio called with a class whose parent is not NPbase");
}

double GammaHtoZZRatio::computeThValue()
{
    return myNPbase->GammaHZZRatio();
}

GammaHtoZgaRatio::GammaHtoZgaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtoZgaRatio called with a class whose parent is not NPbase");
}

double GammaHtoZgaRatio::computeThValue()
{
    return myNPbase->GammaHZgaRatio();
}

GammaHtogagaRatio::GammaHtogagaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtogagaRatio called with a class whose parent is not NPbase");
}

double GammaHtogagaRatio::computeThValue()
{
    return myNPbase->GammaHgagaRatio();
}

GammaHtomumuRatio::GammaHtomumuRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtomumuRatio called with a class whose parent is not NPbase");
}

double GammaHtomumuRatio::computeThValue()
{
    return myNPbase->GammaHmumuRatio();
}

GammaHtotautauRatio::GammaHtotautauRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtotautauRatio called with a class whose parent is not NPbase");
}

double GammaHtotautauRatio::computeThValue()
{
    return myNPbase->GammaHtautauRatio();
}

GammaHtossRatio::GammaHtossRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtossRatio called with a class whose parent is not NPbase");
}

double GammaHtossRatio::computeThValue()
{
    return myNPbase->GammaHssRatio();
}

GammaHtoccRatio::GammaHtoccRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtoccRatio called with a class whose parent is not NPbase");
}

double GammaHtoccRatio::computeThValue()
{
    return myNPbase->GammaHccRatio();
}

GammaHtobbRatio::GammaHtobbRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHtobbRatio called with a class whose parent is not NPbase");
}

double GammaHtobbRatio::computeThValue()
{
    return myNPbase->GammaHbbRatio();
}

GammaHRatio::GammaHRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("GammaHRatio called with a class whose parent is not NPbase");
}

double GammaHRatio::computeThValue()
{
    return myNPbase->computeGammaTotalRatio();
}

// -----------------------------------------------------------------------------
// Branching Ratios
// -----------------------------------------------------------------------------

BrHtoinvRatio::BrHtoinvRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoinvRatio called with a class whose parent is not NPbase");
}

double BrHtoinvRatio::computeThValue()
{
    return myNPbase->BrHtoinvRatio();
}

BrHinvisible::BrHinvisible(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHinvisible called with a class whose parent is not NPbase");
}

double BrHinvisible::computeThValue()
{
    return myNPbase->Br_H_inv();
}

BrHinvisibleNP::BrHinvisibleNP(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHinvisibleNP called with a class whose parent is not NPbase");
}

double BrHinvisibleNP::computeThValue()
{
    return myNPbase->Br_H_inv_NP();
}

BrHexotic::BrHexotic(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHexotic called with a class whose parent is not NPbase");
}

double BrHexotic::computeThValue()
{
    return myNPbase->Br_H_exo();
}

BrHtovisRatio::BrHtovisRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtovisRatio called with a class whose parent is not NPbase");
}

double BrHtovisRatio::computeThValue()
{
    return myNPbase->BrHvisRatio();
}

BrHtoggRatio::BrHtoggRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoggRatio called with a class whose parent is not NPbase");
}

double BrHtoggRatio::computeThValue()
{
    return myNPbase->BrHggRatio();
}

BrHtoWWRatio::BrHtoWWRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoWWRatio called with a class whose parent is not NPbase");
}

double BrHtoWWRatio::computeThValue()
{
    return myNPbase->BrHWWRatio();
}

BrHtoZZRatio::BrHtoZZRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZZRatio called with a class whose parent is not NPbase");
}

double BrHtoZZRatio::computeThValue()
{
    return myNPbase->BrHZZRatio();
}

BrHtoVVRatio::BrHtoVVRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoVVRatio called with a class whose parent is not NPbase");
}

double BrHtoVVRatio::computeThValue()
{
    return myNPbase->BrHVVRatio();
}

BrHtoZgaRatio::BrHtoZgaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZgaRatio called with a class whose parent is not NPbase");
}

double BrHtoZgaRatio::computeThValue()
{
    return myNPbase->BrHZgaRatio();
}

BrHtoZgallRatio::BrHtoZgallRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZgallRatio called with a class whose parent is not NPbase");
}

double BrHtoZgallRatio::computeThValue()
{
    return myNPbase->BrHZgallRatio();
}

BrHtoZgaeeRatio::BrHtoZgaeeRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZgaeeRatio called with a class whose parent is not NPbase");
}

double BrHtoZgaeeRatio::computeThValue()
{
    return myNPbase->BrHZgaeeRatio();
}

BrHtoZgamumuRatio::BrHtoZgamumuRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZgamumuRatio called with a class whose parent is not NPbase");
}

double BrHtoZgamumuRatio::computeThValue()
{
    return myNPbase->BrHZgamumuRatio();
}

BrHtogagaRatio::BrHtogagaRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogagaRatio called with a class whose parent is not NPbase");
}

double BrHtogagaRatio::computeThValue()
{
    return myNPbase->BrHgagaRatio();
}

BrHtomumuRatio::BrHtomumuRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtomumuRatio called with a class whose parent is not NPbase");
}

double BrHtomumuRatio::computeThValue()
{
    return myNPbase->BrHmumuRatio();
}

BrHtotautauRatio::BrHtotautauRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtotautauRatio called with a class whose parent is not NPbase");
}

double BrHtotautauRatio::computeThValue()
{
    return myNPbase->BrHtautauRatio();
}

BrHtoccRatio::BrHtoccRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoccRatio called with a class whose parent is not NPbase");
}

double BrHtoccRatio::computeThValue()
{
    return myNPbase->BrHccRatio();
}

BrHtobbRatio::BrHtobbRatio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtobbRatio called with a class whose parent is not NPbase");
}

double BrHtobbRatio::computeThValue()
{
    return myNPbase->BrHbbRatio();
}

// -----------------------------------------------------------------------------
// More 4 fermion decays
// -----------------------------------------------------------------------------


BrHto2l2vRatio::BrHto2l2vRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto2l2vRatio called with a class whose parent is not NPbase");
}

double BrHto2l2vRatio::computeThValue()
{
    return myNPbase->BrH2l2vRatio();
}


BrHtoevmuvRatio::BrHtoevmuvRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoevmuvRatio called with a class whose parent is not NPbase");
}

double BrHtoevmuvRatio::computeThValue()
{
    return myNPbase->BrHevmuvRatio();
}


BrHto2e2vRatio::BrHto2e2vRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto2e2vRatio called with a class whose parent is not NPbase");
}

double BrHto2e2vRatio::computeThValue()
{
    // SM decay widths (from MG simulations)
    double wH2e2vSM=0.93291e-06, wH2evSM=0.10152e-04;
    
    // Sum
    double wH2e2vTSM=wH2e2vSM+wH2evSM;
    
    return ( (wH2e2vSM*(myNPbase->BrH2e2vRatio()) + wH2evSM*(myNPbase->BrH2evRatio())) / wH2e2vTSM );
}


BrHto2mu2vRatio::BrHto2mu2vRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto2mu2vRatio called with a class whose parent is not NPbase");
}

double BrHto2mu2vRatio::computeThValue()
{
    // SM decay widths (from MG simulations)
    double wH2mu2vSM=0.93288e-06, wH2muvSM=0.10163e-04;
    
    // Sum
    double wH2mu2vTSM=wH2mu2vSM+wH2muvSM;
    
    return ( (wH2mu2vSM*(myNPbase->BrH2mu2vRatio()) + wH2muvSM*(myNPbase->BrH2muvRatio())) / wH2mu2vTSM );    
}


BrHto4lRatio::BrHto4lRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto4lRatio called with a class whose parent is not NPbase");
}

double BrHto4lRatio::computeThValue()
{
    return myNPbase->BrH4lRatio();
}


BrHto4eRatio::BrHto4eRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto4eRatio called with a class whose parent is not NPbase");
}

double BrHto4eRatio::computeThValue()
{
    return myNPbase->BrH4eRatio();
}


BrHto4muRatio::BrHto4muRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto4muRatio called with a class whose parent is not NPbase");
}

double BrHto4muRatio::computeThValue()
{
    return myNPbase->BrH4muRatio();
}


BrHto2e2muRatio::BrHto2e2muRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto2e2muRatio called with a class whose parent is not NPbase");
}

double BrHto2e2muRatio::computeThValue()
{
    return myNPbase->BrH2e2muRatio();
}



BrHtolljjRatio::BrHtolljjRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtolljjRatio called with a class whose parent is not NPbase");
}

double BrHtolljjRatio::computeThValue()
{
    return myNPbase->BrHlljjRatio();
}


BrHtolvjjRatio::BrHtolvjjRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtolvjjRatio called with a class whose parent is not NPbase");
}

double BrHtolvjjRatio::computeThValue()
{
    return myNPbase->BrHlvjjRatio();
}


BrHtolv_lvorjjRatio::BrHtolv_lvorjjRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtolv_lvorjjRatio called with a class whose parent is not NPbase");
}

double BrHtolv_lvorjjRatio::computeThValue()
{
    return myNPbase->BrHlv_lvorjjRatio();
}


BrHtoll_vvorjjRatio::BrHtoll_vvorjjRatio(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoll_vvorjjRatio called with a class whose parent is not NPbase");
}

double BrHtoll_vvorjjRatio::computeThValue()
{
    return myNPbase->BrHll_vvorjjRatio();
}

// -----------------------------------------------------------------------------
// Ratios of BR (ratios with SM)
// -----------------------------------------------------------------------------

BrHtogaga_over_mumu_Ratio::BrHtogaga_over_mumu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogaga_over_mumu_Ratio called with a class whose parent is not NPbase");
}

double BrHtogaga_over_mumu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHgagaRatio()) - (myNPbase->BrHmumuRatio()));
    } else {
        return (myNPbase->BrHgagaRatio()) / (myNPbase->BrHmumuRatio());
    }
}

BrHtoZga_over_mumu_Ratio::BrHtoZga_over_mumu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZga_over_mumu_Ratio called with a class whose parent is not NPbase");
}

double BrHtoZga_over_mumu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHZgaRatio()) - (myNPbase->BrHmumuRatio()));
    } else {
        return (myNPbase->BrHZgaRatio()) / (myNPbase->BrHmumuRatio());
    }
}

BrHtoZmumuga_over_mumu_Ratio::BrHtoZmumuga_over_mumu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZmumuga_over_mumu_Ratio called with a class whose parent is not NPbase");
}

double BrHtoZmumuga_over_mumu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHZgamumuRatio()) - (myNPbase->BrHmumuRatio()));
    } else {
        return (myNPbase->BrHZgamumuRatio()) / (myNPbase->BrHmumuRatio());
    }
}

BrHtoZga_over_4mu_Ratio::BrHtoZga_over_4mu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZga_over_4mu_Ratio called with a class whose parent is not NPbase");
}

double BrHtoZga_over_4mu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHZgaRatio()) - (myNPbase->BrH4muRatio()));
    } else {
        return (myNPbase->BrHZgaRatio()) / (myNPbase->BrH4muRatio());
    }
}

BrHtoZmumuga_over_4mu_Ratio::BrHtoZmumuga_over_4mu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZmumuga_over_4mu_Ratio called with a class whose parent is not NPbase");
}

double BrHtoZmumuga_over_4mu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHZgamumuRatio()) - (myNPbase->BrH4muRatio()));
    } else {
        return (myNPbase->BrHZgamumuRatio()) / (myNPbase->BrH4muRatio());
    }
}

BrHtogaga_over_4l_Ratio::BrHtogaga_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogaga_over_4l_Ratio called with a class whose parent is not NPbase");
}

double BrHtogaga_over_4l_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHgagaRatio()) - (myNPbase->BrH4lRatio()));
    } else {
        return (myNPbase->BrHgagaRatio()) / (myNPbase->BrH4lRatio());
    }
}

BrHtobb_over_4l_Ratio::BrHtobb_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtobb_over_4l_Ratio called with a class whose parent is not NPbase");
}

double BrHtobb_over_4l_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHbbRatio()) - (myNPbase->BrH4lRatio()));
    } else {
        return (myNPbase->BrHbbRatio()) / (myNPbase->BrH4lRatio());
    }
}


BrHto2l2v_over_4l_Ratio::BrHto2l2v_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto2l2v_over_4l_Ratio called with a class whose parent is not NPbase");
}

double BrHto2l2v_over_4l_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrH2l2vRatio()) - (myNPbase->BrH4lRatio()));
    } else {
        return (myNPbase->BrH2l2vRatio()) / (myNPbase->BrH4lRatio());
    }
}


BrHtotautau_over_4l_Ratio::BrHtotautau_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtotautau_over_4l_Ratio called with a class whose parent is not NPbase");
}

double BrHtotautau_over_4l_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHtautauRatio()) - (myNPbase->BrH4lRatio()));
    } else {
        return (myNPbase->BrHtautauRatio()) / (myNPbase->BrH4lRatio());
    }
}


BrHtogaga_over_2e2mu_Ratio::BrHtogaga_over_2e2mu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogaga_over_2e2mu_Ratio called with a class whose parent is not NPbase");
}

double BrHtogaga_over_2e2mu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHgagaRatio()) - (myNPbase->BrH2e2muRatio()));
    } else {
        return (myNPbase->BrHgagaRatio()) / (myNPbase->BrH2e2muRatio());
    }
}

BrHtoZga_over_4l_Ratio::BrHtoZga_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZga_over_4l_Ratio called with a class whose parent is not NPbase");
}

double BrHtoZga_over_4l_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHZgaRatio()) - (myNPbase->BrH4lRatio()));
    } else {
        return (myNPbase->BrHZgaRatio()) / (myNPbase->BrH4lRatio());
    }
}

BrHtomumu_over_4l_Ratio::BrHtomumu_over_4l_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtomumu_over_4l_Ratio called with a class whose parent is not NPbase");
}

double BrHtomumu_over_4l_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHmumuRatio()) - (myNPbase->BrH4lRatio()));
    } else {
        return (myNPbase->BrHmumuRatio()) / (myNPbase->BrH4lRatio());
    }
}

BrHtomumu_over_4mu_Ratio::BrHtomumu_over_4mu_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtomumu_over_4mu_Ratio called with a class whose parent is not NPbase");
}

double BrHtomumu_over_4mu_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHmumuRatio()) - (myNPbase->BrH4muRatio()));
    } else {
        return (myNPbase->BrHmumuRatio()) / (myNPbase->BrH4muRatio());
    }
}

BrHto4l_over_gaga_Ratio::BrHto4l_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto4l_over_gaga_Ratio called with a class whose parent is not NPbase");
}

double BrHto4l_over_gaga_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrH4lRatio()) - (myNPbase->BrHgagaRatio()));
    } else {
        return (myNPbase->BrH4lRatio()) / (myNPbase->BrHgagaRatio());
    }
}

BrHtoZga_over_gaga_Ratio::BrHtoZga_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtoZga_over_gaga_Ratio called with a class whose parent is not NPbase");
}

double BrHtoZga_over_gaga_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHZgaRatio()) - (myNPbase->BrHgagaRatio()));
    } else {
        return (myNPbase->BrHZgaRatio()) / (myNPbase->BrHgagaRatio());
    }
}

BrHtomumu_over_gaga_Ratio::BrHtomumu_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtomumu_over_gaga_Ratio called with a class whose parent is not NPbase");
}

double BrHtomumu_over_gaga_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHmumuRatio()) - (myNPbase->BrHgagaRatio()));
    } else {
        return (myNPbase->BrHmumuRatio()) / (myNPbase->BrHgagaRatio());
    }
}


BrHto2l2v_over_gaga_Ratio::BrHto2l2v_over_gaga_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHto2l2v_over_gaga_Ratio called with a class whose parent is not NPbase");
}

double BrHto2l2v_over_gaga_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrH2l2vRatio()) - (myNPbase->BrHgagaRatio()));
    } else {
        return (myNPbase->BrH2l2vRatio()) / (myNPbase->BrHgagaRatio());
    }
}


BrHtobb_over_cc_Ratio::BrHtobb_over_cc_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtobb_over_cc_Ratio called with a class whose parent is not NPbase");
}

double BrHtobb_over_cc_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHbbRatio()) - (myNPbase->BrHccRatio()));
    } else {
        return (myNPbase->BrHbbRatio()) / (myNPbase->BrHccRatio());
    }
}


BrHtogaga_over_gg_Ratio::BrHtogaga_over_gg_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogaga_over_gg_Ratio called with a class whose parent is not NPbase");
}

double BrHtogaga_over_gg_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHgagaRatio()) - (myNPbase->BrHggRatio()));
    } else {
        return (myNPbase->BrHgagaRatio()) / (myNPbase->BrHggRatio());
    }
}


BrHtogg_over_bb_Ratio::BrHtogg_over_bb_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogg_over_bb_Ratio called with a class whose parent is not NPbase");
}

double BrHtogg_over_bb_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHggRatio()) - (myNPbase->BrHbbRatio()));
    } else {
        return (myNPbase->BrHggRatio()) / (myNPbase->BrHbbRatio());
    }
}


BrHtogg_over_cc_Ratio::BrHtogg_over_cc_Ratio(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("BrHtogg_over_cc_Ratio called with a class whose parent is not NPbase");
}

double BrHtogg_over_cc_Ratio::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (1.0 + (myNPbase->BrHggRatio()) - (myNPbase->BrHccRatio()));
    } else {
        return (myNPbase->BrHggRatio()) / (myNPbase->BrHccRatio());
    }
}


// -----------------------------------------------------------------------------
// Full signal strengths (prod x decay)
// -----------------------------------------------------------------------------

muggHgaga::muggHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHgaga called with a class whose parent is not NPbase");
}

double muggHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHgagaRatio()));
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);
        return ( 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1) );       
    } else {
        return myNPbase->muggHgaga(sqrt_s);
    }
}

muggHgagaInt::muggHgagaInt(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHgagaInt called with a class whose parent is not NPbase");
}

double muggHgagaInt::computeThValue()
{
    return (myNPbase->muggHgagaInt(sqrt_s)) / (myNPbase->BrHgagaRatio());
}

muVBFHgaga::muVBFHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHgaga called with a class whose parent is not NPbase");
}

double muVBFHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHgagaRatio()));
        double muProd1  = myNPbase->delta_muVBF_1(sqrt_s);
        double muProd2  = myNPbase->delta_muVBF_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);
        
        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu;  
    } else {
        return myNPbase->muVBFHgaga(sqrt_s);
    }
}

muZHgaga::muZHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHgaga called with a class whose parent is not NPbase");
}

double muZHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrHgagaRatio()));
        double muProd1  = myNPbase->delta_muZH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muZH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);
                
        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
         
    } else {
        return myNPbase->muZHgaga(sqrt_s);
    }
}

muWHgaga::muWHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHgaga called with a class whose parent is not NPbase");
}

double muWHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHgagaRatio()));
        double muProd1  = myNPbase->delta_muWH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muWH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muWHgaga(sqrt_s);
    }
}

muVHgaga::muVHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHgaga called with a class whose parent is not NPbase");
}

double muVHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHgagaRatio()));
    } else {
        return myNPbase->muVHgaga(sqrt_s);
    }
}

muttHgaga::muttHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHgaga called with a class whose parent is not NPbase");
}

double muttHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHgagaRatio()));
        double muProd1  = myNPbase->delta_muttH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muttH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muttHgaga(sqrt_s);
    }
}

mutHgaga::mutHgaga(const StandardModel& SM_i, const double sqrt_s_i)    //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mutHgaga called with a class whose parent is not NPbase");
}

double mutHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1  = myNPbase->delta_mutH_1(sqrt_s);
        double muProd2  = myNPbase->delta_mutH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmutHgaga = myNPbase->mutHgaga(sqrt_s);
        if(NPmutHgaga==1.0){
            return ((myNPbase->mutH(sqrt_s)) )*(myNPbase->BrHgagaRatio());
        } else{
            return NPmutHgaga;
        }
    }
}

muggHpbbH_Hgaga::muggHpbbH_Hgaga(const StandardModel& SM_i, const double sqrt_s_i)  //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpbbH_Hgaga called with a class whose parent is not NPbase");
}

double muggHpbbH_Hgaga::computeThValue()                                          //AG:added
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        //Also, the bbH is missing here, I'll leave it as it was for the 
        //moment (since bbH is really suppressed).
        double NPmuggHpbbH_Hgaga = myNPbase->muggHpbbH_Hgaga(sqrt_s);
        if (NPmuggHpbbH_Hgaga == 1.0){
            return (myNPbase->muggHgaga(sqrt_s));
        } else{
            return NPmuggHpbbH_Hgaga;
        }
    }
}

muttHptH_Hgaga::muttHptH_Hgaga(const StandardModel& SM_i, const double sqrt_s_i)//AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHptH_Hgaga called with a class whose parent is not NPbase");
}

double muttHptH_Hgaga::computeThValue()                                         //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ttH = 0.499873;
    double xsSM_tH = 0.0821;
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatHq(sqrt_s);

    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double muProd2 = ( xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double dGammaR1 = myNPbase->deltaGammaHgagaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHgagaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuttHptH_Hgaga = myNPbase->muttHptH_Hgaga(sqrt_s);
        if(NPmuttHptH_Hgaga==1.0){
        return ( xsSM_ttH*(myNPbase->muttH(sqrt_s)) + xsSM_tH*(myNPbase->mutH(sqrt_s)) )/(xsSM_ttH+xsSM_tH) * (myNPbase->BrHgagaRatio()) ;
        } else {
            return NPmuttHptH_Hgaga;
        }
    }
}

muggHZga::muggHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHZga called with a class whose parent is not NPbase");
}

double muggHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHZgaRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZgaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZgaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muggHZga(sqrt_s);
    }
}

muggHZgamumu::muggHZgamumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHZgamumu called with a class whose parent is not NPbase");
}

double muggHZgamumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHZgamumuRatio()));
    } else {
        return (myNPbase->muggH(sqrt_s))*(myNPbase->BrHZgamumuRatio());
    }
}

muVBFHZga::muVBFHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHZga called with a class whose parent is not NPbase");
}

double muVBFHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHZgaRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muVBF_1(sqrt_s);
        double muProd2 = myNPbase->delta_muVBF_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZgaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZgaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muVBFHZga(sqrt_s);
    }
}

muZHZga::muZHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHZga called with a class whose parent is not NPbase");
}

double muZHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrHZgaRatio()));
    } else {
        return myNPbase->muZHZga(sqrt_s);
    }
}

muWHZga::muWHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHZga called with a class whose parent is not NPbase");
}

double muWHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHZgaRatio()));
    } else {
        return myNPbase->muWHZga(sqrt_s);
    }
}

muVHZga::muVHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHZga called with a class whose parent is not NPbase");
}

double muVHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHZgaRatio()));
    } else {
        return myNPbase->muVHZga(sqrt_s);
    }
}

muttHZga::muttHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHZga called with a class whose parent is not NPbase");
}

double muttHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHZgaRatio()));
    } else {
        return myNPbase->muttHZga(sqrt_s);
    }
}





muggHpttHptHpbbH_HZga::muggHpttHptHpbbH_HZga(const StandardModel& SM_i, const double sqrt_s_i)    //VM:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpttHptHpbbH_HZga called with a class whose parent is not NPbase");
}

double muggHpttHptHpbbH_HZga::computeThValue()                                           //VM:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this. Furthermore, the bbH is not included. In the SM this is
    //very suppressed (and it's probably also the case in the SMEFT) but 
    //in some NP models it may not be the case.Unfortunately, bbH is not
    //obtained in the SMEFT, we could just set delta_mubbH to zero in the
    //SMEFT and add here the general expression.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ggHbbH = 44.745;
    double xsSM_ttH    = 0.4998;
    double xsSM_tH     = 0.084769;
    //double xsSM_ggH = myNPbase->computeSigmaggH(sqrt_s);
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatH(sqrt_s);
    //double xsSM_bbH = myNPbase->computeSigmabbH(sqrt_s);

    
    
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ggHbbH*(myNPbase->delta_muggH_1(sqrt_s)) + xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ggHbbH+xsSM_ttH+xsSM_tH) ;
        double muProd2 = ( xsSM_ggHbbH*(myNPbase->delta_muggH_2(sqrt_s)) + xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ggHbbH+xsSM_ttH+xsSM_tH) ;
        double dGammaR1 = myNPbase->deltaGammaHZgaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZgaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end). Furthermore, the Hbb
        //is not added in the original formula, fine for the SM (probably also
        //for the SMEFT) but not for all NP models.
        double NPmuggHpttHptHpbbH_HZga = myNPbase->muggHpttHptHpbbH_HZga(sqrt_s);
        if(NPmuggHpttHptHpbbH_HZga==1.0){
            return ( (xsSM_ggHbbH*(myNPbase->muggH(sqrt_s))+xsSM_ttH*(myNPbase->muttH(sqrt_s))+xsSM_tH*(myNPbase->mutH(sqrt_s))) / (xsSM_ggHbbH+xsSM_ttH+xsSM_tH)  )*(myNPbase->BrHZgaRatio()) ;
        } else {
            return NPmuggHpttHptHpbbH_HZga;
        }  
    }
}


muVBFpVH_HZga::muVBFpVH_HZga(const StandardModel& SM_i, const double sqrt_s_i)  //VM:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFpVH_HZga called with a class whose parent is not NPbase");
}

double muVBFpVH_HZga::computeThValue()                                          //VM:added
{
    
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this. 
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_VBF = 3.49948;
    double xsSM_WH  = 1.21539;
    double xsSM_ZH  = 0.795910;
    //double xsSM_VBF = myNPbase->computeSigmaVBF(sqrt_s);
    //double xsSM_WH = myNPbase->computeSigmaWH(sqrt_s);
    //double xsSM_ZH = myNPbase->computeSigmaZH(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_VBF*(myNPbase->delta_muVBF_1(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_1(sqrt_s)) + xsSM_ZH*(myNPbase->delta_muZH_1(sqrt_s)) )/(xsSM_VBF+xsSM_WH+xsSM_ZH);
        double muProd2 = ( xsSM_VBF*(myNPbase->delta_muVBF_2(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_2(sqrt_s)) + xsSM_ZH*(myNPbase->delta_muZH_2(sqrt_s)) )/(xsSM_VBF+xsSM_WH+xsSM_ZH);
        double dGammaR1 = myNPbase->deltaGammaHZgaRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZgaRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end).
        double NPmuVBFpVH_HZga = myNPbase->muVBFpVH_HZga(sqrt_s);
        if(NPmuVBFpVH_HZga==1.0){
            return ( (xsSM_VBF*(myNPbase->muVBF(sqrt_s))+xsSM_WH*(myNPbase->muWH(sqrt_s))+xsSM_ZH*(myNPbase->muZH(sqrt_s))) / (xsSM_VBF+xsSM_WH+xsSM_ZH)  )*(myNPbase->BrHZgaRatio()) ;
        } else {
            return NPmuVBFpVH_HZga;
        } 
    }
}






muggHZZ::muggHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHZZ called with a class whose parent is not NPbase");
}

double muggHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHZZRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muggHZZ(sqrt_s);
    }
}

muVBFHZZ::muVBFHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHZZ called with a class whose parent is not NPbase");
}

double muVBFHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHZZRatio()));
        double muProd1  = myNPbase->delta_muVBF_1(sqrt_s);
        double muProd2  = myNPbase->delta_muVBF_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muVBFHZZ(sqrt_s);
    }
}

muZHZZ::muZHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHZZ called with a class whose parent is not NPbase");
}

double muZHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrHZZRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muZH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muZH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muZHZZ(sqrt_s);
    }
}

muWHZZ::muWHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHZZ called with a class whose parent is not NPbase");
}

double muWHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHZZRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muWH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muWH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muWHZZ(sqrt_s);
    }
}

muVHZZ::muVHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHZZ called with a class whose parent is not NPbase");
}

double muVHZZ::computeThValue()
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_WH  = 1.21539;
    double xsSM_ZH  = 0.795910;
    //double xsSM_WH = myNPbase->computeSigmaWH(sqrt_s);
    //double xsSM_ZH = myNPbase->computeSigmaZH(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHZZRatio()));    
        double muProd1 = ( xsSM_ZH*(myNPbase->delta_muZH_1(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_1(sqrt_s)) )/(xsSM_ZH+xsSM_WH);     
        double muProd2 = ( xsSM_ZH*(myNPbase->delta_muZH_2(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_2(sqrt_s)) )/(xsSM_ZH+xsSM_WH);     
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //return myNPbase->muVHZZ(sqrt_s);
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuVHZZ = myNPbase->muVHZZ(sqrt_s);
        if(NPmuVHZZ==1.0){
            return ( xsSM_ZH*(myNPbase->muZH(sqrt_s)) + xsSM_WH*(myNPbase->muWH(sqrt_s)) )/(xsSM_ZH+xsSM_WH) * (myNPbase->BrHZZRatio());
        } else {
            return NPmuVHZZ;
        }
    }
}

muttHZZ::muttHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHZZ called with a class whose parent is not NPbase");
}

double muttHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHZZRatio()));
    } else {
        return myNPbase->muttHZZ(sqrt_s);
    }
}

muttHptH_HZZ::muttHptH_HZZ(const StandardModel& SM_i, const double sqrt_s_i)    //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHptH_HZZ called with a class whose parent is not NPbase");
}

double muttHptH_HZZ::computeThValue()                                           //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ttH = 0.499873;
    double xsSM_tH = 0.0821;
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatHq(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( 1.0 + muProd1 + (myNPbase->BrHZZRatio()-1.) );
        double muProd1 = ( xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ttH+xsSM_tH);     
        double muProd2 = ( xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ttH+xsSM_tH);     
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuttHptH_HZZ = myNPbase->muttHptH_HZZ(sqrt_s);
        if(NPmuttHptH_HZZ==1.0){
            return ( xsSM_ttH*(myNPbase->muttH(sqrt_s)) + xsSM_tH*(myNPbase->mutH(sqrt_s)) )/(xsSM_ttH+xsSM_tH) * (myNPbase->BrHZZRatio()) ;
        } else {
            return NPmuttHptH_HZZ;
        }  
    }
}

muttHptH_Hmumu::muttHptH_Hmumu(const StandardModel& SM_i, const double sqrt_s_i)//AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHptH_Hmumu called with a class whose parent is not NPbase");
}

double muttHptH_Hmumu::computeThValue()                                         //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ttH = 0.499873;
    double xsSM_tH = 0.0821;
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatHq(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( 1.0 
        //        + ( xsSM_ttH*(myNPbase->muttH(sqrt_s)-1.) + xsSM_tH*(myNPbase->mutH(sqrt_s)-1.) )/(xsSM_ttH+xsSM_tH)
        //        + (myNPbase->BrHmumuRatio()-1.));
        double muProd1 = ( xsSM_ttH*(myNPbase->muttH(sqrt_s)-1.) + xsSM_tH*(myNPbase->mutH(sqrt_s)-1.) )/(xsSM_ttH+xsSM_tH);
        double muProd2 = ( xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuttHptH_Hmumu = myNPbase->muttHptH_Hmumu(sqrt_s);
        if(NPmuttHptH_Hmumu==1.0){
            return ( xsSM_ttH*(myNPbase->muttH(sqrt_s)) + xsSM_tH*(myNPbase->mutH(sqrt_s)) )/(xsSM_ttH+xsSM_tH) * (myNPbase->BrHmumuRatio());
        } else {
            return NPmuttHptH_Hmumu;
        }  
    }
}

muggHpbbH_HZZ::muggHpbbH_HZZ(const StandardModel& SM_i, const double sqrt_s_i)  //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpbbH_HZZ called with a class whose parent is not NPbase");
}

double muggHpbbH_HZZ::computeThValue()                                          //AG:added
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( 1.0 + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHZZRatio()-1.));
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHZZRatio1();
        double dGammaR2 = myNPbase->deltaGammaHZZRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        //Also, the bbH is missing here, I'll leave it as it was for the 
        //moment (since bbH is really suppressed in the SM).
        double NPmuggHpbbH_HZZ = myNPbase->muggHpbbH_HZZ(sqrt_s);
        if (NPmuggHpbbH_HZZ == 1.0){
            return (myNPbase->muggHZZ(sqrt_s));
        } else{
            return NPmuggHpbbH_HZZ;
        } 
    }
}

muggHZZ4l::muggHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHZZ4l called with a class whose parent is not NPbase");
}

double muggHZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrH4lRatio()));
    } else {
        return myNPbase->muggHZZ4l(sqrt_s);
    }
}

muggHZZ4mu::muggHZZ4mu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHZZ4mu called with a class whose parent is not NPbase");
}

double muggHZZ4mu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrH4muRatio()));
    } else {
        return (myNPbase->muggH(sqrt_s))*(myNPbase->BrH4muRatio());
    }
}

muVBFHZZ4l::muVBFHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHZZ4l called with a class whose parent is not NPbase");
}

double muVBFHZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrH4lRatio()));
    } else {
        return myNPbase->muVBFHZZ4l(sqrt_s);
    }
}

muZHZZ4l::muZHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHZZ4l called with a class whose parent is not NPbase");
}

double muZHZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrH4lRatio()));
    } else {
        return myNPbase->muZHZZ4l(sqrt_s);
    }
}

muWHZZ4l::muWHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHZZ4l called with a class whose parent is not NPbase");
}

double muWHZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrH4lRatio()));
    } else {
        return myNPbase->muWHZZ4l(sqrt_s);
    }
}

muVHZZ4l::muVHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHZZ4l called with a class whose parent is not NPbase");
}

double muVHZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrH4lRatio()));
    } else {
        return myNPbase->muVHZZ4l(sqrt_s);
    }
}

muttHZZ4l::muttHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHZZ4l called with a class whose parent is not NPbase");
}

double muttHZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrH4lRatio()));
    } else {
        return myNPbase->muttHZZ4l(sqrt_s);
    }
}

muggHWW::muggHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHWW called with a class whose parent is not NPbase");
}

double muggHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHWWRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHWWRatio1();
        double dGammaR2 = myNPbase->deltaGammaHWWRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muggHWW(sqrt_s);
    }
}

muVBFHWW::muVBFHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHWW called with a class whose parent is not NPbase");
}

double muVBFHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHWWRatio()));
        double muProd1  = myNPbase->delta_muVBF_1(sqrt_s);
        double muProd2  = myNPbase->delta_muVBF_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHWWRatio1();
        double dGammaR2 = myNPbase->deltaGammaHWWRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muVBFHWW(sqrt_s);
    }
}

muZHWW::muZHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHWW called with a class whose parent is not NPbase");
}

double muZHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrHWWRatio()));
        double muProd1  = myNPbase->delta_muZH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muZH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHWWRatio1();
        double dGammaR2 = myNPbase->deltaGammaHWWRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muZHWW(sqrt_s);
    }
}

muWHWW::muWHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHWW called with a class whose parent is not NPbase");
}

double muWHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHWWRatio()));
        double muProd1  = myNPbase->delta_muWH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muWH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHWWRatio1();
        double dGammaR2 = myNPbase->deltaGammaHWWRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muWHWW(sqrt_s);
    }
}

muVHWW::muVHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHWW called with a class whose parent is not NPbase");
}

double muVHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHWWRatio()));
    } else {
        return myNPbase->muVHWW(sqrt_s);
    }
}

muttHWW::muttHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHWW called with a class whose parent is not NPbase");
}

double muttHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHWWRatio()));
    } else {
        return myNPbase->muttHWW(sqrt_s);
    }
}

muttHptH_HWW::muttHptH_HWW(const StandardModel& SM_i, const double sqrt_s_i)    //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHptH_HWW called with a class whose parent is not NPbase");
}

double muttHptH_HWW::computeThValue()                                           //AG:added
{
    
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ttH = 0.499873;
    double xsSM_tH = 0.0821;
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatHq(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double muProd2 = ( xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double dGammaR1 = myNPbase->deltaGammaHWWRatio1();
        double dGammaR2 = myNPbase->deltaGammaHWWRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuttHptH_HWW = myNPbase->muttHptH_HWW(sqrt_s);
        if(NPmuttHptH_HWW==1.0){
            return ( xsSM_ttH*(myNPbase->muttH(sqrt_s)) + xsSM_tH*(myNPbase->mutH(sqrt_s)) )/(xsSM_ttH+xsSM_tH) * (myNPbase->BrHWWRatio()) ;
        } else {
            return NPmuttHptH_HWW;
        }   
    }
}

muggHpbbH_HWW::muggHpbbH_HWW(const StandardModel& SM_i, const double sqrt_s_i)  //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpbbH_HWW called with a class whose parent is not NPbase");
}

double muggHpbbH_HWW::computeThValue()                                          //AG:added
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( 1.0 + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHWWRatio()-1.));
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHWWRatio1();
        double dGammaR2 = myNPbase->deltaGammaHWWRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        //Also, the bbH is missing here, I'll leave it as it was for the 
        //moment (since bbH is really suppressed in the SM).
        double NPmuggHpbbH_HWW = myNPbase->muggHpbbH_HWW(sqrt_s);
        if (NPmuggHpbbH_HWW == 1.0){
            return (myNPbase->muggHWW(sqrt_s));
        } else{
            return NPmuggHpbbH_HWW;
        } 
    }
}

muggHWW2l2v::muggHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHWW2l2v called with a class whose parent is not NPbase");
}

double muggHWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrH2l2vRatio()));
    } else {
        return myNPbase->muggHWW2l2v(sqrt_s);
    }
}

muVBFHWW2l2v::muVBFHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHWW2l2v called with a class whose parent is not NPbase");
}

double muVBFHWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrH2l2vRatio()));
    } else {
        return myNPbase->muVBFHWW2l2v(sqrt_s);
    }
}

muZHWW2l2v::muZHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHWW2l2v called with a class whose parent is not NPbase");
}

double muZHWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrH2l2vRatio()));
    } else {
        return myNPbase->muZHWW2l2v(sqrt_s);
    }
}

muWHWW2l2v::muWHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHWW2l2v called with a class whose parent is not NPbase");
}

double muWHWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrH2l2vRatio()));
    } else {
        return myNPbase->muWHWW2l2v(sqrt_s);
    }
}

muVHWW2l2v::muVHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHWW2l2v called with a class whose parent is not NPbase");
}

double muVHWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrH2l2vRatio()));
    } else {
        return myNPbase->muVHWW2l2v(sqrt_s);
    }
}

muttHWW2l2v::muttHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHWW2l2v called with a class whose parent is not NPbase");
}

double muttHWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrH2l2vRatio()));
    } else {
        return myNPbase->muttHWW2l2v(sqrt_s);
    }
}

muttHVV::muttHVV(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHVV called with a class whose parent is not NPbase");
}

double muttHVV::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHVVRatio()));
    } else {
        return (myNPbase->muttH(sqrt_s)) * (myNPbase->BrHVVRatio());
    }
}

muggHmumu::muggHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHmumu called with a class whose parent is not NPbase");
}

double muggHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHmumuRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muggHmumu(sqrt_s);
    }
}

muVBFHmumu::muVBFHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHmumu called with a class whose parent is not NPbase");
}

double muVBFHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHmumuRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muVBF_1(sqrt_s);
        double muProd2 = myNPbase->delta_muVBF_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muVBFHmumu(sqrt_s);
    }
}

muZHmumu::muZHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHmumu called with a class whose parent is not NPbase");
}

double muZHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = myNPbase->delta_muZH_1(sqrt_s);
        double muProd2 = myNPbase->delta_muZH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muZHmumu(sqrt_s);
    }
}

muWHmumu::muWHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHmumu called with a class whose parent is not NPbase");
}

double muWHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHmumuRatio()));
    } else {
        return myNPbase->muWHmumu(sqrt_s);
    }
}

muVHmumu::muVHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHmumu called with a class whose parent is not NPbase");
}

double muVHmumu::computeThValue()
{
    
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_WH  = 1.21539;
    double xsSM_ZH  = 0.795910;
    //double xsSM_WH = myNPbase->computeSigmaWH(sqrt_s);
    //double xsSM_ZH = myNPbase->computeSigmaZH(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHmumuRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ZH*(myNPbase->delta_muZH_1(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_1(sqrt_s)) )/(xsSM_ZH+xsSM_WH);     
        double muProd2 = ( xsSM_ZH*(myNPbase->delta_muZH_2(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_2(sqrt_s)) )/(xsSM_ZH+xsSM_WH);     
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //return myNPbase->muVHmumu(sqrt_s);
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuVHmumu = myNPbase->muVHmumu(sqrt_s);
        if(NPmuVHmumu==1.0){
            return ( xsSM_ZH*(myNPbase->muZH(sqrt_s)) + xsSM_WH*(myNPbase->muWH(sqrt_s)) )/(xsSM_ZH+xsSM_WH) * (myNPbase->BrHmumuRatio());
        } else {
            return NPmuVHmumu;
        }
        
    }
}

muttHmumu::muttHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHmumu called with a class whose parent is not NPbase");
}

double muttHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHmumuRatio()));
    } else {
        return myNPbase->muttHmumu(sqrt_s);
    }
}

muggHpttHptHpbbH_Hmumu::muggHpttHptHpbbH_Hmumu(const StandardModel& SM_i, const double sqrt_s_i)    //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpttHptHpbbH_Hmumu called with a class whose parent is not NPbase");
}

double muggHpttHptHpbbH_Hmumu::computeThValue()                                           //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this. Furthermore, the bbH is not included. In the SM this is
    //very suppressed (and it's probably also the case in the SMEFT) but 
    //in some NP models it may not be the case.Unfortunately, bbH is not
    //obtained in the SMEFT, we could just set delta_mubbH to zero in the
    //SMEFT and add here the general expression.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ggHbbH = 44.745;
    double xsSM_ttH    = 0.4998;
    double xsSM_tH     = 0.084769;
    //double xsSM_ggH = myNPbase->computeSigmaggH(sqrt_s);
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatH(sqrt_s);
    //double xsSM_bbH = myNPbase->computeSigmabbH(sqrt_s);

    
    
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ggHbbH*(myNPbase->delta_muggH_1(sqrt_s)) + xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ggHbbH+xsSM_ttH+xsSM_tH) ;
        double muProd2 = ( xsSM_ggHbbH*(myNPbase->delta_muggH_2(sqrt_s)) + xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ggHbbH+xsSM_ttH+xsSM_tH) ;
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end). Furthermore, the Hbb
        //is not added in the original formula, fine for the SM (probably also
        //for the SMEFT) but not for all NP models.
        double NPmuggHpttHptHpbbH_Hmumu = myNPbase->muggHpttHptHpbbH_Hmumu(sqrt_s);
        if(NPmuggHpttHptHpbbH_Hmumu==1.0){
            return ( (xsSM_ggHbbH*(myNPbase->muggH(sqrt_s))+xsSM_ttH*(myNPbase->muttH(sqrt_s))+xsSM_tH*(myNPbase->mutH(sqrt_s))) / (xsSM_ggHbbH+xsSM_ttH+xsSM_tH)  )*(myNPbase->BrHmumuRatio()) ;
        } else {
            return NPmuggHpttHptHpbbH_Hmumu;
        }  
    }
}

muVBFpVH_Hmumu::muVBFpVH_Hmumu(const StandardModel& SM_i, const double sqrt_s_i)  //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFpVH_Hmumu called with a class whose parent is not NPbase");
}

double muVBFpVH_Hmumu::computeThValue()                                          //AG:added
{
    
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this. 
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_VBF = 3.49948;
    double xsSM_WH  = 1.21539;
    double xsSM_ZH  = 0.795910;
    //double xsSM_VBF = myNPbase->computeSigmaVBF(sqrt_s);
    //double xsSM_WH = myNPbase->computeSigmaWH(sqrt_s);
    //double xsSM_ZH = myNPbase->computeSigmaZH(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_VBF*(myNPbase->delta_muVBF_1(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_1(sqrt_s)) + xsSM_ZH*(myNPbase->delta_muZH_1(sqrt_s)) )/(xsSM_VBF+xsSM_WH+xsSM_ZH);
        double muProd2 = ( xsSM_VBF*(myNPbase->delta_muVBF_2(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_2(sqrt_s)) + xsSM_ZH*(myNPbase->delta_muZH_2(sqrt_s)) )/(xsSM_VBF+xsSM_WH+xsSM_ZH);
        double dGammaR1 = myNPbase->deltaGammaHmumuRatio1();
        double dGammaR2 = myNPbase->deltaGammaHmumuRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end).
        double NPmuVBFpVH_Hmumu = myNPbase->muVBFpVH_Hmumu(sqrt_s);
        if(NPmuVBFpVH_Hmumu==1.0){
            return ( (xsSM_VBF*(myNPbase->muVBF(sqrt_s))+xsSM_WH*(myNPbase->muWH(sqrt_s))+xsSM_ZH*(myNPbase->muZH(sqrt_s))) / (xsSM_VBF+xsSM_WH+xsSM_ZH)  )*(myNPbase->BrHmumuRatio()) ;
        } else {
            return NPmuVBFpVH_Hmumu;
        } 
    }
}

muggHtautau::muggHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHtautau called with a class whose parent is not NPbase");
}

double muggHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHtautauRatio()));
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muggHtautau(sqrt_s);
    }
}

muVBFHtautau::muVBFHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHtautau called with a class whose parent is not NPbase");
}

double muVBFHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHtautauRatio()));
        double muProd1  = myNPbase->delta_muVBF_1(sqrt_s);
        double muProd2  = myNPbase->delta_muVBF_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muVBFHtautau(sqrt_s);
    }
}



muVBFpVHtautau::muVBFpVHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFpVHtautau called with a class whose parent is not NPbase");
}

double muVBFpVHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVBFpVH(sqrt_s)) + (myNPbase->BrHtautauRatio()));
    } else {
        return (myNPbase->muVBFpVH(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}


muZHtautau::muZHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHtautau called with a class whose parent is not NPbase");
}

double muZHtautau::computeThValue()                                             //AG:modified
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrHtautauRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHtautauRatio()));
        double muProd1  = myNPbase->delta_muZH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muZH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muZHtautau(sqrt_s);
    }
}

muWHtautau::muWHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHtautau called with a class whose parent is not NPbase");
}

double muWHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHtautauRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1  = myNPbase->delta_muWH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muWH_2(sqrt_s);    
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);
        
        /*std::cout << "muProd1="<<muProd1<<std::endl;
        std::cout << "muProd2="<<muProd2<<std::endl;
        std::cout << "dGammaR1="<<dGammaR1<<std::endl;
        std::cout << "dGammaR2="<<dGammaR2<<std::endl;
        std::cout << "dGammaRTot1="<<dGammaRTot1<<std::endl;
        std::cout << "dGammaRTot2="<<dGammaRTot2<<std::endl; 
        
        std::cout<<myNPbase->muWH(sqrt_s)<<std::endl;
        std::cout<<myNPbase->muWH(sqrt_s)<<std::endl;*/
        
        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu;   
        
    } else {
        return myNPbase->muWHtautau(sqrt_s);
    }
}

muVHtautau::muVHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHtautau called with a class whose parent is not NPbase");
}

double muVHtautau::computeThValue()
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_WH  = 1.21539;
    double xsSM_ZH  = 0.795910;
    //double xsSM_WH = myNPbase->computeSigmaWH(sqrt_s);
    //double xsSM_ZH = myNPbase->computeSigmaZH(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHtautauRatio()));    
        double muProd1 = ( xsSM_ZH*(myNPbase->delta_muZH_1(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_1(sqrt_s)) )/(xsSM_ZH+xsSM_WH);
        double muProd2 = ( xsSM_ZH*(myNPbase->delta_muZH_2(sqrt_s)) + xsSM_WH*(myNPbase->delta_muWH_2(sqrt_s)) )/(xsSM_ZH+xsSM_WH);     
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        //return myNPbase->muVHtautau(sqrt_s);
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuVHtautau = myNPbase->muVHtautau(sqrt_s);
        if(NPmuVHtautau==1.0){
            return ( xsSM_ZH*(myNPbase->muZH(sqrt_s)) + xsSM_WH*(myNPbase->muWH(sqrt_s)) )/(xsSM_ZH+xsSM_WH) * (myNPbase->BrHtautauRatio());
        } else {
            return NPmuVHtautau;
        }
    }
}

muttHtautau::muttHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHtautau called with a class whose parent is not NPbase");
}

double muttHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHtautauRatio()));
    } else {
        return myNPbase->muttHtautau(sqrt_s);
    }
}

muttHptH_Htautau::muttHptH_Htautau(const StandardModel& SM_i, const double sqrt_s_i)    //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHptH_Htautau called with a class whose parent is not NPbase");
}

double muttHptH_Htautau::computeThValue()                                           //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ttH = 0.499873;
    double xsSM_tH = 0.0821;
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatHq(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double muProd2 = ( xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ttH+xsSM_tH);  
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuttHptH_Htautau = myNPbase->muttHptH_Htautau(sqrt_s);
        if(NPmuttHptH_Htautau==1.0){
            return ( xsSM_ttH*(myNPbase->muttH(sqrt_s)) + xsSM_tH*(myNPbase->mutH(sqrt_s)) )/(xsSM_ttH+xsSM_tH) * (myNPbase->BrHtautauRatio()) ;
        } else {
            return NPmuttHptH_Htautau;
        }  
        
    }
}

muggHpbbH_Htautau::muggHpbbH_Htautau(const StandardModel& SM_i, const double sqrt_s_i)  //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpbbH_Htautau called with a class whose parent is not NPbase");
}

double muggHpbbH_Htautau::computeThValue()                                          //AG:added
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( 1.0 + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHtautauRatio()-1.));
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHtautauRatio1();
        double dGammaR2 = myNPbase->deltaGammaHtautauRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        //Also, the bbH is missing here, I'll leave it as it was for the 
        //moment (since bbH is really suppressed in the SM).
        double NPmuggHpbbH_Htautau = myNPbase->muggHpbbH_Htautau(sqrt_s);
        if (NPmuggHpbbH_Htautau == 1.0){
            return (myNPbase->muggHtautau(sqrt_s));
        } else{
            return NPmuggHpbbH_Htautau;
        } 
    }
}

muggHbb::muggHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHbb called with a class whose parent is not NPbase");
}

double muggHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        //return ( -1.0 + (myNPbase->muggH(sqrt_s)) + (myNPbase->BrHbbRatio()));
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1  = myNPbase->delta_muggH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muggH_2(sqrt_s);    
        double dGammaR1 = myNPbase->deltaGammaHbbRatio1();
        double dGammaR2 = myNPbase->deltaGammaHbbRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);     

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muggHbb(sqrt_s);
    }
}

muVBFHbb::muVBFHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHbb called with a class whose parent is not NPbase");
}

double muVBFHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHbbRatio()));
    } else {
        return myNPbase->muVBFHbb(sqrt_s);
    }
}

muZHbb::muZHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muZHbb called with a class whose parent is not NPbase");
}

double muZHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muZH(sqrt_s)) + (myNPbase->BrHbbRatio()));
        double muProd1  = myNPbase->delta_muZH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muZH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHbbRatio1();
        double dGammaR2 = myNPbase->deltaGammaHbbRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muZHbb(sqrt_s);
    }
}

muWHbb::muWHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muWHbb called with a class whose parent is not NPbase");
}

double muWHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        //return ( -1.0 + (myNPbase->muWH(sqrt_s)) + (myNPbase->BrHbbRatio()));
        double muProd1  = myNPbase->delta_muWH_1(sqrt_s);
        double muProd2  = myNPbase->delta_muWH_2(sqrt_s);
        double dGammaR1 = myNPbase->deltaGammaHbbRatio1();
        double dGammaR2 = myNPbase->deltaGammaHbbRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
        
    } else {
        return myNPbase->muWHbb(sqrt_s);
    }
}

muVHbb::muVHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHbb called with a class whose parent is not NPbase");
}

double muVHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHbbRatio()));
    } else {
        return myNPbase->muVHbb(sqrt_s);
    }
}

muttHbb::muttHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHbb called with a class whose parent is not NPbase");
}

double muttHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muttH(sqrt_s)) + (myNPbase->BrHbbRatio()));
    } else {
        return myNPbase->muttHbb(sqrt_s);
    }
}

muttHptH_Hbb::muttHptH_Hbb(const StandardModel& SM_i, const double sqrt_s_i)    //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHptH_Hbb called with a class whose parent is not NPbase");
}

double muttHptH_Hbb::computeThValue()                                           //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ttH = 0.499873;
    double xsSM_tH = 0.0821;
    //double xsSM_ttH = myNPbase->computeSigmattH(sqrt_s);
    //double xsSM_tH = myNPbase->computeSigmatHq(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ttH*(myNPbase->delta_muttH_1(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_1(sqrt_s)) )/(xsSM_ttH+xsSM_tH);
        double muProd2 = ( xsSM_ttH*(myNPbase->delta_muttH_2(sqrt_s)) + xsSM_tH*(myNPbase->delta_mutH_2(sqrt_s)) )/(xsSM_ttH+xsSM_tH);
        double dGammaR1 = myNPbase->deltaGammaHbbRatio1();
        double dGammaR2 = myNPbase->deltaGammaHbbRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
           
    } else {
        
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end)
        double NPmuttHptH_Hbb = myNPbase->muttHptH_Hbb(sqrt_s);
        if(NPmuttHptH_Hbb==1.0){
            return ( xsSM_ttH*(myNPbase->muttH(sqrt_s)) + xsSM_tH*(myNPbase->mutH(sqrt_s)) )/(xsSM_ttH+xsSM_tH) * (myNPbase->BrHbbRatio()) ;
        } else {
            return NPmuttHptH_Hbb;
        }   
    }
}

muggHpVBFpbbH_Hbb::muggHpVBFpbbH_Hbb(const StandardModel& SM_i, const double sqrt_s_i)      //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggHpVBFpbbH_Hbb called with a class whose parent is not NPbase");
}

double muggHpVBFpbbH_Hbb::computeThValue()                                      //AG:added
{
    //VM:Note that these values are valid for 13 TeV, they are not general
    //We should access the SM function that has all the values (for the 
    //different energies). The values are slightly different, we should 
    //check this. Furthermore, the bbH is not included. In the SM this is
    //very suppressed (and it's probably also the case in the SMEFT) but 
    //in some NP models it may not be the case.Unfortunately, bbH is not
    //obtained in the SMEFT, we could just set delta_mubbH to zero in the
    //SMEFT and add here the general expression.
    //Ref: https://www.hepdata.net/record/ins2104706 Figure2a (symmetrized)
    double xsSM_ggHbbH = 44.745;
    double xsSM_VBF    = 3.49948;
    //double xsSM_ggH = myNPbase->computeSigmaggH(sqrt_s);
    //double xsSM_VBF = myNPbase->computeSigmaVBF(sqrt_s);
    //double xsSM_bbH = myNPbase->computeSigmabbH(sqrt_s);
    if ((this->getModel()).isModelLinearized() || (this->getModel()).isModelNPquadratic()) {
        //AG: Most general expression including quadratic corrections. 
        //    If FlagQuadraticTerms=false, the quadratic pieces become 0 in NPSMEFTd6General
        double muProd1 = ( xsSM_ggHbbH*(myNPbase->delta_muggH_1(sqrt_s)) + xsSM_VBF*(myNPbase->delta_muVBF_1(sqrt_s)) )/(xsSM_ggHbbH+xsSM_VBF);  
        double muProd2 = ( xsSM_ggHbbH*(myNPbase->delta_muggH_2(sqrt_s)) + xsSM_VBF*(myNPbase->delta_muVBF_2(sqrt_s)) )/(xsSM_ggHbbH+xsSM_VBF);  
        double dGammaR1 = myNPbase->deltaGammaHbbRatio1();
        double dGammaR2 = myNPbase->deltaGammaHbbRatio2();
        
        double dGammaRTot1 = myNPbase->deltaGammaTotalRatio1();
        double dGammaRTot2 = myNPbase->deltaGammaTotalRatio2();
        double Br1 = dGammaR1-dGammaRTot1;
        double Br2 = dGammaR2 -dGammaRTot2 - dGammaR1*dGammaRTot1 + pow(dGammaRTot1,2.0);

        double mu;
        
        mu = 1.0 + (muProd1 + Br1) + (muProd2 + Br2 + muProd1*Br1);
        
        if (mu < 0) return std::numeric_limits<double>::quiet_NaN();
        
        return mu; 
           
    } else {
        
        
        //VM: Just in case someone wants to add directly the production*decay 
        //(which is the observable we fit at the end). Furthermore, the Hbb
        //is not added in the original formula, fine for the SM (probably also
        //for the SMEFT) but not for all NP models.
        double NPmuggHpVBFpbbH_Hbb = myNPbase->muggHpVBFpbbH_Hbb(sqrt_s);
        if(NPmuggHpVBFpbbH_Hbb==1.0){
            return ( xsSM_ggHbbH*(myNPbase->muggH(sqrt_s)) + xsSM_VBF*(myNPbase->muVBF(sqrt_s)) )/(xsSM_ggHbbH+xsSM_VBF) * (myNPbase->BrHbbRatio()) ;
        } else {
            return NPmuggHpVBFpbbH_Hbb;
        }       
    }
}



muVHcc::muVHcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHbb called with a class whose parent is not NPbase");
}

double muVHcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ( -1.0 + (myNPbase->muVH(sqrt_s)) + (myNPbase->BrHccRatio()));
    } else {
        return myNPbase->muVHcc(sqrt_s);
    }
}


muVBFBRinv::muVBFBRinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFBRinv called with a class whose parent is not NPbase");
}

double muVBFBRinv::computeThValue()
{

    return (myNPbase->muVBF(sqrt_s))*(myNPbase->Br_H_inv());

}

muVBFHinv::muVBFHinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVBFHinv called with a class whose parent is not NPbase");
}

double muVBFHinv::computeThValue()
{

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muVBF(sqrt_s)) + (myNPbase->BrHtoinvRatio()) - 1.0);
    } else {
        return (myNPbase->muVBF(sqrt_s))*(myNPbase->BrHtoinvRatio());
    }

}


muVHBRinv::muVHBRinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHBRinv called with a class whose parent is not NPbase");
}

double muVHBRinv::computeThValue()
{

    return (myNPbase->muVH(sqrt_s))*(myNPbase->Br_H_inv());

}

muVHinv::muVHinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muVHinv called with a class whose parent is not NPbase");
}

double muVHinv::computeThValue()
{

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muVH(sqrt_s)) + (myNPbase->BrHtoinvRatio()) - 1.0);
    } else {
        return (myNPbase->muVH(sqrt_s))*(myNPbase->BrHtoinvRatio());
    }

}


muppHmumu::muppHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muppHmumu called with a class whose parent is not NPbase");
}

double muppHmumu::computeThValue()
{
    return myNPbase->muppHmumu(sqrt_s);
}

muppHZga::muppHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muppHZga called with a class whose parent is not NPbase");
}

double muppHZga::computeThValue()
{
    return myNPbase->muppHZga(sqrt_s);
}

muggHH2ga2b::muggHH2ga2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muggH called with a class whose parent is not NPbase");
}

double muggHH2ga2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return (-2.0 + (myNPbase->muggHH(sqrt_s)) + (myNPbase->BrHgagaRatio()) + (myNPbase->BrHbbRatio()));
    } else {
        return (myNPbase->muggHH(sqrt_s))*(myNPbase->BrHgagaRatio())*(myNPbase->BrHbbRatio());
    }
}

muttHZbbboost::muttHZbbboost(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHZbbboost called with a class whose parent is not NPbase");
}

double muttHZbbboost::computeThValue()
{
    return (myNPbase->muttHZbbboost(sqrt_s));
}

muttHgagaZeeboost::muttHgagaZeeboost(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muttHgagaZeeboost called with a class whose parent is not NPbase");
}

double muttHgagaZeeboost::computeThValue()
{
    return (myNPbase->muttHgagaZeeboost(sqrt_s));
}

//AG:begin
ggHgaga::ggHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ggHgaga called with a class whose parent is not NPbase");
}
double ggHgaga::computeThValue()
{
    double SM_prediction = 0.0439;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHgagaRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muggH(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}


ggHZZ::ggHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ggHZZ called with a class whose parent is not NPbase");
}
double ggHZZ::computeThValue()
{
    double SM_prediction = 0.5197;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHZZRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muggH(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}


ggHWW::ggHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ggHWW called with a class whose parent is not NPbase");
}
double ggHWW::computeThValue()
{
    double SM_prediction = 4.1603;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHWWRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muggH(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}


ggHtautau::ggHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ggHtautau called with a class whose parent is not NPbase");
}
double ggHtautau::computeThValue()
{
    double SM_prediction = 1.2215;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muggH(sqrt_s)-1.) + (myNPbase->BrHtautauRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muggH(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}


VBFHgaga::VBFHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("VBFHgaga called with a class whose parent is not NPbase");
}
double VBFHgaga::computeThValue()
{
    double SM_prediction = 0.0037;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muVBF(sqrt_s)-1.) + (myNPbase->BrHgagaRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muVBF(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}


VBFHZZ::VBFHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("VBFHZZ called with a class whose parent is not NPbase");
}
double VBFHZZ::computeThValue()
{
    double SM_prediction = 0.0530;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muVBF(sqrt_s)-1.) + (myNPbase->BrHZZRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muVBF(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}


VBFHWW::VBFHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("VBFHWW called with a class whose parent is not NPbase");
}
double VBFHWW::computeThValue()
{
    double SM_prediction = 0.3494;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muVBF(sqrt_s)-1.) + (myNPbase->BrHWWRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muVBF(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}


VBFHtautau::VBFHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("VBFHtautau called with a class whose parent is not NPbase");
}
double VBFHtautau::computeThValue()
{
    double SM_prediction = 0.1011;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muVBF(sqrt_s)-1.) + (myNPbase->BrHtautauRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muVBF(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}


WHgaga::WHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("WHgaga called with a class whose parent is not NPbase");
}
double WHgaga::computeThValue()
{
    double SM_prediction = 0.0017;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muWH(sqrt_s)-1.) + (myNPbase->BrHgagaRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muWH(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}


WHWW::WHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("WHWW called with a class whose parent is not NPbase");
}
double WHWW::computeThValue()
{
    double SM_prediction = 0.1614;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muWH(sqrt_s)-1.) + (myNPbase->BrHWWRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muWH(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}


WHtautau::WHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("WWHtautau called with a class whose parent is not NPbase");
}
double WHtautau::computeThValue()
{
    double SM_prediction = 0.0462;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muWH(sqrt_s)-1.) + (myNPbase->BrHtautauRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muWH(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}


WHbb::WHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("WHbb called with a class whose parent is not NPbase");
}
double WHbb::computeThValue()
{
    double SM_prediction = 0.4090;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muWH(sqrt_s)-1.) + (myNPbase->BrHbbRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muWH(sqrt_s))*(myNPbase->BrHbbRatio());
    }    
}


ZHgaga::ZHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ZHgaga called with a class whose parent is not NPbase");
}
double ZHgaga::computeThValue()
{
    double SM_prediction = 0.0011;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muZH(sqrt_s)-1.) + (myNPbase->BrHgagaRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muZH(sqrt_s))*(myNPbase->BrHgagaRatio());
    } 
}


ZHWW::ZHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ZHWW called with a class whose parent is not NPbase");
}
double ZHWW::computeThValue()
{
    double SM_prediction = 0.0996;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muZH(sqrt_s)-1.) + (myNPbase->BrHWWRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muZH(sqrt_s))*(myNPbase->BrHWWRatio());
    } 
}


ZHtautau::ZHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ZHtautau called with a class whose parent is not NPbase");
}
double ZHtautau::computeThValue()
{
    double SM_prediction = 0.0304;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muZH(sqrt_s)-1.) + (myNPbase->BrHtautauRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muZH(sqrt_s))*(myNPbase->BrHtautauRatio());
    } 
}


ZHbb::ZHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ZHbb called with a class whose parent is not NPbase");
}
double ZHbb::computeThValue()
{
    double SM_prediction = 0.2410;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muZH(sqrt_s)-1.) + (myNPbase->BrHbbRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muZH(sqrt_s))*(myNPbase->BrHbbRatio());
    } 
}


ttHgaga::ttHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ttHgaga called with a class whose parent is not NPbase");
}
double ttHgaga::computeThValue()
{
    double SM_prediction = 0.0004;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muttH(sqrt_s)-1.) + (myNPbase->BrHgagaRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muttH(sqrt_s))*(myNPbase->BrHgagaRatio());
    } 
}


ttHWW::ttHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ttHWW called with a class whose parent is not NPbase");
}
double ttHWW::computeThValue()
{
    double SM_prediction = 0.0281;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muttH(sqrt_s)-1.) + (myNPbase->BrHWWRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muttH(sqrt_s))*(myNPbase->BrHWWRatio());
    } 
}


ttHtautau::ttHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ttHtautau called with a class whose parent is not NPbase");
}
double ttHtautau::computeThValue()
{
    double SM_prediction = 0.0106;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muttH(sqrt_s)-1.) + (myNPbase->BrHtautauRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muttH(sqrt_s))*(myNPbase->BrHtautauRatio());
    } 
}


ttHbb::ttHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("ttHbb called with a class whose parent is not NPbase");
}
double ttHbb::computeThValue()
{
    double SM_prediction = 0.0751;   //Ref:https://www.hepdata.net/record/ins1468068, Table1 (after symmetrizing)
    if ((this->getModel()).isModelLinearized()) {
        return SM_prediction*( 1. + (myNPbase->muttH(sqrt_s)-1.) + (myNPbase->BrHbbRatio()-1.) );
    } else {
        return SM_prediction*(myNPbase->muttH(sqrt_s))*(myNPbase->BrHbbRatio());
    } 
}

//AG:end

UpperLimit_ppHZgammaA::UpperLimit_ppHZgammaA(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("UpperLimit_ppHZgammaA called with a class whose parent is not NPbase");
}

double UpperLimit_ppHZgammaA::computeThValue()
{
    return myNPbase->UpperLimitZgammaA(sqrt_s);
}

UpperLimit_ppHZgammaA13::UpperLimit_ppHZgammaA13(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("UpperLimit_ppHZgammaA13 called with a class whose parent is not NPbase");
}

double UpperLimit_ppHZgammaA13::computeThValue()
{
    return myNPbase->UpperLimitZgammaA13(sqrt_s);
}

UpperLimit_ppHZgammaC13::UpperLimit_ppHZgammaC13(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("UpperLimit_ppHZgammaC13 called with a class whose parent is not NPbase");
}

double UpperLimit_ppHZgammaC13::computeThValue()
{
    return myNPbase->UpperLimitZgammaC13(sqrt_s);
}

UpperLimit_ppHZgammaC::UpperLimit_ppHZgammaC(const StandardModel& SM_i, const double sqrt_s_i) : ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("UpperLimit_ppHZgammaC called with a class whose parent is not NPbase");
}

double UpperLimit_ppHZgammaC::computeThValue()
{
    return myNPbase->UpperLimitZgammaC(sqrt_s);
}

cg_plus_ct::cg_plus_ct(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cg_plus_ct called with a class whose parent is not NPbase");
}

double cg_plus_ct::computeThValue()
{
    return myNPbase->cgplusct();
}

cga_plus_ct::cga_plus_ct(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cga_plus_ct called with a class whose parent is not NPbase");
}

double cga_plus_ct::computeThValue()
{
    return myNPbase->cgaplusct();
}

cg_minus_cga::cg_minus_cga(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cg_minus_cga called with a class whose parent is not NPbase");
}

double cg_minus_cga::computeThValue()
{
    return myNPbase->cgminuscga();
}

cV_plus_cb::cV_plus_cb(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cV_plus_cb called with a class whose parent is not NPbase");
}
\
    

double cV_plus_cb::computeThValue()
{
    return myNPbase->cVpluscb();
}

cV_plus_ctau::cV_plus_ctau(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cV_plus_ctau called with a class whose parent is not NPbase");
}

double cV_plus_ctau::computeThValue()
{
    return myNPbase->cVplusctau();
}

cb_minus_cc::cb_minus_cc(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cb_minus_cc called with a class whose parent is not NPbase");
}

double cb_minus_cc::computeThValue()
{
    return myNPbase->cbminuscc();
}

cb_minus_ctau::cb_minus_ctau(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cb_minus_ctau called with a class whose parent is not NPbase");
}

double cb_minus_ctau::computeThValue()
{
    return myNPbase->cbminusctau();
}

cc_minus_ctau::cc_minus_ctau(const StandardModel& SM_i) : ThObservable(SM_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("cc_minus_ctau called with a class whose parent is not NPbase");
}

double cc_minus_ctau::computeThValue()
{
    return myNPbase->ccminusctau();
}


//  Full signal strengths at e+ e- colliders
//  ----------------------------------------

mueeZHbb::mueeZHbb(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHbb called with a class whose parent is not NPbase");
}

double mueeZHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Hbb)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
    }
}

mueeZHcc::mueeZHcc(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHcc called with a class whose parent is not NPbase");
}

double mueeZHcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHccRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Hcc)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHccRatio());
    }
}

mueeZHss::mueeZHss(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHss called with a class whose parent is not NPbase");
}

double mueeZHss::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHssRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Hss)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHssRatio());
    }
}

mueeZHgg::mueeZHgg(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHgg called with a class whose parent is not NPbase");
}

double mueeZHgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHggRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Hgg)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHggRatio());
    }
}

mueeZHWW::mueeZHWW(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHWW called with a class whose parent is not NPbase");
}

double mueeZHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHWWRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1HWW)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHWWRatio());
    }
}

mueeZHtautau::mueeZHtautau(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHtautau called with a class whose parent is not NPbase");
}

double mueeZHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHtautauRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Htautau)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHtautauRatio());
    }
}

mueeZHZZ::mueeZHZZ(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHZZ called with a class whose parent is not NPbase");
}

double mueeZHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHZZRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1HZZ)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHZZRatio());
    }
}

mueeZHZga::mueeZHZga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHZga called with a class whose parent is not NPbase");
}

double mueeZHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHZgaRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1HZga)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHZgaRatio());
    }
}

mueeZHgaga::mueeZHgaga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHgaga called with a class whose parent is not NPbase");
}

double mueeZHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHgagaRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Hgaga)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHgagaRatio());
    }
}

mueeZHmumu::mueeZHmumu(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHmumu called with a class whose parent is not NPbase");
}

double mueeZHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHmumuRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1Hmumu)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHmumuRatio());
    }
}

mueeZHBRinv::mueeZHBRinv(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHBRinv called with a class whose parent is not NPbase");
}

double mueeZHBRinv::computeThValue()
{

    return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->Br_H_inv());

}

mueeZHinv::mueeZHinv(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZHinv called with a class whose parent is not NPbase");
}

double mueeZHinv::computeThValue()
{

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHtoinvRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZH(sqrt_s), myNPbase->C1HZZ)) - 1.0);
    } else {
        return (myNPbase->mueeZH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHtoinvRatio());
    }

}

mueeWBFbb::mueeWBFbb(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFbb called with a class whose parent is not NPbase");

}

double mueeWBFbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1Hbb)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
    }
}


mueeWBFcc::mueeWBFcc(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFcc called with a class whose parent is not NPbase");

}

double mueeWBFcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHccRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1Hcc)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHccRatio());
    }
}

mueeWBFgg::mueeWBFgg(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFgg called with a class whose parent is not NPbase");

}

double mueeWBFgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHggRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1Hgg)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHggRatio());
    }
}

mueeWBFWW::mueeWBFWW(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFWW called with a class whose parent is not NPbase");

}

double mueeWBFWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHWWRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1HWW)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHWWRatio());
    }
}

mueeWBFtautau::mueeWBFtautau(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFtautau called with a class whose parent is not NPbase");

}

double mueeWBFtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHtautauRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1Htautau)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHtautauRatio());
    }
}

mueeWBFZZ::mueeWBFZZ(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFZZ called with a class whose parent is not NPbase");

}

double mueeWBFZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHZZRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1HZZ)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHZZRatio());
    }
}

mueeWBFZga::mueeWBFZga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFZga called with a class whose parent is not NPbase");

}

double mueeWBFZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHZgaRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1HZga)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHZgaRatio());
    }
}

mueeWBFgaga::mueeWBFgaga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFgaga called with a class whose parent is not NPbase");

}

double mueeWBFgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHgagaRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1Hgaga)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHgagaRatio());
    }
}

mueeWBFmumu::mueeWBFmumu(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeWBFmumu called with a class whose parent is not NPbase");

}

double mueeWBFmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHmumuRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeWBF(sqrt_s), myNPbase->C1Hmumu)) - 1.0);
    } else {
        return (myNPbase->mueeWBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHmumuRatio());
    }
}

mueeHvvbb::mueeHvvbb(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvbb called with a class whose parent is not NPbase");

}

double mueeHvvbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Hbb)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
    }
}


mueeHvvcc::mueeHvvcc(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvcc called with a class whose parent is not NPbase");

}

double mueeHvvcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHccRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Hcc)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHccRatio());
    }
}


mueeHvvss::mueeHvvss(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvss called with a class whose parent is not NPbase");

}

double mueeHvvss::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHssRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Hss)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHssRatio());
    }
}


mueeHvvgg::mueeHvvgg(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvgg called with a class whose parent is not NPbase");

}

double mueeHvvgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHggRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Hgg)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHggRatio());
    }
}


mueeHvvWW::mueeHvvWW(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvWW called with a class whose parent is not NPbase");

}

double mueeHvvWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHWWRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1HWW)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHWWRatio());
    }
}


mueeHvvtautau::mueeHvvtautau(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvtautau called with a class whose parent is not NPbase");

}

double mueeHvvtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHtautauRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Htautau)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHtautauRatio());
    }
}


mueeHvvZZ::mueeHvvZZ(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvZZ called with a class whose parent is not NPbase");

}

double mueeHvvZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHZZRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1HZZ)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHZZRatio());
    }
}


mueeHvvZga::mueeHvvZga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvZga called with a class whose parent is not NPbase");

}

double mueeHvvZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHZgaRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1HZga)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHZgaRatio());
    }
}


mueeHvvgaga::mueeHvvgaga(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvgaga called with a class whose parent is not NPbase");

}

double mueeHvvgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHgagaRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Hgaga)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHgagaRatio());
    }
}


mueeHvvmumu::mueeHvvmumu(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeHvvmumu called with a class whose parent is not NPbase");

}

double mueeHvvmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHmumuRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeHvv(sqrt_s), myNPbase->C1Hmumu)) - 1.0);
    } else {
        return (myNPbase->mueeHvv(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHmumuRatio());
    }
}


mueeZBFbb::mueeZBFbb(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueeZBFbb called with a class whose parent is not NPbase");
}

double mueeZBFbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueeZBF(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eeZBF(sqrt_s), myNPbase->C1Hbb)) - 1.0);
    } else {
        return (myNPbase->mueeZBF(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
    }
}


mueettHbb::mueettHbb(const StandardModel& SM_i, const double sqrt_s_i, const double Pol_em_i, const double Pol_ep_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), Pol_em(Pol_em_i), Pol_ep(Pol_ep_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mueettHbb called with a class whose parent is not NPbase");
}

double mueettHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mueettH(sqrt_s, Pol_em, Pol_ep)) + (myNPbase->BrHbbRatio()) - 1.0); // + (myNPbase->delta2sBRH3(myNPbase->C1eettH(sqrt_s), myNPbase->C1Hbb)) - 1.0);
    } else {
        return (myNPbase->mueettH(sqrt_s, Pol_em, Pol_ep))*(myNPbase->BrHbbRatio());
    }
}



//  Production signal strengths at mu+ mu- colliders
//  ------------------------------------------------

mummZH::mummZH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZH called with a class whose parent is not NPbase");
}

double mummZH::computeThValue()
{
    return myNPbase->mummZH(sqrt_s);
}


mummHvv::mummHvv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvv called with a class whose parent is not NPbase");
}

double mummHvv::computeThValue()
{
    return myNPbase->mummHvv(sqrt_s);
}


mummHmm::mummHmm(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmm called with a class whose parent is not NPbase");
}

double mummHmm::computeThValue()
{
    return myNPbase->mummHmm(sqrt_s);
}


mummttH::mummttH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttH called with a class whose parent is not NPbase");
}

double mummttH::computeThValue()
{
    return myNPbase->mummttH(sqrt_s);
}


//  Full signal strengths at mu+ mu- colliders
//  -------------------------------------------


mummHbb::mummHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHbb called with a class whose parent is not NPbase");
}

double mummHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

mummHcc::mummHcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHcc called with a class whose parent is not NPbase");
}

double mummHcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

mummHgg::mummHgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHgg called with a class whose parent is not NPbase");
}

double mummHgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

mummHWW::mummHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHcc called with a class whose parent is not NPbase");
}

double mummHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}

mummHtautau::mummHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHtautau called with a class whose parent is not NPbase");
}

double mummHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

mummHZZ::mummHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHZZ called with a class whose parent is not NPbase");
}

double mummHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}

mummHZga::mummHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHZga called with a class whose parent is not NPbase");
}

double mummHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHZgaRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHZgaRatio());
    }
}

mummHgaga::mummHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHgaga called with a class whose parent is not NPbase");
}

double mummHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

mummHmumu::mummHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmumu called with a class whose parent is not NPbase");
}

double mummHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummH(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
    } else {
        return (myNPbase->mummH(sqrt_s))*(myNPbase->BrHmumuRatio());
    }
}



mummZHbb::mummZHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHbb called with a class whose parent is not NPbase");
}

double mummZHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

mummZHcc::mummZHcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHcc called with a class whose parent is not NPbase");
}

double mummZHcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

mummZHgg::mummZHgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHgg called with a class whose parent is not NPbase");
}

double mummZHgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

mummZHWW::mummZHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHcc called with a class whose parent is not NPbase");
}

double mummZHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}

mummZHtautau::mummZHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHtautau called with a class whose parent is not NPbase");
}

double mummZHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

mummZHZZ::mummZHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHZZ called with a class whose parent is not NPbase");
}

double mummZHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}

mummZHZga::mummZHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHZga called with a class whose parent is not NPbase");
}

double mummZHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHZgaRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHZgaRatio());
    }
}

mummZHgaga::mummZHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHgaga called with a class whose parent is not NPbase");
}

double mummZHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

mummZHmumu::mummZHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummZHmumu called with a class whose parent is not NPbase");
}

double mummZHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummZH(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
    } else {
        return (myNPbase->mummZH(sqrt_s))*(myNPbase->BrHmumuRatio());
    }
}



mummHvvbb::mummHvvbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvbb called with a class whose parent is not NPbase");
}

double mummHvvbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

mummHvvcc::mummHvvcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvcc called with a class whose parent is not NPbase");
}

double mummHvvcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

mummHvvgg::mummHvvgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvgg called with a class whose parent is not NPbase");
}

double mummHvvgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

mummHvvWW::mummHvvWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvcc called with a class whose parent is not NPbase");
}

double mummHvvWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}

mummHvvtautau::mummHvvtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvtautau called with a class whose parent is not NPbase");
}

double mummHvvtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

mummHvvZZ::mummHvvZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvZZ called with a class whose parent is not NPbase");
}

double mummHvvZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}

mummHvvZga::mummHvvZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvZga called with a class whose parent is not NPbase");
}

double mummHvvZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHZgaRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHZgaRatio());
    }
}

mummHvvgaga::mummHvvgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvgaga called with a class whose parent is not NPbase");
}

double mummHvvgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

mummHvvmumu::mummHvvmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHvvmumu called with a class whose parent is not NPbase");
}

double mummHvvmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHvv(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
    } else {
        return (myNPbase->mummHvv(sqrt_s))*(myNPbase->BrHmumuRatio());
    }
}



mummHmmbb::mummHmmbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmbb called with a class whose parent is not NPbase");
}

double mummHmmbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

mummHmmcc::mummHmmcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmcc called with a class whose parent is not NPbase");
}

double mummHmmcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

mummHmmgg::mummHmmgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmgg called with a class whose parent is not NPbase");
}

double mummHmmgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

mummHmmWW::mummHmmWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmcc called with a class whose parent is not NPbase");
}

double mummHmmWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}

mummHmmtautau::mummHmmtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmtautau called with a class whose parent is not NPbase");
}

double mummHmmtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

mummHmmZZ::mummHmmZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmZZ called with a class whose parent is not NPbase");
}

double mummHmmZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}

mummHmmZga::mummHmmZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmZga called with a class whose parent is not NPbase");
}

double mummHmmZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHZgaRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHZgaRatio());
    }
}

mummHmmgaga::mummHmmgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmgaga called with a class whose parent is not NPbase");
}

double mummHmmgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

mummHmmmumu::mummHmmmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmmmumu called with a class whose parent is not NPbase");
}

double mummHmmmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHmm(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
    } else {
        return (myNPbase->mummHmm(sqrt_s))*(myNPbase->BrHmumuRatio());
    }
}


mummttHbb::mummttHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHbb called with a class whose parent is not NPbase");
}

double mummttHbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

mummttHcc::mummttHcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHcc called with a class whose parent is not NPbase");
}

double mummttHcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

mummttHgg::mummttHgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHgg called with a class whose parent is not NPbase");
}

double mummttHgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

mummttHWW::mummttHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHcc called with a class whose parent is not NPbase");
}

double mummttHWW::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}

mummttHtautau::mummttHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHtautau called with a class whose parent is not NPbase");
}

double mummttHtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

mummttHZZ::mummttHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHZZ called with a class whose parent is not NPbase");
}

double mummttHZZ::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}

mummttHZga::mummttHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHZga called with a class whose parent is not NPbase");
}

double mummttHZga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHZgaRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHZgaRatio());
    }
}

mummttHgaga::mummttHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHgaga called with a class whose parent is not NPbase");
}

double mummttHgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

mummttHmumu::mummttHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummttHmumu called with a class whose parent is not NPbase");
}

double mummttHmumu::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummttH(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
    } else {
        return (myNPbase->mummttH(sqrt_s))*(myNPbase->BrHmumuRatio());
    }
}


// The same in the narrow width approximation

mummHbbNWA::mummHbbNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHbbNWA called with a class whose parent is not NPbase");
}

double mummHbbNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

mummHccNWA::mummHccNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHccNWA called with a class whose parent is not NPbase");
}

double mummHccNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

mummHggNWA::mummHggNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHggNWA called with a class whose parent is not NPbase");
}

double mummHggNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

mummHWWNWA::mummHWWNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHccNWA called with a class whose parent is not NPbase");
}

double mummHWWNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHWWRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHWWRatio());
    }
}

mummHtautauNWA::mummHtautauNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHtautauNWA called with a class whose parent is not NPbase");
}

double mummHtautauNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

mummHZZNWA::mummHZZNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHZZNWA called with a class whose parent is not NPbase");
}

double mummHZZNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHZZRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHZZRatio());
    }
}

mummHZgaNWA::mummHZgaNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHZgaNWA called with a class whose parent is not NPbase");
}

double mummHZgaNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHZgaRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHZgaRatio());
    }
}

mummHgagaNWA::mummHgagaNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHgagaNWA called with a class whose parent is not NPbase");
}

double mummHgagaNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

mummHmumuNWA::mummHmumuNWA(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("mummHmumuNWA called with a class whose parent is not NPbase");
}

double mummHmumuNWA::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->mummHNWA(sqrt_s)) + (myNPbase->BrHmumuRatio()) - 1.0);
    } else {
        return (myNPbase->mummHNWA(sqrt_s))*(myNPbase->BrHmumuRatio());
    }
}

//  Full signal strengths at ep colliders
//  -------------------------------------

muepWBFbb::muepWBFbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFbb called with a class whose parent is not NPbase");

}

double muepWBFbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

muepWBFcc::muepWBFcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFcc called with a class whose parent is not NPbase");

}

double muepWBFcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

muepWBFgg::muepWBFgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFgg called with a class whose parent is not NPbase");

}

double muepWBFgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

muepWBFWW2l2v::muepWBFWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFWW2l2v called with a class whose parent is not NPbase");

}

double muepWBFWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrH2l2vRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrH2l2vRatio());
    }
}

muepWBFZZ4l::muepWBFZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFZZ4l called with a class whose parent is not NPbase");

}

double muepWBFZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

muepWBFgaga::muepWBFgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFgaga called with a class whose parent is not NPbase");

}

double muepWBFgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

muepWBFtautau::muepWBFtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepWBFtautau called with a class whose parent is not NPbase");

}

double muepWBFtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepWBF(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->muepWBF(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}

muepZBFbb::muepZBFbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFbb called with a class whose parent is not NPbase");

}

double muepZBFbb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

muepZBFcc::muepZBFcc(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFcc called with a class whose parent is not NPbase");

}

double muepZBFcc::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrHccRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrHccRatio());
    }
}

muepZBFgg::muepZBFgg(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFgg called with a class whose parent is not NPbase");

}

double muepZBFgg::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrHggRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrHggRatio());
    }
}

muepZBFWW2l2v::muepZBFWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFWW2l2v called with a class whose parent is not NPbase");

}

double muepZBFWW2l2v::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrH2l2vRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrH2l2vRatio());
    }
}

muepZBFZZ4l::muepZBFZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFZZ4l called with a class whose parent is not NPbase");

}

double muepZBFZZ4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

muepZBFgaga::muepZBFgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFgaga called with a class whose parent is not NPbase");

}

double muepZBFgaga::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrHgagaRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrHgagaRatio());
    }
}

muepZBFtautau::muepZBFtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muepZBFtautau called with a class whose parent is not NPbase");

}

double muepZBFtautau::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->muepZBF(sqrt_s)) + (myNPbase->BrHtautauRatio()) - 1.0);
    } else {
        return (myNPbase->muepZBF(sqrt_s))*(myNPbase->BrHtautauRatio());
    }
}


// -----------------------------------------------------------------------------
// STXS bins
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Stage 0
// -----------------------------------------------------------------------------

STXS_0_qqH::STXS_0_qqH(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS_0_qqH called with a class whose parent is not NPbase");

}

double STXS_0_qqH::computeThValue()
{
    return myNPbase->STXS0_qqH(sqrt_s);
}


// -----------------------------------------------------------------------------
// Stage 1
// -----------------------------------------------------------------------------

STXSggH_VBFtopo_j3v_4l::STXSggH_VBFtopo_j3v_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH_VBFtopo_j3v_4l called with a class whose parent is not NPbase");

}

double STXSggH_VBFtopo_j3v_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH_VBFtopo_j3v(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH_VBFtopo_j3v(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSggH_VBFtopo_j3_4l::STXSggH_VBFtopo_j3_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH_VBFtopo_j3_4l called with a class whose parent is not NPbase");

}

double STXSggH_VBFtopo_j3_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH_VBFtopo_j3(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH_VBFtopo_j3(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSggH0j4l::STXSggH0j4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH0j4l called with a class whose parent is not NPbase");

}

double STXSggH0j4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH0j(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH0j(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSggH1j_pTH_0_60_4l::STXSggH1j_pTH_0_60_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH1j_pTH_0_60_4l called with a class whose parent is not NPbase");

}

double STXSggH1j_pTH_0_60_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH1j_pTH_0_60(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH1j_pTH_0_60(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSggH1j_pTH_60_120_4l::STXSggH1j_pTH_60_120_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH1j_pTH_60_120_4l called with a class whose parent is not NPbase");

}

double STXSggH1j_pTH_60_120_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH1j_pTH_60_120(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH1j_pTH_60_120(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSggH1j_pTH_120_200_4l::STXSggH1j_pTH_120_200_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH1j_pTH_120_200_4l called with a class whose parent is not NPbase");

}

double STXSggH1j_pTH_120_200_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH1j_pTH_120_200(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH1j_pTH_120_200(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSggH1j_pTH_200_4l::STXSggH1j_pTH_200_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH1j_pTH_200_4l called with a class whose parent is not NPbase");

}

double STXSggH1j_pTH_200_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH1j_pTH_200(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH1j_pTH_200(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSggH2j_pTH_0_200_4l::STXSggH2j_pTH_0_200_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH2j_pTH_0_200_4l called with a class whose parent is not NPbase");

}

double STXSggH2j_pTH_0_200_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH2j_pTH_0_200(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH2j_pTH_0_200(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSggH2j_pTH_0_60_4l::STXSggH2j_pTH_0_60_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH2j_pTH_0_60_4l called with a class whose parent is not NPbase");

}

double STXSggH2j_pTH_0_60_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH2j_pTH_0_60(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH2j_pTH_0_60(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSggH2j_pTH_60_120_4l::STXSggH2j_pTH_60_120_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH2j_pTH_60_120_4l called with a class whose parent is not NPbase");

}

double STXSggH2j_pTH_60_120_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH2j_pTH_60_120(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH2j_pTH_60_120(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSggH2j_pTH_120_200_4l::STXSggH2j_pTH_120_200_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH2j_pTH_120_200_4l called with a class whose parent is not NPbase");

}

double STXSggH2j_pTH_120_200_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH2j_pTH_120_200(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH2j_pTH_120_200(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSggH2j_pTH_200_4l::STXSggH2j_pTH_200_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSggH2j_pTH_200_4l called with a class whose parent is not NPbase");

}

double STXSggH2j_pTH_200_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ggH2j_pTH_200(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ggH2j_pTH_200(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHqq_VBFtopo_Rest_4l::STXSqqHqq_VBFtopo_Rest_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_VBFtopo_Rest_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_VBFtopo_Rest_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_VBFtopo_Rest(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_VBFtopo_Rest(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHqq_VBFtopo_j3v_4l::STXSqqHqq_VBFtopo_j3v_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_VBFtopo_j3v_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_VBFtopo_j3v_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_VBFtopo_j3v(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_VBFtopo_j3v(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}



STXSqqHqq_VBFtopo_j3_4l::STXSqqHqq_VBFtopo_j3_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_VBFtopo_j3_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_VBFtopo_j3_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_VBFtopo_j3(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_VBFtopo_j3(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}



STXSqqHqq_nonVHtopo_4l::STXSqqHqq_nonVHtopo_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_nonVHtopo_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_nonVHtopo_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_nonVHtopo(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_nonVHtopo(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}



STXSqqHqq_VHtopo_4l::STXSqqHqq_VHtopo_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_VHtopo_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_VHtopo_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_VHtopo(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_VHtopo(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHqq_Rest_4l::STXSqqHqq_Rest_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_Rest_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_Rest_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_Rest(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_Rest(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHqq_pTj_200_4l::STXSqqHqq_pTj_200_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHqq_pTj_200_4l called with a class whose parent is not NPbase");

}

double STXSqqHqq_pTj_200_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHqq_pTj_200(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHqq_pTj_200(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSqqHlv_pTV_0_250_4l::STXSqqHlv_pTV_0_250_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_0_250_4l called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_0_250_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_0_250(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_0_250(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHlv_pTV_0_150_4l::STXSqqHlv_pTV_0_150_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_0_150_4l called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_0_150_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_0_150(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_0_150(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}



STXSqqHlv_pTV_150_250_0j_4l::STXSqqHlv_pTV_150_250_0j_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_150_250_0j_4l called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_150_250_0j_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_150_250_0j(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_150_250_0j(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}



STXSqqHlv_pTV_150_250_1j_4l::STXSqqHlv_pTV_150_250_1j_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_150_250_1j_4l called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_150_250_1j_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_150_250_1j(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_150_250_1j(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHlv_pTV_250_4l::STXSqqHlv_pTV_250_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_250_4l called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_250_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_250(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_250(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSqqHll_pTV_0_150_4l::STXSqqHll_pTV_0_150_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_0_150_4l called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_0_150_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_0_150(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_0_150(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSqqHll_pTV_150_250_4l::STXSqqHll_pTV_150_250_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_150_250_4l called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_150_250_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_150_250(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_150_250(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHll_pTV_150_250_0j_4l::STXSqqHll_pTV_150_250_0j_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_150_250_0j_4l called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_150_250_0j_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_150_250_0j(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_150_250_0j(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHll_pTV_150_250_1j_4l::STXSqqHll_pTV_150_250_1j_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_150_250_1j_4l called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_150_250_1j_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_150_250_1j(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_150_250_1j(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHll_pTV_250_4l::STXSqqHll_pTV_250_4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_250_4l called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_250_4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_250(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_250(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}

STXSttHtH4l::STXSttHtH4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSttHtH4l called with a class whose parent is not NPbase");

}

double STXSttHtH4l::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ttHtH(sqrt_s)) + (myNPbase->BrH4lRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ttHtH(sqrt_s))*(myNPbase->BrH4lRatio());
    }
}


STXSqqHlv_pTV_0_250_bb::STXSqqHlv_pTV_0_250_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_0_250_bb called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_0_250_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_0_250(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_0_250(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSqqHlv_pTV_0_150_bb::STXSqqHlv_pTV_0_150_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_0_150_bb called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_0_150_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_0_150(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_0_150(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}



STXSqqHlv_pTV_150_250_0j_bb::STXSqqHlv_pTV_150_250_0j_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_150_250_0j_bb called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_150_250_0j_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_150_250_0j(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_150_250_0j(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}



STXSqqHlv_pTV_150_250_1j_bb::STXSqqHlv_pTV_150_250_1j_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_150_250_1j_bb called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_150_250_1j_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_150_250_1j(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_150_250_1j(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSqqHlv_pTV_250_bb::STXSqqHlv_pTV_250_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHlv_pTV_250_bb called with a class whose parent is not NPbase");

}

double STXSqqHlv_pTV_250_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHlv_pTV_250(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHlv_pTV_250(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

STXSqqHll_pTV_0_150_bb::STXSqqHll_pTV_0_150_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_0_150_bb called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_0_150_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_0_150(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_0_150(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}

STXSqqHll_pTV_150_250_bb::STXSqqHll_pTV_150_250_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_150_250_bb called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_150_250_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_150_250(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_150_250(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSqqHll_pTV_150_250_0j_bb::STXSqqHll_pTV_150_250_0j_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_150_250_0j_bb called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_150_250_0j_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_150_250_0j(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_150_250_0j(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSqqHll_pTV_150_250_1j_bb::STXSqqHll_pTV_150_250_1j_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_150_250_1j_bb called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_150_250_1j_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_150_250_1j(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_150_250_1j(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSqqHll_pTV_250_bb::STXSqqHll_pTV_250_bb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSqqHll_pTV_250_bb called with a class whose parent is not NPbase");

}

double STXSqqHll_pTV_250_bb::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_qqHll_pTV_250(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_qqHll_pTV_250(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSWHqqHqq_VBFtopo_j3v_2b::STXSWHqqHqq_VBFtopo_j3v_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSWHqqHqq_VBFtopo_j3v_2b called with a class whose parent is not NPbase");

}

double STXSWHqqHqq_VBFtopo_j3v_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_WHqqHqq_VBFtopo_j3v(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_WHqqHqq_VBFtopo_j3v(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSWHqqHqq_VBFtopo_j3_2b::STXSWHqqHqq_VBFtopo_j3_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSWHqqHqq_VBFtopo_j3_2b called with a class whose parent is not NPbase");

}

double STXSWHqqHqq_VBFtopo_j3_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_WHqqHqq_VBFtopo_j3(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_WHqqHqq_VBFtopo_j3(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSWHqqHqq_VH2j_2b::STXSWHqqHqq_VH2j_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSWHqqHqq_VH2j_2b called with a class whose parent is not NPbase");

}

double STXSWHqqHqq_VH2j_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_WHqqHqq_VH2j(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_WHqqHqq_VH2j(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSWHqqHqq_Rest_2b::STXSWHqqHqq_Rest_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSWHqqHqq_Rest_2b called with a class whose parent is not NPbase");

}

double STXSWHqqHqq_Rest_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_WHqqHqq_Rest(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_WHqqHqq_Rest(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSWHqqHqq_pTj1_200_2b::STXSWHqqHqq_pTj1_200_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSWHqqHqq_pTj1_200_2b called with a class whose parent is not NPbase");

}

double STXSWHqqHqq_pTj1_200_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_WHqqHqq_pTj1_200(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_WHqqHqq_pTj1_200(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSZHqqHqq_VBFtopo_j3v_2b::STXSZHqqHqq_VBFtopo_j3v_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSZHqqHqq_VBFtopo_j3v_2b called with a class whose parent is not NPbase");

}

double STXSZHqqHqq_VBFtopo_j3v_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ZHqqHqq_VBFtopo_j3v(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ZHqqHqq_VBFtopo_j3v(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSZHqqHqq_VBFtopo_j3_2b::STXSZHqqHqq_VBFtopo_j3_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSZHqqHqq_VBFtopo_j3_2b called with a class whose parent is not NPbase");

}

double STXSZHqqHqq_VBFtopo_j3_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ZHqqHqq_VBFtopo_j3(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ZHqqHqq_VBFtopo_j3(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSZHqqHqq_VH2j_2b::STXSZHqqHqq_VH2j_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSZHqqHqq_VH2j_2b called with a class whose parent is not NPbase");

}

double STXSZHqqHqq_VH2j_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ZHqqHqq_VH2j(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ZHqqHqq_VH2j(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSZHqqHqq_Rest_2b::STXSZHqqHqq_Rest_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSZHqqHqq_Rest_2b called with a class whose parent is not NPbase");

}

double STXSZHqqHqq_Rest_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ZHqqHqq_Rest(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ZHqqHqq_Rest(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}


STXSZHqqHqq_pTj1_200_2b::STXSZHqqHqq_pTj1_200_2b(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXSZHqqHqq_pTj1_200_2b called with a class whose parent is not NPbase");

}

double STXSZHqqHqq_pTj1_200_2b::computeThValue()
{
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS_ZHqqHqq_pTj1_200(sqrt_s)) + (myNPbase->BrHbbRatio()) - 1.0);
    } else {
        return (myNPbase->STXS_ZHqqHqq_pTj1_200(sqrt_s))*(myNPbase->BrHbbRatio());
    }
}



// -----------------------------------------------------------------------------
// Stage 1.2
// -----------------------------------------------------------------------------

STXS12_ggH_pTH200_300_Nj01::STXS12_ggH_pTH200_300_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH200_300_Nj01 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH200_300_Nj01::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_pTH200_300_Nj01 called with invalid argument for final state in fstate_i");
    }        
    
    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_pTH200_300_Nj01(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_pTH200_300_Nj01(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggH_pTH300_450_Nj01::STXS12_ggH_pTH300_450_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH300_450_Nj01 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH300_450_Nj01::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_pTH300_450_Nj01 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_pTH300_450_Nj01(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_pTH300_450_Nj01(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggH_pTH450_650_Nj01::STXS12_ggH_pTH450_650_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH450_650_Nj01 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH450_650_Nj01::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_pTH450_650_Nj01 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_pTH450_650_Nj01(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_pTH450_650_Nj01(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggH_pTH650_Inf_Nj01::STXS12_ggH_pTH650_Inf_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH650_Inf_Nj01 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH650_Inf_Nj01::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_pTH650_Inf_Nj01 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_pTH650_Inf_Nj01(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_pTH650_Inf_Nj01(sqrt_s))*(BrHXXRatio);
    }
}


// -----------------------------------------------------------------------------

STXS12_ggH_pTH10_Inf_Nj0::STXS12_ggH_pTH10_Inf_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH10_Inf_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH10_Inf_Nj0::computeThValue()
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    
    double muProd   = myNPbase->STXS12_ggH_pTH10_Inf_Nj0(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;
    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;

    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH10_Inf_Nj0 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
}





// -----------------------------------------------------------------------------

//VM:STXS2024;
STXS12_ggH_pTH0_10_Nj0::STXS12_ggH_pTH0_10_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH0_10_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH0_10_Nj0::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH0_10_Nj0(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (6.6369); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH0_10_Nj0 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}



// -----------------------------------------------------------------------------

//VM:STXS2024;
STXS12_ggH_pTH10_200_Nj0::STXS12_ggH_pTH10_200_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH10_200_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH10_200_Nj0::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH10_200_Nj0(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (20.642); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH10_200_Nj0 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}




// -----------------------------------------------------------------------------

//VM:STXS2024;
STXS12_ggH_pTH0_200_Nj0::STXS12_ggH_pTH0_200_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH0_200_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH0_200_Nj0::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = ( 6.6369*(myNPbase->STXS12_ggH_pTH0_10_Nj0(sqrt_s)) +
                       20.642*(myNPbase->STXS12_ggH_pTH10_200_Nj0(sqrt_s)) )
                        /(6.6369 + 20.642);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (6.6369 + 20.642); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH0_200_Nj0 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}




// -----------------------------------------------------------------------------

//VM:STXS2024;
STXS12_ggH_mjj0_350_pTH0_60_Nj1::STXS12_ggH_mjj0_350_pTH0_60_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH0_60_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj0_350_pTH0_60_Nj1::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.

    double muProd   = ((myNPbase->STXS12_ggH_mjj0_350_pTH0_60_Nj1(sqrt_s))  );
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (1.); //NOT in Ref: https://www.hepdata.net/record/ins2104706 Figure7. Missing SM prediction!!!
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH0_120_Nj1 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}







// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ggH_pTH0_60_Nj1::STXS12_ggH_pTH0_60_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH0_60_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH0_60_Nj1::computeThValue()                                 //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ggH_pTH0_60_Nj1(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 6.50045; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH0_60_Nj1 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ggH_pTH60_120_Nj1::STXS12_ggH_pTH60_120_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH60_120_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH60_120_Nj1::computeThValue()                               //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ggH_pTH60_120_Nj1(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 4.50294; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH60_120_Nj1 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ggH_pTH120_200_Nj1::STXS12_ggH_pTH120_200_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH120_200_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH120_200_Nj1::computeThValue()                              //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ggH_pTH120_200_Nj1(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.74712; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH120_200_Nj1 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}



// -----------------------------------------------------------------------------

//VM:STXS2025
STXS12_ggH_pTH60_200_Nj1::STXS12_ggH_pTH60_200_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH60_200_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH60_200_Nj1::computeThValue()                              //VM:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = (4.50294*(myNPbase->STXS12_ggH_pTH60_120_Nj1(sqrt_s))+
                      0.74712*(myNPbase->STXS12_ggH_pTH120_200_Nj1(sqrt_s))
                      )/(4.50294+0.74712);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (4.50294+0.74712); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH60_200_Nj1 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}


// -----------------------------------------------------------------------------

STXS12_ggH_mjj0_350_pTH0_60_Nj2::STXS12_ggH_mjj0_350_pTH0_60_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH0_60_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj0_350_pTH0_60_Nj2::computeThValue()
{
    
    
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    
    double muProd   = ((myNPbase->STXS12_ggH_mjj0_350_pTH0_60_Nj2(sqrt_s)));
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    
    
    
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (1.24); //Ref: CMS note HIG-21-018-PAS-V3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH0_60_Nj2 called with invalid argument for final state in fstate_i");
    } 

    
    
    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
}

// -----------------------------------------------------------------------------

STXS12_ggH_mjj0_350_pTH60_120_Nj2::STXS12_ggH_mjj0_350_pTH60_120_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH60_120_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj0_350_pTH60_120_Nj2::computeThValue()
{
    
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    
    
    double muProd   = ((myNPbase->STXS12_ggH_mjj0_350_pTH60_120_Nj2(sqrt_s)));
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    
    
    
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (2.00); //Ref: CMS note HIG-21-018-PAS-V3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH60_120_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
}





// -----------------------------------------------------------------------------

//VM:STXS2024;
STXS12_ggH_mjj0_350_pTH0_120_Nj2::STXS12_ggH_mjj0_350_pTH0_120_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH0_120_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj0_350_pTH0_120_Nj2::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    
    double muProd   = ((1.24)*(myNPbase->STXS12_ggH_mjj0_350_pTH0_60_Nj2(sqrt_s)) +
                       (2.00)*(myNPbase->STXS12_ggH_mjj0_350_pTH60_120_Nj2(sqrt_s)) )
                        /(1.24+2.00);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (2.9636); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
        //This weight is a bit smaller than the sum of the weights considered for the two bins that
        //generate this bin when combined 1.24+2.00=3.24, This is because the values have been obtained
        //from different papers and they give slightly different predictions. However, the uncertainty
        //given by the ATLAS paper that combines the bins is of 0.6 pb for the SM prediction, above the
        //difference that we find with the other SM prediction so it is fine.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH0_120_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}




// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ggH_mjj0_350_pTH120_200_Nj2::STXS12_ggH_mjj0_350_pTH120_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH120_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj0_350_pTH120_200_Nj2::computeThValue()                     //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ggH_mjj0_350_pTH120_200_Nj2(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.94325; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_mjj0_350_pTH120_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}




//VM:STXS2024;
STXS12_ggH_mjj350_Inf_pTH0_200_Nj2::STXS12_ggH_mjj350_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj350_Inf_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj350_Inf_pTH0_200_Nj2::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    //We should weight each xsection with the SM prediction. We need to check
    //this values, nevertheless, the difference between the two parametrisations
    //is extremely small (way beyond our current precision)
    double muProd   = ((myNPbase->STXS12_ggH_mjj350_700_pTH0_200_Nj2(sqrt_s)) +
                       (myNPbase->STXS12_ggH_mjj700_Inf_pTH0_200_Nj2(sqrt_s)) )
                        /(2.);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.87751); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_mjj350_Inf_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}


// -----------------------------------------------------------------------------

STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2::STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_mjj350_700_pTH0_200_ptHjj0_25_Nj2(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2::STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{

    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_mjj350_700_pTH0_200_ptHjj25_Inf_Nj2(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2::STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj0_25_Nj2(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2::STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggH_mjj700_Inf_pTH0_200_ptHjj25_Inf_Nj2(sqrt_s))*(BrHXXRatio);
    }
}





// -----------------------------------------------------------------------------

//VM:STXS2024;
STXS12_ggH_pTH0_200_Nj2::STXS12_ggH_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH0_200_Nj2::computeThValue()                                 
{
    
    
    
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    //Ref: Also CMS note HIG-21-018-PAS-v3
    //VM: In the reference above we have the SM predictions for less bins. Angelica should have
    //the ones for Madgraph and we could use those. In any case the parametrisation is very stable
    //we give the same weight to the bins for which we don't have the SM prediction (for the moment)!
    //The possible error that we could be introducing here is way below our precision!
    double muProd   = (
    1.24*(myNPbase->STXS12_ggH_mjj0_350_pTH0_60_Nj2(sqrt_s))+
    2.00*(myNPbase->STXS12_ggH_mjj0_350_pTH60_120_Nj2(sqrt_s))+
    0.94321*(myNPbase->STXS12_ggH_mjj0_350_pTH120_200_Nj2(sqrt_s))+
    0.87751/2*(myNPbase->STXS12_ggH_mjj350_700_pTH0_200_Nj2(sqrt_s)+myNPbase->STXS12_ggH_mjj700_Inf_pTH0_200_Nj2(sqrt_s))
                        )/(1.24+2.00+0.94321+0.87751);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (1.24+2.00+0.94321+0.87751); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}







// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH200_Inf::STXS12_ggH_pTH200_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH200_Inf called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH200_Inf::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = ( 0.45825*(myNPbase->STXS12_ggH_pTH200_300(sqrt_s)) +
                       0.10632*(myNPbase->STXS12_ggH_pTH300_450(sqrt_s)) +
                       0.017974*(myNPbase->STXS12_ggH_pTH450_Inf(sqrt_s)))
                        /(0.45825 + 0.10632 + 0.017974);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.45825 + 0.10632 + 0.017974); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH200_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}




// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH300_Inf::STXS12_ggH_pTH300_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH300_Inf called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH300_Inf::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = (0.10632*(myNPbase->STXS12_ggH_pTH300_450(sqrt_s)) +
                       0.017974*(myNPbase->STXS12_ggH_pTH450_Inf(sqrt_s)))
                        /(0.10632 + 0.017974);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.10632 + 0.017974); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH200_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}


// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH200_300::STXS12_ggH_pTH200_300(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH200_300 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH200_300::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH200_300(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.45825 ); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH200_300 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}



// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH300_450::STXS12_ggH_pTH300_450(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH300_450 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH300_450::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH300_450(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.10632; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH300_450 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}




// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH450_650::STXS12_ggH_pTH450_650(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH450_650 called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH450_650::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH450_650(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 1.; //Not available in Ref: https://www.hepdata.net/record/ins2104706 Figure7. 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH450_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}



// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH650_Inf::STXS12_ggH_pTH650_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH650_Inf called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH650_Inf::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH650_Inf(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 1.; //Not available in Ref: https://www.hepdata.net/record/ins2104706 Figure7. 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH450_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}



// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_ggH_pTH450_Inf::STXS12_ggH_pTH450_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggH_pTH450_Inf called with a class whose parent is not NPbase");

}

double STXS12_ggH_pTH450_Inf::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_ggH_pTH450_Inf(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.017974; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH450_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}


// -----------------------------------------------------------------------------

STXS12_ggHll_pTV0_75::STXS12_ggHll_pTV0_75(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggHll_pTV0_75 called with a class whose parent is not NPbase");

}

double STXS12_ggHll_pTV0_75::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggHll_pTV0_75 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggHll_pTV0_75(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggHll_pTV0_75(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggHll_pTV75_150::STXS12_ggHll_pTV75_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggHll_pTV75_150 called with a class whose parent is not NPbase");

}

double STXS12_ggHll_pTV75_150::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggHll_pTV75_150 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggHll_pTV75_150(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggHll_pTV75_150(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggHll_pTV150_250_Nj0::STXS12_ggHll_pTV150_250_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggHll_pTV150_250_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_ggHll_pTV150_250_Nj0::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggHll_pTV150_250_Nj0 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggHll_pTV150_250_Nj0(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggHll_pTV150_250_Nj0(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggHll_pTV150_250_Nj1::STXS12_ggHll_pTV150_250_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggHll_pTV150_250_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_ggHll_pTV150_250_Nj1::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggHll_pTV150_250_Nj1 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggHll_pTV150_250_Nj1(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggHll_pTV150_250_Nj1(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_ggHll_pTV250_Inf::STXS12_ggHll_pTV250_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ggHll_pTV250_Inf called with a class whose parent is not NPbase");

}

double STXS12_ggHll_pTV250_Inf::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_ggHll_pTV250_Inf called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_ggHll_pTV250_Inf(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_ggHll_pTV250_Inf(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_Nj0::STXS12_qqHqq_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_Nj0::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHqq_Nj called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHqq_Nj0(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHqq_Nj0(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_Nj1::STXS12_qqHqq_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_Nj1::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHqq_Nj1 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHqq_Nj1(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHqq_Nj1(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_mjj0_60_Nj2::STXS12_qqHqq_mjj0_60_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj0_60_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj0_60_Nj2::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj0_60_Nj2 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHqq_mjj0_60_Nj2(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHqq_mjj0_60_Nj2(sqrt_s))*(BrHXXRatio);
    }
}



// -----------------------------------------------------------------------------


//VM:STXS2024;
STXS12_qqHqq_VH_veto_Nj01::STXS12_qqHqq_VH_veto_Nj01(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_VH_veto_Nj01 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_VH_veto_Nj01::computeThValue()                                 
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    double muProd   = myNPbase->STXS12_qqHqq_VH_veto_Nj01(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (1.); //Not available in Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ggH_pTH200_300 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}






// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj60_120_Nj2::STXS12_qqHqq_mjj60_120_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj60_120_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj60_120_Nj2::computeThValue()
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj60_120_Nj2(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 1.0; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj60_120_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_mjj120_350_Nj2::STXS12_qqHqq_mjj120_350_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj120_350_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj120_350_Nj2::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj120_350_Nj2 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHqq_mjj120_350_Nj2(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHqq_mjj120_350_Nj2(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2::STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2::computeThValue()
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 1.0; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_pTH200_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 


    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
    
}

// -----------------------------------------------------------------------------







STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2::STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2::computeThValue()
{

    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.3948; //Ref: CMS-HIG-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj0_25_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
    
    
    
    
    
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2::STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2::computeThValue()
{
    
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.175166; //Ref: CMS-HIG-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
    
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2::STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2::computeThValue()
{
    
    
    
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.571956; //Ref: CMS-HIG-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj0_25_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
    
    
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2::STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2::computeThValue()
{

    
    
    
    
        //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.1647; //Ref: CMS-HIG-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
    
    
    
    
}






// -----------------------------------------------------------------------------



//VM:STXS2025
STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2::STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2::computeThValue()
{

    
    
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = (
                        0.175166*(myNPbase->STXS12_qqHqq_mjj350_700_pTH0_200_pTHjj25_Inf_Nj2(sqrt_s))+
                        0.1647*(myNPbase->STXS12_qqHqq_mjj700_Inf_pTH0_200_pTHjj25_Inf_Nj2(sqrt_s)) 
                        )/(0.175166+0.1647);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.1647; //Ref: CMS-HIG-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_pTH0_200_pTHjj25_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
    
}







// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj350_700_pTH0_200_Nj2::STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_700_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_700_pTH0_200_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.53537; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_700_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2::STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight =  0.25614; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2::STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight =  0.22408; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2::STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight =  0.21578; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2::STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.22408*(myNPbase->STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(sqrt_s))
                      + 0.21578*(myNPbase->STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(sqrt_s))
                      ) / (0.22408+0.21578);
    double muProd1  = muProd - 1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.22408+0.21578) ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj1000_Inf_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}



// -----------------------------------------------------------------------------

//VM:STXS2025
STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2::STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2::computeThValue()                   //VM:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.25614*(myNPbase->STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(sqrt_s))
                    + 0.22408*(myNPbase->STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(sqrt_s))
                    + 0.21578*(myNPbase->STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(sqrt_s))
                    )/(0.25614+0.22408+0.21578);
    double muProd1  = (muProd-1.0);
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.25614+0.22408+0.21578); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj700_Inf_pTH0_200_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}




// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2::STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.07372 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2::STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.07315 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

STXS12_qqHqq_mjj350_Inf_Nj2::STXS12_qqHqq_mjj350_Inf_Nj2(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)         //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_Nj2 called with a class whose parent is not NPbase");

}

double STXS12_qqHqq_mjj350_Inf_Nj2::computeThValue()                   //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
   
    double muProd   = (0.53537*myNPbase->STXS12_qqHqq_mjj350_700_pTH0_200_Nj2(sqrt_s)
                    + 0.25614*myNPbase->STXS12_qqHqq_mjj700_1000_pTH0_200_Nj2(sqrt_s)
                    + 0.22408*myNPbase->STXS12_qqHqq_mjj1000_1500_pTH0_200_Nj2(sqrt_s)
                    + 0.21578*myNPbase->STXS12_qqHqq_mjj1500_Inf_pTH0_200_Nj2(sqrt_s)
                    + 0.073727*myNPbase->STXS12_qqHqq_mjj350_1000_pTH200_Inf_Nj2(sqrt_s)
                    + 0.07315*myNPbase->STXS12_qqHqq_mjj1000_Inf_pTH200_Inf_Nj2(sqrt_s)) 
                    / (0.53537+0.25614+0.22408+0.21578+0.073727+0.07315) ;

    double muProd1  = muProd - 1.;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.53537+0.25614+0.22408+0.21578+0.073727+0.07315); //Ref: https://www.hepdata.net/record/ins2104706 Figure7.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    
    
    
    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
//        std::cout<<"\033[1;31m CHECK IF LINEAR  \033[0m" << std::endl;
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
//        std::cout<<"\033[1;31m CHECK IF QUADRATIC  \033[0m" << std::endl;
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
//        std::cout<<"\033[1;31m CHECK IF NEITHER LINEAR NOR QUADRATIC  \033[0m" << std::endl;
        return weight*(muProd)*(BrHXXRatio);
    } 
    
}

// -----------------------------------------------------------------------------

STXS12_qqHlv_pTV0_75::STXS12_qqHlv_pTV0_75(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV0_75 called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV0_75::computeThValue()                                   //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHlv_pTV0_75(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.21509; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

STXS12_qqHlv_pTV75_150::STXS12_qqHlv_pTV75_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV75_150 called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV75_150::computeThValue()                                 //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHlv_pTV75_150(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.13440; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV75_150 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHlv_pTV150_250_Nj0::STXS12_qqHlv_pTV150_250_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV150_250_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV150_250_Nj0::computeThValue()                            //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHlv_pTV150_250_Nj0(sqrt_s);
    double muProd1  = myNPbase->STXS12_qqHlv_pTV150_250_Nj0(sqrt_s) - 1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.04117; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHqq_mjj350_Inf_Nj2 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------






STXS12_qqHlv_pTV150_250_Nj1::STXS12_qqHlv_pTV150_250_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV150_250_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV150_250_Nj1::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV150_250_Nj1 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHlv_pTV150_250_Nj1(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHlv_pTV150_250_Nj1(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

STXS12_qqHlv_pTV250_Inf::STXS12_qqHlv_pTV250_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV250_Inf called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV250_Inf::computeThValue()
{
       
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (myNPbase->STXS12_qqHlv_pTV250_Inf(sqrt_s));
    double muProd1  = (muProd -1.);
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.01127+0.00339) ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV250_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
    
    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHlv_pTV0_150::STXS12_qqHlv_pTV0_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV0_150 called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV0_150::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.21509*(myNPbase->STXS12_qqHlv_pTV0_75(sqrt_s))
                      + 0.13440*(myNPbase->STXS12_qqHlv_pTV75_150(sqrt_s))/(0.21509+0.13440));
    double muProd1  = (0.21509*(myNPbase->STXS12_qqHlv_pTV0_75(sqrt_s)-1.0)
                      + 0.13440*(myNPbase->STXS12_qqHlv_pTV75_150(sqrt_s)-1.0)/(0.21509+0.13440));
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.21509+0.13440) ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV0_150 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHlv_pTV150_Inf::STXS12_qqHlv_pTV150_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV150_Inf called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV150_Inf::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.04117 * myNPbase->STXS12_qqHlv_pTV150_250_Nj0(sqrt_s) 
                    + 0.01004 * myNPbase->STXS12_qqHlv_pTV250_400(sqrt_s)
                    + 0.00214 * myNPbase->STXS12_qqHlv_pTV400_Inf(sqrt_s) )/(0.04117+0.01004+0.00214);
    double muProd1  = (0.04117 * (myNPbase->STXS12_qqHlv_pTV150_250_Nj0(sqrt_s) -1.0)
                    + 0.01004 * (myNPbase->STXS12_qqHlv_pTV250_400(sqrt_s)-1.0)
                    + 0.00214 * (myNPbase->STXS12_qqHlv_pTV400_Inf(sqrt_s)-1.0) )/(0.04117+0.01004+0.00214);
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.04117+0.01004+0.00214) ; 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV150_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHlv_pTV250_400::STXS12_qqHlv_pTV250_400(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV250_400 called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV250_400::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHlv_pTV250_400(sqrt_s) ;
    double muProd1  = myNPbase->STXS12_qqHlv_pTV250_400(sqrt_s)-1.0 ;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight =0.01004 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV250_400 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHlv_pTV400_Inf::STXS12_qqHlv_pTV400_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHlv_pTV400_Inf called with a class whose parent is not NPbase");

}

double STXS12_qqHlv_pTV400_Inf::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHlv_pTV400_Inf(sqrt_s) ;
    double muProd1  = muProd-1.0 ;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight =0.00214 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHlv_pTV400_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}


// -----------------------------------------------------------------------------

STXS12_qqHll_pTV0_75::STXS12_qqHll_pTV0_75(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV0_75 called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV0_75::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV0_75 called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHll_pTV0_75(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHll_pTV0_75(sqrt_s))*(BrHXXRatio);
    }
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHll_pTV75_150::STXS12_qqHll_pTV75_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV75_150 called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV75_150::computeThValue()
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHll_pTV75_150(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 1.0 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV75_150 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHll_pTV150_250_Nj0::STXS12_qqHll_pTV150_250_Nj0(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV150_250_Nj0 called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV150_250_Nj0::computeThValue()                            //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHll_pTV150_250_Nj0(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.0147 ; //Ref: CMS-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV150_250_Nj0 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

STXS12_qqHll_pTV150_250_Nj1::STXS12_qqHll_pTV150_250_Nj1(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV150_250_Nj1 called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV150_250_Nj1::computeThValue()
{
        
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHll_pTV150_250_Nj1(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.0168) ; //Ref: CMS-21-018-PAS-v3
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV150_250_Nj0 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 
    
    
    
    
    
}

// -----------------------------------------------------------------------------

STXS12_qqHll_pTV250_Inf::STXS12_qqHll_pTV250_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV250_Inf called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV250_Inf::computeThValue()
{
    double BrHXXRatio = 1.0;
    if (fstate == 1){
        BrHXXRatio = (myNPbase->STXS12_BrH4lRatio());
    } else if (fstate == 2){
        BrHXXRatio = (myNPbase->STXS12_BrHgagaRatio());
    } else if (fstate == 3){
        BrHXXRatio = (myNPbase->STXS12_BrHbbRatio());
    } else if (fstate == 4){
        BrHXXRatio = (myNPbase->STXS12_BrHevmuvRatio());
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV250_Inf called with invalid argument for final state in fstate_i");
    } 

    if ((this->getModel()).isModelLinearized()) {
        return ((myNPbase->STXS12_qqHll_pTV250_Inf(sqrt_s)) + (BrHXXRatio) - 1.0);
    } else {
        return (myNPbase->STXS12_qqHll_pTV250_Inf(sqrt_s))*(BrHXXRatio);
    }
}


// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHll_pTV0_150::STXS12_qqHll_pTV0_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV0_150 called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV0_150::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHll_pTV0_150(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.19845 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV0_150 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHll_pTV250_400::STXS12_qqHll_pTV250_400(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV250_400 called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV250_400::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHll_pTV250_400(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.00715 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV250_400 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHll_pTV400_Inf::STXS12_qqHll_pTV400_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV400_Inf called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV400_Inf::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_qqHll_pTV400_Inf(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.00126 ; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV400_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_qqHll_pTV150_Inf::STXS12_qqHll_pTV150_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll_pTV150_Inf called with a class whose parent is not NPbase");

}

double STXS12_qqHll_pTV150_Inf::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.0147*(myNPbase->STXS12_qqHll_pTV150_250_Nj0(sqrt_s))
                       +0.01683*(myNPbase->STXS12_qqHll_pTV150_250_Nj1(sqrt_s))
                       + 0.00715*(myNPbase->STXS12_qqHll_pTV250_400(sqrt_s))
                       + 0.00126*(myNPbase->STXS12_qqHll_pTV400_Inf(sqrt_s)) )/(0.0147+0.01683+0.00715+0.00126);
    double muProd1  = muProd-1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = (0.0147+0.01683+0.00715+0.00126);
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll_pTV150_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}








//VM:STXS2025
STXS12_qqHll::STXS12_qqHll(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_qqHll called with a class whose parent is not NPbase");

}

double STXS12_qqHll::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.19845*(myNPbase->STXS12_qqHll_pTV0_150(sqrt_s))
                       +0.0147*(myNPbase->STXS12_qqHll_pTV150_250_Nj0(sqrt_s))
                       +0.01683*(myNPbase->STXS12_qqHll_pTV150_250_Nj1(sqrt_s))
                       + 0.00715*(myNPbase->STXS12_qqHll_pTV250_400(sqrt_s))
                       + 0.00126*(myNPbase->STXS12_qqHll_pTV400_Inf(sqrt_s)) 
                        )/(0.19845+0.0147+0.01683+0.00715+0.00126);
    double muProd1  = muProd - 1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = (0.19845+0.0147+0.01683+0.00715+0.00126);
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_qqHll called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}




// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_VHlep::STXS12_VHlep(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_VHlep called with a class whose parent is not NPbase");

}

double STXS12_VHlep::computeThValue()                                //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = ( 0.71256*myNPbase->STXS12_qqHlv_pTV0_75(sqrt_s)
                    + 0.06739*myNPbase->STXS12_qqHlv_pTV75_150(sqrt_s)
                    + 0.03943*myNPbase->STXS12_qqHlv_pTV150_250_Nj0(sqrt_s)
                    + 0.01127*myNPbase->STXS12_qqHlv_pTV250_400(sqrt_s)
                    + 0.00339*myNPbase->STXS12_qqHlv_pTV400_Inf(sqrt_s)
                    + 0.07934*myNPbase->STXS12_qqHll_pTV0_150(sqrt_s) 
                    +0.0147*(myNPbase->STXS12_qqHll_pTV150_250_Nj0(sqrt_s))
                    +0.01683*(myNPbase->STXS12_qqHll_pTV150_250_Nj1(sqrt_s))
                    + 0.00746*myNPbase->STXS12_qqHll_pTV250_400(sqrt_s)
                    + 0.00043*myNPbase->STXS12_qqHll_pTV400_Inf(sqrt_s)) 
                    /(0.71256 + 0.06739 + 0.03943 + 0.01127 + 0.00339 + 0.07934 
                    + 0.0147 + 0.01683 + 0.00746 + 0.00043);

    double muProd1  = muProd - 1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = (0.71256 + 0.06739 + 0.03943 + 0.01127 + 0.00339 + 0.07934 
                    + 0.0147 + 0.01683 + 0.00746 + 0.00043); 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_VHlep called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------






//VM:STXS2025
STXS12_VHlep_pTV0_150::STXS12_VHlep_pTV0_150(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //VM:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_VHlep_pTV0_150 called with a class whose parent is not NPbase");

}

double STXS12_VHlep_pTV0_150::computeThValue()                                //VM:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = ( 0.71256*myNPbase->STXS12_qqHlv_pTV0_75(sqrt_s)
                    + 0.06739*myNPbase->STXS12_qqHlv_pTV75_150(sqrt_s)
                    + 0.07934*myNPbase->STXS12_qqHll_pTV0_150(sqrt_s) 
                    ) 
                    /(0.71256 + 0.06739 + 0.07934);
    double muProd1  = muProd - 1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = 1.0; 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_VHlep_pTV0_150 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------





//VM:STXS2024
STXS12_VHlep_pTV150_Inf::STXS12_VHlep_pTV150_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)       //VM:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_VHlep_pTV150_Inf called with a class whose parent is not NPbase");

}

double STXS12_VHlep_pTV150_Inf::computeThValue()                                //VM:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.03943*myNPbase->STXS12_qqHlv_pTV150_250_Nj0(sqrt_s)
                    + 0.01127*myNPbase->STXS12_qqHlv_pTV250_400(sqrt_s)
                    + 0.00339*myNPbase->STXS12_qqHlv_pTV400_Inf(sqrt_s) 
                    +0.0147*(myNPbase->STXS12_qqHll_pTV150_250_Nj0(sqrt_s))
                    +0.01683*(myNPbase->STXS12_qqHll_pTV150_250_Nj1(sqrt_s))
                    + 0.00746*myNPbase->STXS12_qqHll_pTV250_400(sqrt_s)
                    + 0.00043*myNPbase->STXS12_qqHll_pTV400_Inf(sqrt_s)) 
                    /(0.03943 + 0.01127 + 0.00339 + 0.0147 + 0.01683 + 0.00746 + 0.00043);
    double muProd1  = (muProd - 1.0);
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = (0.03943 + 0.01127 + 0.00339 + 0.0147 + 0.01683 + 0.00746 + 0.00043); 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_VHlep_pTV150_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------






//AG:STXS2024
STXS12_ttH_pTH0_60::STXS12_ttH_pTH0_60(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH0_60 called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH0_60::computeThValue()                                     //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ttH_pTH0_60(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.11821; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH0_60 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH60_120::STXS12_ttH_pTH60_120(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH60_120 called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH60_120::computeThValue()                                   //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ttH_pTH60_120(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.17813; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH60_120 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH0_120::STXS12_ttH_pTH0_120(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH0_120 called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH0_120::computeThValue()                                   //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = ( 0.11821*(myNPbase->STXS12_ttH_pTH0_60(sqrt_s) )
                       + 0.17813*(myNPbase->STXS12_ttH_pTH60_120(sqrt_s) ) )/(0.11821+0.17813);
    double muProd1  = ( 0.11821*(myNPbase->STXS12_ttH_pTH0_60(sqrt_s) -1.0)
                       + 0.17813*(myNPbase->STXS12_ttH_pTH60_120(sqrt_s)-1.0 ) )/(0.11821+0.17813);;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.17813; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH0_120 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH120_200::STXS12_ttH_pTH120_200(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH120_200 called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH120_200::computeThValue()                                  //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ttH_pTH120_200(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.12647; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH120_200 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH200_300::STXS12_ttH_pTH200_300(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH200_300 called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH200_300::computeThValue()                                  //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ttH_pTH200_300(sqrt_s) ;
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.05263; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH200_300 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

STXS12_ttH_pTH300_Inf::STXS12_ttH_pTH300_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH300_Inf called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH300_Inf::computeThValue()
{
    //VM: Modified to take the thinner bins
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //double muProd   = myNPbase->STXS12_ttH_pTH300_Inf(sqrt_s) ;
    double muProd   = (0.01903*(myNPbase->STXS12_ttH_pTH300_450(sqrt_s))
                       +0.00538*(myNPbase->STXS12_ttH_pTH450_Inf(sqrt_s))
                        )/(0.01903+0.00538);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.01903+0.00538); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH300_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH300_450::STXS12_ttH_pTH300_450(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)   //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH300_450 called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH300_450::computeThValue()                                  //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ttH_pTH300_450(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.01903; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH300_450 called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}
// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH450_Inf::STXS12_ttH_pTH450_Inf(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)   //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH450_Inf called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH450_Inf::computeThValue()                                  //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_ttH_pTH450_Inf(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.00538; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH450_Inf called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH_pTH300_Inf_add::STXS12_ttH_pTH300_Inf_add(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)   //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_pTH300_Inf_add called with a class whose parent is not NPbase");

}

double STXS12_ttH_pTH300_Inf_add::computeThValue()                                  //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Since adding bins, include partial weigths of SM_predictions 
    //(https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.)
    double muProd   = (0.01903*myNPbase->STXS12_ttH_pTH300_450(sqrt_s)
                    + 0.00538*myNPbase->STXS12_ttH_pTH450_Inf(sqrt_s)) / (0.01903+0.00538);
    double muProd1  = muProd - 1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = (0.01903+0.00538); //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH_pTH300_Inf_add called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_ttH::STXS12_ttH(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)   //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH called with a class whose parent is not NPbase");

}

double STXS12_ttH::computeThValue()                                  //AG:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    //Since adding bins, include partial weigths (SM_predictions)
    double muProd   = ( 0.09779*myNPbase->STXS12_ttH_pTH0_60(sqrt_s)
                    + 0.13884*myNPbase->STXS12_ttH_pTH60_120(sqrt_s)
                    + 0.05439*myNPbase->STXS12_ttH_pTH120_200(sqrt_s)
                    + 0.05303*myNPbase->STXS12_ttH_pTH200_300(sqrt_s)
                    + 0.00498*myNPbase->STXS12_ttH_pTH300_450(sqrt_s)
                    + 0.00060*myNPbase->STXS12_ttH_pTH450_Inf(sqrt_s) ) /(0.09779 + 0.13884 + 0.05439 + 0.05303 + 0.00498 + 0.00060 );
    double muProd1  = ( 0.09779*(myNPbase->STXS12_ttH_pTH0_60(sqrt_s)-1.0)
                    + 0.13884*(myNPbase->STXS12_ttH_pTH60_120(sqrt_s)-1.0)
                    + 0.05439*(myNPbase->STXS12_ttH_pTH120_200(sqrt_s)-1.0)
                    + 0.05303*(myNPbase->STXS12_ttH_pTH200_300(sqrt_s)-1.0)
                    + 0.00498*(myNPbase->STXS12_ttH_pTH300_450(sqrt_s)-1.0)
                    + 0.00060*(myNPbase->STXS12_ttH_pTH450_Inf(sqrt_s)-1.0) ) /(0.09779 + 0.13884 + 0.05439 + 0.05303 + 0.00498 + 0.00060 );
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = (0.09779 + 0.13884 + 0.05439 + 0.05303 + 0.00498 + 0.00060 ); 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    }  else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------

//AG:STXS2024
STXS12_tH::STXS12_tH(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_tH called with a class whose parent is not NPbase");

}

double STXS12_tH::computeThValue()                                              //AG:modified
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    double muProd   = myNPbase->STXS12_tH(sqrt_s);
    double muProd1  = muProd -1.0;
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        // Use for Cross-section [pb] with no Higgs-boson decay
        weight = 0.08207; //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    }  else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_tH called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

    
}



// -----------------------------------------------------------------------------



//VM:STXS2025
STXS12_ttH_tH::STXS12_ttH_tH(const StandardModel& SM_i, const double sqrt_s_i, unsigned int fstate_i)   //AG:added
: ThObservable(SM_i), sqrt_s(sqrt_s_i), fstate(fstate_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("STXS12_ttH_tH called with a class whose parent is not NPbase");

}

double STXS12_ttH_tH::computeThValue()                                  //VM:added
{
    //-- Production:
    double weight = 1.0; //If normalized to the SM
    //Ref: https://www.hepdata.net/record/ins2104706 Figure7. After symmetrizing.
    //Since adding bins, include partial weigths (SM_predictions)
    double muProd   = ( 0.09779*myNPbase->STXS12_ttH_pTH0_60(sqrt_s)
                    + 0.13884*myNPbase->STXS12_ttH_pTH60_120(sqrt_s)
                    + 0.05439*myNPbase->STXS12_ttH_pTH120_200(sqrt_s)
                    + 0.05303*myNPbase->STXS12_ttH_pTH200_300(sqrt_s)
                    + 0.00498*myNPbase->STXS12_ttH_pTH300_450(sqrt_s)
                    + 0.00060*myNPbase->STXS12_ttH_pTH450_Inf(sqrt_s) 
                    + 0.08207*myNPbase->STXS12_tH(sqrt_s) 
                    ) /(0.09779 + 0.13884 + 0.05439 + 0.05303 + 0.00498 + 0.00060 + 0.08207);
    double muProd1  = ( muProd - 1.0 );
    double muProd2  = 0.0;

    //-- Decay:
    double BrHXXRatio   = 1.0;
    double dBrHXXRatio1 = 0.0;
    double dBrHXXRatio2 = 0.0;
    if (fstate==0){
        weight = (0.09779 + 0.13884 + 0.05439 + 0.05303 + 0.00498 + 0.00060 + 0.08207); 
    } else if (fstate == 1){
        BrHXXRatio   = (myNPbase->STXS12_BrH4lRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 2){
        BrHXXRatio   = (myNPbase->BrHgagaRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 3){
        BrHXXRatio   = (myNPbase->BrHbbRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 4){
        BrHXXRatio   = (myNPbase->STXS12_BrHevmuvRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 5){
        BrHXXRatio   = (myNPbase->BrHtautauRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    }  else if (fstate == 6){
        BrHXXRatio   = (myNPbase->BrHWWRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else if (fstate == 7){
        BrHXXRatio   = (myNPbase->BrHZZRatio());
        dBrHXXRatio1 = BrHXXRatio - 1.0;
        dBrHXXRatio2 = 0.0;
    } else {
        throw std::runtime_error("STXS12_ttH called with invalid argument for final state in fstate_i");
    } 

    //-- Production x Decay:    
    if ((this->getModel()).isModelLinearized()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1));
    } else if((this->getModel()).isModelNPquadratic()){
        return weight*( 1.0 + (muProd1 + dBrHXXRatio1) + (muProd2 + dBrHXXRatio2 + muProd1*dBrHXXRatio1) );    
    } else {
        return weight*(muProd)*(BrHXXRatio);
    } 

}

// -----------------------------------------------------------------------------


//-----------------------------------------------------------------------------------------
//-- Special Hadron collider signal strengths with separate full TH unc U(prod x decay) ---
//-----------------------------------------------------------------------------------------


muTHUggHgaga::muTHUggHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHgaga called with a class whose parent is not NPbase");
}

double muTHUggHgaga::computeThValue()
{
        return myNPbase->muTHUggHgaga(sqrt_s);
}


muTHUVBFHgaga::muTHUVBFHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHgaga called with a class whose parent is not NPbase");
}

double muTHUVBFHgaga::computeThValue()
{
        return myNPbase->muTHUVBFHgaga(sqrt_s);
}

muTHUZHgaga::muTHUZHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHgaga called with a class whose parent is not NPbase");
}

double muTHUZHgaga::computeThValue()
{
        return myNPbase->muTHUZHgaga(sqrt_s);
}

muTHUWHgaga::muTHUWHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHgaga called with a class whose parent is not NPbase");
}

double muTHUWHgaga::computeThValue()
{
        return myNPbase->muTHUWHgaga(sqrt_s);
}

muTHUVHgaga::muTHUVHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHgaga called with a class whose parent is not NPbase");
}

double muTHUVHgaga::computeThValue()
{
        return myNPbase->muTHUVHgaga(sqrt_s);
}

muTHUttHgaga::muTHUttHgaga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHgaga called with a class whose parent is not NPbase");
}

double muTHUttHgaga::computeThValue()
{
        return myNPbase->muTHUttHgaga(sqrt_s);
}

muTHUggHZga::muTHUggHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHZga called with a class whose parent is not NPbase");
}

double muTHUggHZga::computeThValue()
{
        return myNPbase->muTHUggHZga(sqrt_s);
}

muTHUggHZgamumu::muTHUggHZgamumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHZgamumu called with a class whose parent is not NPbase");
}

double muTHUggHZgamumu::computeThValue()
{
        return (myNPbase->muTHUggHZgamumu(sqrt_s));
}

muTHUVBFHZga::muTHUVBFHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHZga called with a class whose parent is not NPbase");
}

double muTHUVBFHZga::computeThValue()
{

        return myNPbase->muTHUVBFHZga(sqrt_s);
}

muTHUZHZga::muTHUZHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHZga called with a class whose parent is not NPbase");
}

double muTHUZHZga::computeThValue()
{
        return myNPbase->muTHUZHZga(sqrt_s);
}

muTHUWHZga::muTHUWHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHZga called with a class whose parent is not NPbase");
}

double muTHUWHZga::computeThValue()
{
        return myNPbase->muTHUWHZga(sqrt_s);
}

muTHUVHZga::muTHUVHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHZga called with a class whose parent is not NPbase");
}

double muTHUVHZga::computeThValue()
{
        return myNPbase->muTHUVHZga(sqrt_s);
}

muTHUttHZga::muTHUttHZga(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHZga called with a class whose parent is not NPbase");
}

double muTHUttHZga::computeThValue()
{
        return myNPbase->muTHUttHZga(sqrt_s);
}

muTHUggHZZ::muTHUggHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHZZ called with a class whose parent is not NPbase");
}

double muTHUggHZZ::computeThValue()
{
        return myNPbase->muTHUggHZZ(sqrt_s);
}

muTHUVBFHZZ::muTHUVBFHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHZZ called with a class whose parent is not NPbase");
}

double muTHUVBFHZZ::computeThValue()
{
        return myNPbase->muTHUVBFHZZ(sqrt_s);
}

muTHUZHZZ::muTHUZHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHZZ called with a class whose parent is not NPbase");
}

double muTHUZHZZ::computeThValue()
{
        return myNPbase->muTHUZHZZ(sqrt_s);
}

muTHUWHZZ::muTHUWHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHZZ called with a class whose parent is not NPbase");
}

double muTHUWHZZ::computeThValue()
{
        return myNPbase->muTHUWHZZ(sqrt_s);
}

muTHUVHZZ::muTHUVHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHZZ called with a class whose parent is not NPbase");
}

double muTHUVHZZ::computeThValue()
{
        return myNPbase->muTHUVHZZ(sqrt_s);
}

muTHUttHZZ::muTHUttHZZ(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHZZ called with a class whose parent is not NPbase");
}

double muTHUttHZZ::computeThValue()
{
        return myNPbase->muTHUttHZZ(sqrt_s);
}

muTHUggHZZ4l::muTHUggHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHZZ4l called with a class whose parent is not NPbase");
}

double muTHUggHZZ4l::computeThValue()
{
        return myNPbase->muTHUggHZZ4l(sqrt_s);
}

muTHUggHZZ4mu::muTHUggHZZ4mu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHZZ4mu called with a class whose parent is not NPbase");
}

double muTHUggHZZ4mu::computeThValue()
{
        return (myNPbase->muTHUggHZZ4mu(sqrt_s));
}

muTHUVBFHZZ4l::muTHUVBFHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHZZ4l called with a class whose parent is not NPbase");
}

double muTHUVBFHZZ4l::computeThValue()
{
        return myNPbase->muTHUVBFHZZ4l(sqrt_s);
}

muTHUZHZZ4l::muTHUZHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHZZ4l called with a class whose parent is not NPbase");
}

double muTHUZHZZ4l::computeThValue()
{
        return myNPbase->muTHUZHZZ4l(sqrt_s);
}

muTHUWHZZ4l::muTHUWHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHZZ4l called with a class whose parent is not NPbase");
}

double muTHUWHZZ4l::computeThValue()
{
        return myNPbase->muTHUWHZZ4l(sqrt_s);
}

muTHUVHZZ4l::muTHUVHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHZZ4l called with a class whose parent is not NPbase");
}

double muTHUVHZZ4l::computeThValue()
{
        return myNPbase->muTHUVHZZ4l(sqrt_s);
}

muTHUttHZZ4l::muTHUttHZZ4l(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHZZ4l called with a class whose parent is not NPbase");
}

double muTHUttHZZ4l::computeThValue()
{
        return myNPbase->muTHUttHZZ4l(sqrt_s);
}

muTHUggHWW::muTHUggHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHWW called with a class whose parent is not NPbase");
}

double muTHUggHWW::computeThValue()
{
        return myNPbase->muTHUggHWW(sqrt_s);
}

muTHUVBFHWW::muTHUVBFHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHWW called with a class whose parent is not NPbase");
}

double muTHUVBFHWW::computeThValue()
{
        return myNPbase->muTHUVBFHWW(sqrt_s);
}

muTHUZHWW::muTHUZHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHWW called with a class whose parent is not NPbase");
}

double muTHUZHWW::computeThValue()
{
        return myNPbase->muTHUZHWW(sqrt_s);
}

muTHUWHWW::muTHUWHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHWW called with a class whose parent is not NPbase");
}

double muTHUWHWW::computeThValue()
{
        return myNPbase->muTHUWHWW(sqrt_s);
}

muTHUVHWW::muTHUVHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHWW called with a class whose parent is not NPbase");
}

double muTHUVHWW::computeThValue()
{
        return myNPbase->muTHUVHWW(sqrt_s);
}

muTHUttHWW::muTHUttHWW(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHWW called with a class whose parent is not NPbase");
}

double muTHUttHWW::computeThValue()
{
        return myNPbase->muTHUttHWW(sqrt_s);
}

muTHUggHWW2l2v::muTHUggHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHWW2l2v called with a class whose parent is not NPbase");
}

double muTHUggHWW2l2v::computeThValue()
{
        return myNPbase->muTHUggHWW2l2v(sqrt_s);
}

muTHUVBFHWW2l2v::muTHUVBFHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHWW2l2v called with a class whose parent is not NPbase");
}

double muTHUVBFHWW2l2v::computeThValue()
{
        return myNPbase->muTHUVBFHWW2l2v(sqrt_s);
}

muTHUZHWW2l2v::muTHUZHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHWW2l2v called with a class whose parent is not NPbase");
}

double muTHUZHWW2l2v::computeThValue()
{
        return myNPbase->muTHUZHWW2l2v(sqrt_s);
}

muTHUWHWW2l2v::muTHUWHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHWW2l2v called with a class whose parent is not NPbase");
}

double muTHUWHWW2l2v::computeThValue()
{
        return myNPbase->muTHUWHWW2l2v(sqrt_s);
}

muTHUVHWW2l2v::muTHUVHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHWW2l2v called with a class whose parent is not NPbase");
}

double muTHUVHWW2l2v::computeThValue()
{
        return myNPbase->muTHUVHWW2l2v(sqrt_s);
}

muTHUttHWW2l2v::muTHUttHWW2l2v(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHWW2l2v called with a class whose parent is not NPbase");
}

double muTHUttHWW2l2v::computeThValue()
{
        return myNPbase->muTHUttHWW2l2v(sqrt_s);
}

muTHUggHmumu::muTHUggHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHmumu called with a class whose parent is not NPbase");
}

double muTHUggHmumu::computeThValue()
{
        return myNPbase->muTHUggHmumu(sqrt_s);
}

muTHUVBFHmumu::muTHUVBFHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHmumu called with a class whose parent is not NPbase");
}

double muTHUVBFHmumu::computeThValue()
{
        return myNPbase->muTHUVBFHmumu(sqrt_s);
}

muTHUZHmumu::muTHUZHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHmumu called with a class whose parent is not NPbase");
}

double muTHUZHmumu::computeThValue()
{
        return myNPbase->muTHUZHmumu(sqrt_s);
}

muTHUWHmumu::muTHUWHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHmumu called with a class whose parent is not NPbase");
}

double muTHUWHmumu::computeThValue()
{
        return myNPbase->muTHUWHmumu(sqrt_s);
}

muTHUVHmumu::muTHUVHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHmumu called with a class whose parent is not NPbase");
}

double muTHUVHmumu::computeThValue()
{
        return myNPbase->muTHUVHmumu(sqrt_s);
}

muTHUttHmumu::muTHUttHmumu(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHmumu called with a class whose parent is not NPbase");
}

double muTHUttHmumu::computeThValue()
{
        return myNPbase->muTHUttHmumu(sqrt_s);
}

muTHUggHtautau::muTHUggHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHtautau called with a class whose parent is not NPbase");
}

double muTHUggHtautau::computeThValue()
{
        return myNPbase->muTHUggHtautau(sqrt_s);
}

muTHUVBFHtautau::muTHUVBFHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHtautau called with a class whose parent is not NPbase");
}

double muTHUVBFHtautau::computeThValue()
{
        return myNPbase->muTHUVBFHtautau(sqrt_s);
}

muTHUZHtautau::muTHUZHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHtautau called with a class whose parent is not NPbase");
}

double muTHUZHtautau::computeThValue()
{
        return myNPbase->muTHUZHtautau(sqrt_s);
}

muTHUWHtautau::muTHUWHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHtautau called with a class whose parent is not NPbase");
}

double muTHUWHtautau::computeThValue()
{
        return myNPbase->muTHUWHtautau(sqrt_s);
}

muTHUVHtautau::muTHUVHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHtautau called with a class whose parent is not NPbase");
}

double muTHUVHtautau::computeThValue()
{
        return myNPbase->muTHUVHtautau(sqrt_s);
}

muTHUttHtautau::muTHUttHtautau(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHtautau called with a class whose parent is not NPbase");
}

double muTHUttHtautau::computeThValue()
{
        return myNPbase->muTHUttHtautau(sqrt_s);
}

muTHUggHbb::muTHUggHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUggHbb called with a class whose parent is not NPbase");
}

double muTHUggHbb::computeThValue()
{
        return myNPbase->muTHUggHbb(sqrt_s);
}

muTHUVBFHbb::muTHUVBFHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHbb called with a class whose parent is not NPbase");
}

double muTHUVBFHbb::computeThValue()
{
        return myNPbase->muTHUVBFHbb(sqrt_s);
}

muTHUZHbb::muTHUZHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUZHbb called with a class whose parent is not NPbase");
}

double muTHUZHbb::computeThValue()
{
        return myNPbase->muTHUZHbb(sqrt_s);
}

muTHUWHbb::muTHUWHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUWHbb called with a class whose parent is not NPbase");
}

double muTHUWHbb::computeThValue()
{
        return myNPbase->muTHUWHbb(sqrt_s);
}

muTHUVHbb::muTHUVHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHbb called with a class whose parent is not NPbase");
}

double muTHUVHbb::computeThValue()
{
        return myNPbase->muTHUVHbb(sqrt_s);
}

muTHUttHbb::muTHUttHbb(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUttHbb called with a class whose parent is not NPbase");
}

double muTHUttHbb::computeThValue()
{
        return myNPbase->muTHUttHbb(sqrt_s);
}


muTHUVBFBRinv::muTHUVBFBRinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFBRinv called with a class whose parent is not NPbase");
}

double muTHUVBFBRinv::computeThValue()
{

    return (myNPbase->muTHUVBFBRinv(sqrt_s));

}

muTHUVBFHinv::muTHUVBFHinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVBFHinv called with a class whose parent is not NPbase");
}

double muTHUVBFHinv::computeThValue()
{

        return (myNPbase->muTHUVBFHinv(sqrt_s));

}


muTHUVHBRinv::muTHUVHBRinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHBRinv called with a class whose parent is not NPbase");
}

double muTHUVHBRinv::computeThValue()
{

    return (myNPbase->muTHUVHBRinv(sqrt_s));

}

muTHUVHinv::muTHUVHinv(const StandardModel& SM_i, const double sqrt_s_i)
: ThObservable(SM_i), sqrt_s(sqrt_s_i)
{
    if ((myNPbase = dynamic_cast<const NPbase*> (&SM)) == NULL)
        throw std::runtime_error("muTHUVHinv called with a class whose parent is not NPbase");
}

double muTHUVHinv::computeThValue()
{
        return (myNPbase->muTHUVHinv(sqrt_s));
}
