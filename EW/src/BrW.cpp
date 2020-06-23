/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BrW.h"
#include "StandardModel.h"

double BrWlepton::computeThValue()
{
    double BRe = SM.BrW(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON));
    double BRmu = SM.BrW(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU));
    double BRtau = SM.BrW(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU));
    
    return (1./3.)*(BRe + BRmu + BRtau);
}

double BrWelectron::computeThValue()
{
    return SM.BrW(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON));
}

double BrWmuon::computeThValue()
{
    return SM.BrW(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU));
}

double BrWtau::computeThValue()
{
    return SM.BrW(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU));
}

double BrWhadrons::computeThValue()
{
    double Brhad = 0.;
    
//  Current SM formula asummes diagonal CKM
    Brhad += SM.BrW(SM.getQuarks(StandardModel::UP), SM.getQuarks(StandardModel::DOWN));
    Brhad += SM.BrW(SM.getQuarks(StandardModel::CHARM), SM.getQuarks(StandardModel::STRANGE));

    return Brhad;
}
