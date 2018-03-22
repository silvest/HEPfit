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
    double GammW = SM.GammaW();
    double Gamme = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON));
    double Gammmu = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU));
    double Gammtau = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU));
    
    return (1./3.)*(Gamme + Gammmu + Gammtau)/GammW;
}

double BrWelectron::computeThValue()
{
    double GammW = SM.GammaW();
    double Gamme = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON));

    return Gamme/GammW;
}

double BrWmuon::computeThValue()
{
    double GammW = SM.GammaW();
    double Gammmu = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU));

    return Gammmu/GammW;
}

double BrWtau::computeThValue()
{
    double GammW = SM.GammaW();
    double Gammtau = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU));

    return Gammtau/GammW;
}

double BrWhadrons::computeThValue()
{
    double GammW = SM.GammaW();
    double Gammhad = 0.;
    
//  Current SM formula asummes diagonal CKM
    Gammhad += SM.GammaW(SM.getQuarks(StandardModel::UP), SM.getQuarks(StandardModel::DOWN));
    Gammhad += SM.GammaW(SM.getQuarks(StandardModel::CHARM), SM.getQuarks(StandardModel::STRANGE));

    return Gammhad/GammW;
}
