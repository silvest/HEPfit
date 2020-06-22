/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "RVllp.h"
#include "StandardModel.h"

double RWmue::computeThValue()
{    
    double Gamml1 = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU));
    double Gamml2 = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON));

    return Gamml1/Gamml2;
}


double RWtaue::computeThValue()
{    
    double Gamml1 = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU));
    double Gamml2 = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_1), SM.getLeptons(StandardModel::ELECTRON));

    return Gamml1/Gamml2;
}

double RWtaumu::computeThValue()
{    
    double Gamml1 = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_3), SM.getLeptons(StandardModel::TAU));
    double Gamml2 = SM.GammaW(SM.getLeptons(StandardModel::NEUTRINO_2), SM.getLeptons(StandardModel::MU));

    return Gamml1/Gamml2;
}




double RZmue::computeThValue()
{    
    double Gamml1 = SM.GammaZ(SM.getLeptons(StandardModel::MU));
    double Gamml2 = SM.GammaZ(SM.getLeptons(StandardModel::ELECTRON));

    return Gamml1/Gamml2;
}


double RZtaue::computeThValue()
{    
    double Gamml1 = SM.GammaZ(SM.getLeptons(StandardModel::TAU));
    double Gamml2 = SM.GammaZ(SM.getLeptons(StandardModel::ELECTRON));

    return Gamml1/Gamml2;
}

double RZtaumu::computeThValue()
{    
    double Gamml1 = SM.GammaZ(SM.getLeptons(StandardModel::TAU));
    double Gamml2 = SM.GammaZ(SM.getLeptons(StandardModel::MU));

    return Gamml1/Gamml2;
}


