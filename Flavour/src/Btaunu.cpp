/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Btaunu.h"

Btaunu::Btaunu(const StandardModel& SM_i)
: ThObservable(SM_i)
{
    if (SM.ModelName().compare("StandardModel") != 0) std::cout << "\nWARNING: B to tau nu not implemented in: " + SM.ModelName() + " model, returning Standard Model value.\n" << std::endl;
};

double Btaunu::computeThValue()
{
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getMyFlavour()->ComputeCoeffbtaunu();
    double mtau = SM.getLeptons(StandardModel::TAU).getMass();
    double mB = SM.getMesons(QCD::B_P).getMass();
    double mb = SM.getQuarks(QCD::BOTTOM).getMass();
    double fact = 0.989;
    return 1./(64. * M_PI) * mtau * mtau * pow(fact * SM.getMesons(QCD::B_D).getDecayconst(), 2.) * mB * pow(1. - mtau * mtau / mB / mB, 2.) / SM.getMesons(QCD::B_P).computeWidth() * ((*(allcoeff[LO]))(0) + mB * mB/mb/mtau * ((*(allcoeff[LO]))(1) + (*(allcoeff[LO]))(2))).abs2();// PLEASE NOTE THE DECAY CONST
}