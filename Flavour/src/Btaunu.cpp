/*
 * Copyright (C) 2015 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Btaunu.h"
#include "StandardModel.h"

Btaunu::Btaunu(const StandardModel& SM_i, QCD::meson meson_i)
: ThObservable(SM_i)
{
    meson = meson_i;
    SM.initializeMeson(meson);
};

double Btaunu::computeThValue()
{
    gslpp::vector<gslpp::complex> ** allcoeff = SM.getFlavour().ComputeCoeffbtaunu(meson);
    double mtau = SM.getLeptons(StandardModel::TAU).getMass();
    double mB = SM.getMesons(meson).getMass();
    double mb = SM.getQuarks(QCD::BOTTOM).getMass();
    double fact = 1.; /*factor introduced to scale the decay constant from that of the neutral B to the charged B.*/
    //double fact = 0.989;
    return 1. / (64. * M_PI) * mtau * mtau * pow(fact * SM.getMesons(meson).getDecayconst(), 2.) * mB * pow(1. - mtau * mtau / mB / mB, 2.) / SM.getMesons(meson).computeWidth() * ((*(allcoeff[LO]))(0) 
               + mB * mB / mb / mtau * ((*(allcoeff[LO]))(1) + (*(allcoeff[LO]))(2))).abs2(); // PLEASE NOTE THE DECAY CONST
}