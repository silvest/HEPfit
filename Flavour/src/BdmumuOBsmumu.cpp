/*
 * Copyright (C) 2016 SusyFit Collaboration
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BdmumuOBsmumu.h"
#include "StandardModel.h"

BdmumuOBsmumu::BdmumuOBsmumu(const StandardModel& SM_i) : 
ThObservable(SM_i) { };

double BdmumuOBsmumu::computeThValue() {
    double ratio = SM.getMesons(QCD::B_D).getLifetime()/SM.getMesons(QCD::B_S).getLifetime();
    ratio *= SM.getMesons(QCD::B_D).getMass()/SM.getMesons(QCD::B_S).getMass();
    ratio *= SM.getMesons(QCD::B_D).getDecayconst()/SM.getMesons(QCD::B_S).getDecayconst();
    ratio *= SM.getMesons(QCD::B_D).getDecayconst()/SM.getMesons(QCD::B_S).getDecayconst();
    ratio *= sqrt(1.-4.*SM.getLeptons(StandardModel::MU).getMass()*SM.getLeptons(StandardModel::MU).getMass()/
            SM.getMesons(QCD::B_D).getMass()/SM.getMesons(QCD::B_D).getMass())/
            sqrt(1.-4.*SM.getLeptons(StandardModel::MU).getMass()*SM.getLeptons(StandardModel::MU).getMass()/
            SM.getMesons(QCD::B_S).getMass()/SM.getMesons(QCD::B_S).getMass());
    ratio *= SM.getVCKM()(2,0).abs2()/SM.getVCKM()(2,1).abs2();
    
    return ratio;
}


