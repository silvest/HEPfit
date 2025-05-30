/*
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "Phis_JPsiPhi.h"
#include "StandardModel.h"
#include "AmpDB2.h"

Phis_JPsiPhi::Phis_JPsiPhi(const StandardModel& SM_i) : ThObservable(SM_i){
            SM.getFlavour().getDB2(1);
};

double Phis_JPsiPhi::computeThValue() 
{
    return (-remainder((SM.getFlavour().getDB2(1).getM21(FULLNLO).arg() - 2. * SM.getCKM().computelamc_s().arg() + 2.*SM.getPhiBs() ),2.*M_PI));
}
