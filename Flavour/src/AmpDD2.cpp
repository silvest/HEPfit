/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDD2.h"

AmpDD2::AmpDD2(Flavour& Flavour) : myFlavour(Flavour) {

}

complex AmpDD2::AmpDD(orders order) {
    if (myFlavour.getHDF2().getCoeffDD().getOrder() < order)
        throw std::runtime_error("DmD::computeThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff =  myFlavour.ComputeCoeffdd(
            myFlavour.getModel().getBD().getMu(), 
            myFlavour.getModel().getBD().getScheme()); 
    vector<double> me(myFlavour.getModel().getBD().getBpars()); 
    double MD = myFlavour.getModel().getMesons(QCD::D_0).getMass();
    double Mc = myFlavour.getModel().getQuarks(QCD::CHARM).getMass();
    double Mu = myFlavour.getModel().getQuarks(QCD::UP).getMass();
    double KD = MD/(Mc+Mu)*MD/(Mc+Mu);
    double FD = myFlavour.getModel().getMesons(QCD::D_0).getDecayconst();
    me(0) *= 1./3.*MD*FD*FD;
    me(1) *= -5./24.*KD*MD*FD*FD;
    me(2) *= 1./24.*KD*MD*FD*FD;
    me(3) *= 1./4.*KD*MD*FD*FD;
    me(4) *= 1./12.*KD*MD*FD*FD;
   
    switch(order) {
        case NLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT + myFlavour.getModel().getSM_M12D());
        case LO:
            return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDD2::AmpDD(): order not implemented"); 
    }
}

