/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDK2.h"

AmpDK2::AmpDK2(Flavour& Flavour) : myFlavour(Flavour) {
}

complex AmpDK2::AmpDK(orders order) {
    if (myFlavour.getHDF2().getCoeffK().getOrder() < order)
        throw std::runtime_error("AmpDK::getThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff = myFlavour.ComputeCoeffK( 
            myFlavour.getModel().getBK().getMu(),
            myFlavour.getModel().getBK().getScheme());
            
    vector<double> me(myFlavour.getModel().getBK().getBpars());
    double MK = myFlavour.getModel().getMesons(QCD::K_0).getMass();
    double Ms = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass(), 3., NNLO);
    double Md = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass(), 3., NNLO);
    double KK = MK/(Ms+Md)*MK/(Ms+Md);
    double FK = myFlavour.getModel().getMesons(QCD::K_0).getDecayconst();
    me(0) *= 1./3.*MK*FK*FK;
    me(1) *= -5./24.*KK*MK*FK*FK;
    me(2) *= 1./24.*KK*MK*FK*FK;
    me(3) *= 1./4.*KK*MK*FK*FK;
    me(4) *= 1./12.*KK*MK*FK*FK;

    switch(order) {
        case NLO:
           return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me);
        //case LO:
        //    return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDK2::AmpDK(): order not implemented"); 
    }
}

complex AmpDK2::AmpMK(orders order) {
    if (myFlavour.getHDF2().getCoeffmK().getOrder() < order)
        throw "AmpDK::getThValue(): requires cofficient of order not computed";

    vector<complex> ** allcoeff = myFlavour.ComputeCoeffmK( 
            myFlavour.getModel().getBK().getMu(),
            myFlavour.getModel().getBK().getScheme());

    vector<double> me(myFlavour.getModel().getBK().getBpars());
    double MK = myFlavour.getModel().getMesons(QCD::K_0).getMass();
    double Ms = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass(), 3);
    double Md = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass(), 3);
    double KK = MK/(Ms+Md)*MK/(Ms+Md);
    double FK = myFlavour.getModel().getMesons(QCD::K_0).getDecayconst();
    me(0) *= 1./3.*MK*FK*FK;
    me(1) *= -5./24.*KK*MK*FK*FK;
    me(2) *= 1./24.*KK*MK*FK*FK;
    me(3) *= 1./4.*KK*MK*FK*FK;
    me(4) *= 1./12.*KK*MK*FK*FK;

    switch(order) {
        case NLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me);
        case LO:
            return((*(allcoeff[LO])) * me);
        default:
            throw "AmpDM2::AmpDK(): order not implemented";
    }
}

