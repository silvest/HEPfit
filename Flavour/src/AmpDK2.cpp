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
    if (myFlavour.getHDF2().getCoeffK().getOrder() < order % 3)
        throw std::runtime_error("AmpDK::getThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff = myFlavour.ComputeCoeffK( 
            myFlavour.getModel().getBK().getMu(),
            myFlavour.getModel().getBK().getScheme());
            
    vector<double> me(myFlavour.getModel().getBK().getBpars());

    double MK = myFlavour.getModel().getMesons(QCD::K_0).getMass();
    double Ms = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KK = MK/(Ms+Md)*MK/(Ms+Md);
    double FK = myFlavour.getModel().getMesons(QCD::K_0).getDecayconst();
    double mm = MK*FK*FK;
    KK *= mm;
    me(0) *= 1./3.*mm;
    me(1) *= -5./24.*KK;
    me(2) *= 1./24.*KK;
    me(3) *= 1./4.*KK;
    me(4) *= 1./12.*KK;

    switch(order) {
        case FULLNLO:
           return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me);
        case LO:
            return((*(allcoeff[LO])) * me);
        default:
            throw std::runtime_error("AmpDK2::AmpDK(): order not implemented"); 
    }
}

complex AmpDK2::AmpMK(orders order) {
    if (myFlavour.getHDF2().getCoeffmK().getOrder() < order % 3)
        throw std::runtime_error("AmpDK::getThValue(): requires cofficient of order not computed");

    vector<complex> ** allcoeff = myFlavour.ComputeCoeffmK( 
            myFlavour.getModel().getBK().getMu(),
            myFlavour.getModel().getBK().getScheme());

    vector<double> me(myFlavour.getModel().getBK().getBpars());
    double MK = myFlavour.getModel().getMesons(QCD::K_0).getMass();
    double Ms = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = myFlavour.getModel().Mrun(myFlavour.getModel().getBK().getMu(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KK = MK/(Ms+Md)*MK/(Ms+Md);
    double FK = myFlavour.getModel().getMesons(QCD::K_0).getDecayconst();
    double mm = MK*FK*FK;
    KK *= mm;
    me(0) *= 1./3.*mm;
    me(1) *= -5./24.*KK;
    me(2) *= 1./24.*KK;
    me(3) *= 1./4.*KK;
    me(4) *= 1./12.*KK;

    switch(order) {
        case FULLNLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me);
        case LO:
            return((*(allcoeff[LO])) * me);
        default:
            throw "AmpDM2::AmpDK(): order not implemented";
    }
}

