/* 
 * File:   AmpDK2.cpp
 * Author: stefano
 * 
 * Created on June 14, 2011, 2:40 PM
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
    double Ms = myFlavour.getModel().getQuarks(QCD::STRANGE).getMass();
    double Md = myFlavour.getModel().getQuarks(QCD::DOWN).getMass();
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

