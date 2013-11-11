/* 
 * Copyright (C) 2012 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDB2.h"

AmpDB2::AmpDB2(Flavour& Flavour) : myFlavour(Flavour) {
}

complex AmpDB2::AmpBd(orders order) {
    if (myFlavour.getHDF2().getCoeffBd().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff = myFlavour.ComputeCoeffBd( 
            myFlavour.getModel().getBBd().getMu(), 
        myFlavour.getModel().getBBd().getScheme()); 
    
    vector<double> me(myFlavour.getModel().getBBd().getBpars()); 
    double MBd = myFlavour.getModel().getMesons(QCD::B_D).getMass();
    double Mb = myFlavour.getModel().Mrun(myFlavour.getModel().getBBd().getMu(),
                myFlavour.getModel().getQuarks(QCD::BOTTOM).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double Md = myFlavour.getModel().Mrun(myFlavour.getModel().getBBd().getMu(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KBd = MBd/(Mb+Md)*MBd/(Mb+Md);
    double Fb = myFlavour.getModel().getMesons(QCD::B_D).getDecayconst();
    me(0) *= 1./3.*MBd*Fb*Fb;
    me(1) *= -5./24.*KBd*MBd*Fb*Fb;
    me(2) *= 1./24.*KBd*MBd*Fb*Fb;
    me(3) *= 1./4.*KBd*MBd*Fb*Fb;
    me(4) *= 1./12.*KBd*MBd*Fb*Fb;
    
#if SUSYFIT_DEBUG & 1
    std::cout << "Bd: me(0) = " << me(0)  << std::endl;
#endif
#if SUSYFIT_DEBUG & 2
    std::cout << "M: " << me << std::endl;
    std::cout << "M.U: " << myFlavour.getHDF2().getUDF2().Df2Evol(4.2,1.e6,LO).transpose()*me << std::endl;
#endif

    switch(order) {
        case FULLNLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
        case LO:
            return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDB2::AmpBd(): order not implemented"); 
    }
}

complex AmpDB2::AmpBs(orders order) {
    if (myFlavour.getHDF2().getCoeffBs().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff = myFlavour.ComputeCoeffBs(
            myFlavour.getModel().getBBs().getMu(),
            myFlavour.getModel().getBBs().getScheme());

    vector<double> me(myFlavour.getModel().getBBs().getBpars());
    double MBs = myFlavour.getModel().getMesons(QCD::B_S).getMass();
    double Mb = myFlavour.getModel().getQuarks(QCD::BOTTOM).getMass();
    double Ms = myFlavour.getModel().Mrun(myFlavour.getModel().getBBs().getMu(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass_scale(),
                myFlavour.getModel().getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double KBs = MBs/(Mb+Ms)*MBs/(Mb+Ms);
    double Fbs = myFlavour.getModel().getMesons(QCD::B_S).getDecayconst();
    me(0) *= 1./3.*MBs*Fbs*Fbs;
    me(1) *= -5./24.*KBs*MBs*Fbs*Fbs;
    me(2) *= 1./24.*KBs*MBs*Fbs*Fbs;
    me(3) *= 1./4.*KBs*MBs*Fbs*Fbs;
    me(4) *= 1./12.*KBs*MBs*Fbs*Fbs;
#if SUSYFIT_DEBUG & 1
    std::cout << "Bs: me(0) = " << me(0)  << std::endl;
#endif

    
    switch(order) {
        case FULLNLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT);
        case LO:
            return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDB2::AmpBs(): order not implemented"); 
    }
}
