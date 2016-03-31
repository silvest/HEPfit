/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDB2.h"

using namespace gslpp;

AmpDB2::AmpDB2(const StandardModel& SM_i) 
: mySM(SM_i) 
{}

gslpp::complex AmpDB2::AmpBd(orders order) 
{
    if (mySM.getMyFlavour()->getHDF2().getCoeffBd().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffBd( 
            mySM.getBBd().getMu(),
        mySM.getBBd().getScheme());
    
    gslpp::vector<double> me(mySM.getBBd().getBpars());
    double MBd = mySM.getMesons(QCD::B_D).getMass();
    double Mb = mySM.Mrun(mySM.getBBd().getMu(),
                mySM.getQuarks(QCD::BOTTOM).getMass_scale(),
                mySM.getQuarks(QCD::BOTTOM).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBBd().getMu(),
                mySM.getQuarks(QCD::DOWN).getMass_scale(),
                mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KBd = MBd/(Mb+Md)*MBd/(Mb+Md);
    double Fb = mySM.getMesons(QCD::B_D).getDecayconst();
    me(0) *= 1./3.*MBd*Fb*Fb;
    me(1) *= -5./24.*KBd*MBd*Fb*Fb;
    me(2) *= 1./24.*KBd*MBd*Fb*Fb;
    me(3) *= 1./4.*KBd*MBd*Fb*Fb;
    me(4) *= 1./12.*KBd*MBd*Fb*Fb;
    
#if SUSYFIT_DEBUG & 1
    std::cout << "Bd: me(0) = " << me(0)  << std::endl;
#endif
#if SUSYFIT_DEBUG & 2
    std::cout << "coefficient Bd: " << (*(allcoeff[LO]) + *(allcoeff[NLO]))(0) << std::endl;
    std::cout << "M: " << me << std::endl;
    std::cout << "mu : " << mySM.getBBd().getMu() << ", mut: " << mySM.getMut() << ", scheme: " << mySM.getBBd().getScheme() << ", B par.: " <<  mySM.getBBd().getBpars()(0) << std::endl;
    std::cout << "U (mut): " << (mySM.getMyFlavour()->getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(),mySM.getMut(),LO)(0,0) +  
            mySM.getMyFlavour()->getHDF2().getUDF2().Df2Evol(mySM.getBBd().getMu(),mySM.getMut(),NLO)(0,0))<< std::endl;
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

gslpp::complex AmpDB2::AmpBs(orders order) 
{
    if (mySM.getMyFlavour()->getHDF2().getCoeffBs().getOrder() < order % 3)
        throw std::runtime_error("DmBd::computeThValue(): requires cofficient of order not computed"); 

    vector<complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffBs(
            mySM.getBBs().getMu(),
            mySM.getBBs().getScheme());

    gslpp::vector<double> me(mySM.getBBs().getBpars());
    double MBs = mySM.getMesons(QCD::B_S).getMass();
    double Mb = mySM.getQuarks(QCD::BOTTOM).getMass();
    double Ms = mySM.Mrun(mySM.getBBs().getMu(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double KBs = MBs/(Mb+Ms)*MBs/(Mb+Ms);
    double Fbs = mySM.getMesons(QCD::B_S).getDecayconst();
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
