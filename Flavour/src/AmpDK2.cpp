/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDK2.h"

using namespace gslpp;

AmpDK2::AmpDK2(const StandardModel& SM_i) 
: mySM(SM_i) 
{}

gslpp::complex AmpDK2::AmpDK(orders order) 
{
    if (mySM.getMyFlavour()->getHDF2().getCoeffK().getOrder() < order % 3)
        throw std::runtime_error("AmpDK::computeThValue(): requires cofficient of order not computed"); 

    gslpp::vector<complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffK(
            mySM.getBK().getMu(),
            mySM.getBK().getScheme());
            
    gslpp::vector<double> me(mySM.getBK().getBpars());
#if SUSYFIT_DEBUG & 2
    std::cout << "B-parameter: " << me(0) << std::endl;
#endif

    double MK = mySM.getMesons(QCD::K_0).getMass();
    double Ms = mySM.Mrun(mySM.getBK().getMu(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBK().getMu(),
                mySM.getQuarks(QCD::DOWN).getMass_scale(),
                mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KK = MK/(Ms+Md)*MK/(Ms+Md);
    double FK = mySM.getMesons(QCD::K_0).getDecayconst();
    double mm = MK*FK*FK;
    KK *= mm;
    me(0) *= 1./3.*mm;
    me(1) *= -5./24.*KK;
    me(2) *= 1./24.*KK;
    me(3) *= 1./4.*KK;
    me(4) *= 1./12.*KK;

#if SUSYFIT_DEBUG & 2
    std::cout << "matrix element: " << me(0) << std::endl;
#endif
    switch(order) {
        case FULLNLO:
           return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me);
        case LO:
            return((*(allcoeff[LO])) * me);
        default:
            throw std::runtime_error("AmpDK2::AmpDK(): order not implemented"); 
    }
}

gslpp::complex AmpDK2::AmpMK(orders order) 
{
    if (mySM.getMyFlavour()->getHDF2().getCoeffmK().getOrder() < order % 3)
        throw std::runtime_error("AmpDK::computeThValue(): requires cofficient of order not computed");

    vector<complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffmK(
            mySM.getBK().getMu(),
            mySM.getBK().getScheme());

    gslpp::vector<double> me(mySM.getBK().getBpars());
    double MK = mySM.getMesons(QCD::K_0).getMass();
    double Ms = mySM.Mrun(mySM.getBK().getMu(),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBK().getMu(),
                mySM.getQuarks(QCD::DOWN).getMass_scale(),
                mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    double KK = MK/(Ms+Md)*MK/(Ms+Md);
    double FK = mySM.getMesons(QCD::K_0).getDecayconst();
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

