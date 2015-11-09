/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDK2.h"

AmpDK2::AmpDK2(const StandardModel& SM_i) 
: mySM(SM_i) 
{}

gslpp::complex AmpDK2::AmpDK(orders order) 
{
    if (mySM.getMyFlavour()->getHDF2().getCoeffK().getOrder() < order % 3)
        throw std::runtime_error("AmpDK::computeThValue(): requires cofficient of order not computed"); 

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffK(
            mySM.getBK().getMu()(0),
            mySM.getBK().getScheme());
            
    gslpp::vector<double> me(mySM.getBK().getBpars());

    double MK = mySM.getMesons(QCD::K_0).getMass();
    double Ms = mySM.Mrun(mySM.getBK().getMu()(0),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBK().getMu()(0),
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
            throw std::runtime_error("AmpDK2::AmpDK(): order not implemented"); 
    }
}

gslpp::complex AmpDK2::AmpMK(orders order) 
{
    if (mySM.getMyFlavour()->getHDF2().getCoeffmK().getOrder() < order % 3)
        throw std::runtime_error("AmpDK::computeThValue(): requires cofficient of order not computed");

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffmK(
            mySM.getBK().getMu()(0),
            mySM.getBK().getScheme());

    gslpp::vector<double> me(mySM.getBK().getBpars());
    double MK = mySM.getMesons(QCD::K_0).getMass();
    double Ms = mySM.Mrun(mySM.getBK().getMu()(0),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBK().getMu()(0),
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

