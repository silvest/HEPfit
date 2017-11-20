/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDD2.h"
#include "StandardModel.h"
#include "HeffDF2.h"

AmpDD2::AmpDD2(const StandardModel& SM_i) 
: mySM(SM_i) 
{
    mySM.initializeBParameter("BD");
}

gslpp::complex AmpDD2::AmpDD(orders order) 
{
    if (mySM.getFlavour().getHDF2().getCoeffDD().getOrder() < order)
        throw std::runtime_error("DmD::computeThValue(): requires cofficient of order not computed"); 

    gslpp::vector<gslpp::complex> ** allcoeff =  mySM.getFlavour().ComputeCoeffdd(
            mySM.getBD().getMu(),
            mySM.getBD().getScheme());
    gslpp::vector<double> me(mySM.getBD().getBpars());
    double MD = mySM.getMesons(QCD::D_0).getMass();
    double Mc = mySM.getQuarks(QCD::CHARM).getMass();
    double Mu = mySM.getQuarks(QCD::UP).getMass();
    double KD = MD/(Mc+Mu)*MD/(Mc+Mu);
    double FD = mySM.getMesons(QCD::D_0).getDecayconst();
    me(0) *= 1./3.*MD*FD*FD;
    me(1) *= -5./24.*KD*MD*FD*FD;
    me(2) *= 1./24.*KD*MD*FD*FD;
    me(3) *= 1./4.*KD*MD*FD*FD;
    me(4) *= 1./12.*KD*MD*FD*FD;
   
    switch(order) {
        case NLO:
            return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me / HCUT + mySM.getOptionalParameter("SM_M12D"));
        case LO:
            return((*(allcoeff[LO])) * me / HCUT);
        default:
            throw std::runtime_error("AmpDD2::AmpDD(): order not implemented"); 
    }
}

