/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDS1.h"
#include "StandardModel.h"
#include "HeffDS1.h"
#include <sstream>

AmpDS1::AmpDS1(const StandardModel& SM_i) 
: mySM(SM_i) 
{
    mySM.initializeBParameter("BKd1");
    mySM.initializeBParameter("BKd3");
}

gslpp::complex AmpDS1::AmpDS1pp0(orders order) 
{
    if (mySM.getFlavour().getHDS1().getCoeffDS1PP().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("AmpDK1::computeThValue(): requires cofficient of order" 
                                 + out.str() + "not computed");
    }

    gslpp::vector<gslpp::complex> ** allcoeffv = mySM.getFlavour().ComputeCoeffDS1PPv(
            mySM.getBKd1().getMu(),
            mySM.getBKd1().getScheme());
    
    gslpp::vector<gslpp::complex> ** allcoeffz = mySM.getFlavour().ComputeCoeffDS1PPz(
            mySM.getBKd1().getMu(),
            mySM.getBKd1().getScheme());
        
    gslpp::vector<gslpp::complex> allcoeffzLO = (*allcoeffz[LO]) + (*allcoeffz[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffzNLO = (*allcoeffz[NLO]) + (*allcoeffz[NLO_QED11]);
    gslpp::vector<gslpp::complex> allcoeffyLO = (*allcoeffv[LO]) + (*allcoeffv[LO_QED]) - (*allcoeffz[LO]) - (*allcoeffz[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffyNLO = (*allcoeffv[NLO]) + (*allcoeffv[NLO_QED11]) - (*allcoeffz[NLO]) - (*allcoeffz[NLO_QED11]);
    
    gslpp::vector<double> meBKd1(mySM.getBKd1().getBpars());

    switch(order) {
        case NLO:
           return (allcoeffzLO + allcoeffzNLO) * meBKd1 + (allcoeffyLO + allcoeffyNLO) * meBKd1;
        case LO:
            return (allcoeffzLO) * meBKd1 + (allcoeffyLO) * meBKd1;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("AmpDK1::AmpDK(): order " + out.str() + "not implemented");
    }
}

gslpp::complex AmpDS1::AmpDS1pp2(orders order) 
{
    if (mySM.getFlavour().getHDS1().getCoeffDS1PP().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("AmpDK1::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeffv = mySM.getFlavour().ComputeCoeffDS1PPv(
            mySM.getBKd3().getMu(),
            mySM.getBKd3().getScheme());
    
    gslpp::vector<gslpp::complex> ** allcoeffz = mySM.getFlavour().ComputeCoeffDS1PPz(
            mySM.getBKd3().getMu(),
            mySM.getBKd3().getScheme());

    gslpp::vector<gslpp::complex> allcoeffzLO = (*allcoeffz[LO]) + (*allcoeffz[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffzNLO = (*allcoeffz[NLO]) + (*allcoeffz[NLO_QED11]);
    gslpp::vector<gslpp::complex> allcoeffyLO = (*allcoeffv[LO]) + (*allcoeffv[LO_QED]) - (*allcoeffz[LO]) - (*allcoeffz[LO_QED]);
    gslpp::vector<gslpp::complex> allcoeffyNLO = (*allcoeffv[NLO]) + (*allcoeffv[NLO_QED11]) - (*allcoeffz[NLO]) - (*allcoeffz[NLO_QED11]);
    
    gslpp::vector<double> meBKd3(mySM.getBKd3().getBpars());
    
    switch(order) {
        case NLO:
           return (allcoeffzLO + allcoeffzNLO) * meBKd3 + (allcoeffyLO + allcoeffyNLO) * meBKd3;
        case LO:
            return (allcoeffzLO) * meBKd3 + (allcoeffyLO) * meBKd3;
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("AmpDK1::AmpDK(): order " + out.str() + "not implemented");;
    }
}
