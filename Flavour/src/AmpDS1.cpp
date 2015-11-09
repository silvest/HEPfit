/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "AmpDS1.h"
#include <sstream>

AmpDS1::AmpDS1(const StandardModel& SM_i) 
: mySM(SM_i) 
{}

gslpp::complex AmpDS1::AmpDS1pp0(orders order) 
{
    if (mySM.getMyFlavour()->getHDS1().getCoeffDS1PP().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("AmpDK1::computeThValue(): requires cofficient of order" 
                                 + out.str() + "not computed");
    }

    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffDS1PP(
            mySM.getBKd1().getMu()(0),
            mySM.getBKd1().getScheme());

    gslpp::vector<double> me1(mySM.getBKd1().getBpars());
    
    double MK = mySM.getMesons(QCD::K_0).getMass();
    double MP = mySM.getMesons(QCD::P_0).getMass();
    double FK = mySM.getMesons(QCD::K_0).getDecayconst();
    double FP = mySM.getMesons(QCD::P_0).getDecayconst();
    double Ms = mySM.Mrun(mySM.getBKd1().getMu()(0),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBKd1().getMu()(0),
                mySM.getQuarks(QCD::DOWN).getMass_scale(),
                mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    
    double X = sqrt(1.5)*FP*(MK*MK-MP*MP);
    double KK = FP/(FK-FP);
    double Q6b = -4.*sqrt(1.5)*(MK*MK/(Ms+Md))*(MK*MK/(Ms+Md))*FP/KK;
    
    me1(0) *= -1./9.*X;
    me1(1) *= 5./9.*X;
    me1(2) *= 1./3.*X;
    me1(3) = me1(2)+me1(1)-me1(0);
    me1(4) *= Q6b/3.;
    me1(5) *= Q6b;
    me1(6) *= -(1./6.*Q6b*(KK+1)-0.5*X);
    me1(7) *= -(0.5*Q6b*(KK+1)-1./6.*X);
    me1(8) = 1.5*me1(0)-0.5*me1(2);
    me1(9) = me1(1)+0.5*me1(0)-0.5*me1(2);
    me1(10) *= 0.;
    me1(11) *= 0.;

    switch(order) {
        case NLO:
           return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me1);
        case LO:
            return((*(allcoeff[LO])) * me1);
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("AmpDK1::AmpDK(): order " + out.str() + "not implemented");
    }
}

gslpp::complex AmpDS1::AmpDS1pp2(orders order) 
{
    if (mySM.getMyFlavour()->getHDS1().getCoeffDS1PP().getOrder() < order){
        std::stringstream out;
        out << order;
        throw std::runtime_error("AmpDK1::computeThValue(): requires cofficient of "
                                 "order" + out.str() + "not computed");
    }
    
    gslpp::vector<gslpp::complex> ** allcoeff = mySM.getMyFlavour()->ComputeCoeffDS1PP(
            mySM.getBKd3().getMu()(0),
            mySM.getBKd3().getScheme());
    
    gslpp::vector<double> me2(mySM.getBKd3().getBpars());
    
    double MK = mySM.getMesons(QCD::K_0).getMass();
    double MP = mySM.getMesons(QCD::P_0).getMass();
    double FK = mySM.getMesons(QCD::K_0).getDecayconst();
    double FP = mySM.getMesons(QCD::P_0).getDecayconst();
    double Ms = mySM.Mrun(mySM.getBKd3().getMu()(0),
                mySM.getQuarks(QCD::STRANGE).getMass_scale(),
                mySM.getQuarks(QCD::STRANGE).getMass(), FULLNNLO);
    double Md = mySM.Mrun(mySM.getBKd3().getMu()(0),
                mySM.getQuarks(QCD::DOWN).getMass_scale(),
                mySM.getQuarks(QCD::DOWN).getMass(), FULLNNLO);
    
    double X = sqrt(1.5)*FP*(MK*MK-MP*MP);
    double KK = FP/(FK-FP);
    double Q6b = -4.*sqrt(1.5)*(MK*MK/(Ms+Md))*(MK*MK/(Ms+Md))*FP/KK;
    
    me2(0) *= 4./9.*sqrt(2.)*X;
    me2(1) *= me2(0);
    me2(2) *= 0.;
    me2(3) *= 0.;
    me2(4) *= 0.;
    me2(5) *= 0.;
    me2(6) *= -(KK/6./sqrt(2.)*Q6b*(KK+1)-X/sqrt(2.));
    me2(7) *= -(KK/2./sqrt(2.)*Q6b*(KK+1)-sqrt(2.)/6.*X);
    me2(8) = 1.5*me2(0);
    me2(9) = 1.5*me2(0);
    me2(10) = 0.;
    me2(11) = 0.;
    
    switch(order) {
        case NLO:
           return((*(allcoeff[LO]) + *(allcoeff[NLO])) * me2);
        case LO:
            return((*(allcoeff[LO])) * me2);
        default:
            std::stringstream out;
            out << order;
            throw std::runtime_error("AmpDK1::AmpDK(): order " + out.str() + "not implemented");;
    }
}
