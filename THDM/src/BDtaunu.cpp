/* 
 * Copyright (C) 2015 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "BDtaunu.h"
#include "StandardModel.h"

BDtaunu::BDtaunu(const StandardModel& SM_i): 

        ThObservable(SM_i), 
        myTHDM(static_cast<const THDM*> (&SM_i))
{
}

BDtaunu::~BDtaunu()
{}

double BDtaunu::computeThValue()
{
    std::string modelflag=myTHDM->getModelTypeflag();
    if (modelflag != "type2") {
        throw std::runtime_error("flag_model in BDtaunu::computeThValue() can only be \"type2\" at the moment");
    }
   return 0.;
}

RBDtaunu::RBDtaunu(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double RBDtaunu::computeThValue()
{
    double RBDtaunu_SM=myTHDM->getBDtaunu_SM();
    double AD=myTHDM->getBDtaunu_A();
    double BD=myTHDM->getBDtaunu_B();
    double mHp2=myTHDM->getmHp2();
    double tanb=myTHDM->gettanb();
    return RBDtaunu_SM + ( AD + BD * tanb*tanb/(mHp2) ) * tanb*tanb/(mHp2);
}

RBDstartaunu::RBDstartaunu(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double RBDstartaunu::computeThValue()
{
    double RBDstartaunu_SM=myTHDM->getBDstartaunu_SM();
    double ADstar=myTHDM->getBDstartaunu_A();
    double BDstar=myTHDM->getBDstartaunu_B();
    double mHp2=myTHDM->getmHp2();
    double tanb=myTHDM->gettanb();
    return RBDstartaunu_SM + ( ADstar + BDstar * tanb*tanb/(mHp2) ) * tanb*tanb/(mHp2);
}

obsBDtaunu_SM::obsBDtaunu_SM(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double obsBDtaunu_SM::computeThValue()
{
    return myTHDM->getBDtaunu_SM();
}

obsBDtaunu_A::obsBDtaunu_A(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double obsBDtaunu_A::computeThValue()
{
    return myTHDM->getBDtaunu_A();
}

obsBDtaunu_B::obsBDtaunu_B(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double obsBDtaunu_B::computeThValue()
{
    return myTHDM->getBDtaunu_B();
}

obsBDstartaunu_SM::obsBDstartaunu_SM(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double obsBDstartaunu_SM::computeThValue()
{
    return myTHDM->getBDstartaunu_SM();
}

obsBDstartaunu_A::obsBDstartaunu_A(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double obsBDstartaunu_A::computeThValue()
{
    return myTHDM->getBDstartaunu_A();
}

obsBDstartaunu_B::obsBDstartaunu_B(const StandardModel& SM_i)
: BDtaunu(SM_i)
{}

double obsBDstartaunu_B::computeThValue()
{
    return myTHDM->getBDstartaunu_B();
}
