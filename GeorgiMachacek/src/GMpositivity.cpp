/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GMpositivity.h"
#include "StandardModel.h"

GMpositivity1::GMpositivity1(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity1::computeThValue()
{
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
//    std::cout<<"p1="<<lambda2+lambda3<<std::endl;
    return lambda2+lambda3;
}



GMpositivity2::GMpositivity2(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity2::computeThValue()
{
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;

//    std::cout<<"p2="<<2.0*lambda2+lambda3<<std::endl;
    return 2.0*lambda2+lambda3;
}



GMpositivity3::GMpositivity3(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity3::computeThValue()
{
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double pos3 = -1.0;
    if(lambda1>0 && lambda2+lambda3>0)
    {
        pos3 = -fabs(lambda4)+2.0*sqrt(lambda1*(lambda2+lambda3));
    }
    else
    {
        pos3 = -fabs(lambda1*(lambda2+lambda3));
    }
//    std::cout<<"p3="<<pos3<<std::endl;
    return pos3;
}



GMpositivity4::GMpositivity4(const StandardModel& SM_i)
: ThObservable(SM_i), myGM(static_cast<const GeorgiMachacek&> (SM_i))
{}

double GMpositivity4::computeThValue()
{
    double lambda1 = myGM.getMyGMCache()->lambda1;
    double lambda2 = myGM.getMyGMCache()->lambda2;
    double lambda3 = myGM.getMyGMCache()->lambda3;
    double lambda4 = myGM.getMyGMCache()->lambda4;
    double lambda5 = myGM.getMyGMCache()->lambda5;
    double pos4 = -1.0;
    if(lambda1>0 && 2.0*lambda2+lambda3>0)
    {
        pos4 = lambda4-0.25*fabs(lambda5)+sqrt(2.0*lambda1*(2.0*lambda2+lambda3));
    }
    else
    {
        pos4 = -fabs(lambda1*(2.0*lambda2+lambda3));
    }
//    std::cout<<"p4="<<pos4<<std::endl;
    return pos4;
}
