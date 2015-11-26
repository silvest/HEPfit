/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "mueconversion.h"
#include "StandardModel.h"

//mueconversion::mueconversion(const StandardModel& SM_i): ThObservable(SM_i)
//{
//};
//
//double mueconversion::computeThValue()
//{
//    return 0.0;
//}

mueconversion_Ti::mueconversion_Ti(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double mueconversion_Ti::computeThValue()
{
    double alph = mySM.getAle();
//    double mE = mySM.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
//    double GammaMU = 2.99598e-19;
    gslpp::vector<gslpp::complex> ** allcoeff_mueconv = mySM.getMyLeptonFlavour()->ComputeCoeffmueconversion();
    double ZTi=22.0;
    double NTi=26.0;
    double Zeff=17.6;
    double Fq=0.54;
    double Gammamue_Ti = 4.0*pow(alph,5)*pow(Zeff,4)/ZTi*Fq*Fq*pow(mMU,5)
                         *((ZTi*((*(allcoeff_mueconv[LO]))(1)-(*(allcoeff_mueconv[LO]))(2))
                            -(2.0*ZTi+NTi)*(*(allcoeff_mueconv[LO]))(5)
                            -(ZTi+2.0*NTi)*(*(allcoeff_mueconv[LO]))(7)).abs2()
                           +(ZTi*((*(allcoeff_mueconv[LO]))(0)-(*(allcoeff_mueconv[LO]))(3))
                             -(2.0*ZTi+NTi)*(*(allcoeff_mueconv[LO]))(4)
                             -(ZTi+2.0*NTi)*(*(allcoeff_mueconv[LO]))(6)).abs2());

    return Gammamue_Ti;
}
