/*
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_lj_gamma.h"
#include "StandardModel.h"
// #include "SUSYMatching.h"

li_lj_gamma::li_lj_gamma(const StandardModel& SM_i): ThObservable(SM_i), mySM(SM_i){
//li_lj_gamma::li_lj_gamma(const StandardModel& SM_i, int obsFlag): ThObservable(SM_i), mySM(SM_i){
//    if (obsFlag > 0 and obsFlag < 4) obs = obsFlag;
//    else throw std::runtime_error("li_to_lj can only be"
//            "1 (mu to electron) or"
//            "2 (tau to mu) or"
//            "3 (tau to electron)");
};

double li_lj_gamma::computeThValue(){

    double alph = mySM.getAle();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
 
//    if (obs == 1) {
        gslpp::vector<complex> ** allcoeff = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma();
        return (alph*pow(mMU,5.0) * ((*(allcoeff[LO])) * (*(allcoeff[LO])).conjugate()).abs() / (3.0e-19) );
//        return (alph*ml * ((*(allcoeff[LO])) * (*(allcoeff[LO])).conjugate()).abs());

//        gslpp::vector<complex> C7m2e(2, 0.);
//        C7m2e.assign(0, mySUSYMatching.C7_Lepton(1)(0));
//        C7m2e.assign(1, mySUSYMatching.C7_Lepton(1)(1));
////        C7(obs).assign( 0, mySUSYMatching.C7_Lepton(obs));
////        C7(obs).assign( 1, mySUSYMatching.C7_Lepton(obs));
//        return (alph*pow(mMU,5.0) * ( C7m2e(0).abs2() + C7m2e(1).abs2() ));
//    }
//    if (obs == 2) return(?);
//    if (obs == 3) return(?);
    
//    throw std::runtime_error("EWPO::computeThValue(): Observable type not defined. Can be only any of (1) or (2) or (3)");
//    return (EXIT_FAILURE);

}