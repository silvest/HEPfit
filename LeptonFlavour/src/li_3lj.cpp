/*
 * Copyright (C) 2015 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "li_3lj.h"
#include "StandardModel.h"

//li_3lj::li_3lj(const StandardModel& SM_i): ThObservable(SM_i)
//{
//};
//
//double li_3lj::computeThValue()
//{
//    return 0.0;
//}

mu_3e::mu_3e(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

double mu_3e::computeThValue()
{
    double alph = mySM.getAle();
    double mE = mySM.getLeptons(StandardModel::ELECTRON).getMass();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    double GammaMU = 2.99598e-19;
    gslpp::vector<gslpp::complex> ** allcoeff_m3e = mySM.getMyLeptonFlavour()->ComputeCoeffli_3lj(1);
//    C_7 = (*(allcoeff_m3e[LO]))(0);
//    C_7p = (*(allcoeff_m3e[LO]))(1);
//    C_9 = (*(allcoeff_m3e[LO]))(2);
//    C_9p = (*(allcoeff_m3e[LO]))(3);
//    C_10 = (*(allcoeff_m3e[LO]))(4);
//    C_10p = (*(allcoeff_m3e[LO]))(5);
//    C_S = (*(allcoeff_m3e[LO]))(6);
//    C_Sp = (*(allcoeff_m3e[LO]))(7);
//    C_P = (*(allcoeff_m3e[LO]))(8);
//    C_Pp = (*(allcoeff_m3e[LO]))(9);
//    C_T = (*(allcoeff_m3e[LO]))(10);
//    C_T5 = (*(allcoeff_m3e[LO]))(11);

//std::cout<<"0"<< (*(allcoeff_m3e[LO]))(0)<<std::endl;
//std::cout<<"1"<< (*(allcoeff_m3e[LO]))(1)<<std::endl;
//std::cout<<"2"<< (*(allcoeff_m3e[LO]))(2)<<std::endl;
//std::cout<<"3"<< (*(allcoeff_m3e[LO]))(3)<<std::endl;
//std::cout<<"4"<< (*(allcoeff_m3e[LO]))(4)<<std::endl;
//std::cout<<"5"<< (*(allcoeff_m3e[LO]))(5)<<std::endl;
//std::cout<<"6"<< (*(allcoeff_m3e[LO]))(6)<<std::endl;
//std::cout<<"7"<< (*(allcoeff_m3e[LO]))(7)<<std::endl;
//std::cout<<"8"<< (*(allcoeff_m3e[LO]))(8)<<std::endl;
//std::cout<<"9"<< (*(allcoeff_m3e[LO]))(9)<<std::endl;
//std::cout<<"10"<< (*(allcoeff_m3e[LO]))(10)<<std::endl;
//std::cout<<"11"<< (*(allcoeff_m3e[LO]))(11)<<std::endl;
//std::cout<<"12"<< (*(allcoeff_m3e[LO]))(12)<<std::endl;
//std::cout<<"13"<< (*(allcoeff_m3e[LO]))(13)<<std::endl;
//std::cout<<"14"<< (*(allcoeff_m3e[LO]))(14)<<std::endl;
//std::cout<<"15"<< (*(allcoeff_m3e[LO]))(15)<<std::endl;
//std::cout<<"16"<< (*(allcoeff_m3e[LO]))(16)<<std::endl;
//std::cout<<"17"<< (*(allcoeff_m3e[LO]))(17)<<std::endl;
//std::cout<<"18"<< (*(allcoeff_m3e[LO]))(18)<<std::endl;
//std::cout<<"19"<< (*(allcoeff_m3e[LO]))(19)<<std::endl;

//    double Brmu3e = ((*(allcoeff_m3e[LO]))(0)* (*(allcoeff_m3e[LO]))(0).conjugate()).abs();
    double Brmu3e = alph*alph/(32.0*M_PI) * pow(mMU,5.0)
                    * ((*(allcoeff_m3e[LO]))(1).abs2()+(*(allcoeff_m3e[LO]))(0).abs2()
                       -4.0*((*(allcoeff_m3e[LO]))(1)*(*(allcoeff_m3e[LO]))(2).conjugate()+(*(allcoeff_m3e[LO]))(3)*(*(allcoeff_m3e[LO]))(0).conjugate()).real()
                       +((*(allcoeff_m3e[LO]))(3).abs2()+(*(allcoeff_m3e[LO]))(2).abs2())*(16.0*log(mMU/mE)-22.0)/3.0
                       +((*(allcoeff_m3e[LO]))(5).abs2()+(*(allcoeff_m3e[LO]))(4).abs2())/6.0
                       +(((*(allcoeff_m3e[LO]))(7)+(*(allcoeff_m3e[LO]))(13)).abs2()+((*(allcoeff_m3e[LO]))(6)+(*(allcoeff_m3e[LO]))(12)).abs2())/3.0
                       +(((*(allcoeff_m3e[LO]))(9)+(*(allcoeff_m3e[LO]))(15)).abs2()+((*(allcoeff_m3e[LO]))(8)+(*(allcoeff_m3e[LO]))(14)).abs2())/24.0
                       +((*(allcoeff_m3e[LO]))(11).abs2()+(*(allcoeff_m3e[LO]))(10).abs2())*6.0
                       -(((*(allcoeff_m3e[LO]))(9)+(*(allcoeff_m3e[LO]))(15))*(*(allcoeff_m3e[LO]))(11).conjugate()+((*(allcoeff_m3e[LO]))(8)+(*(allcoeff_m3e[LO]))(14))*(*(allcoeff_m3e[LO]))(10).conjugate()).real()
                       +((*(allcoeff_m3e[LO]))(1)*(*(allcoeff_m3e[LO]))(5).conjugate()+(*(allcoeff_m3e[LO]))(0)*(*(allcoeff_m3e[LO]))(4).conjugate()
                         +(*(allcoeff_m3e[LO]))(1)*((*(allcoeff_m3e[LO]))(7)+(*(allcoeff_m3e[LO]))(13)).conjugate()+(*(allcoeff_m3e[LO]))(0)*((*(allcoeff_m3e[LO]))(6)+(*(allcoeff_m3e[LO]))(12)).conjugate()).real()*2.0/3.0
                       -4.0*((*(allcoeff_m3e[LO]))(2)*(*(allcoeff_m3e[LO]))(5).conjugate()+(*(allcoeff_m3e[LO]))(3)*(*(allcoeff_m3e[LO]))(4).conjugate()
                             +(*(allcoeff_m3e[LO]))(3)*((*(allcoeff_m3e[LO]))(6)+(*(allcoeff_m3e[LO]))(12)).conjugate()+(*(allcoeff_m3e[LO]))(2)*((*(allcoeff_m3e[LO]))(7)+(*(allcoeff_m3e[LO]))(13)).conjugate()).real()/3.0
                       +2.0*((*(allcoeff_m3e[LO]))(19).abs2()+(*(allcoeff_m3e[LO]))(16).abs2())/3.0
                       +((*(allcoeff_m3e[LO]))(18).abs2()+(*(allcoeff_m3e[LO]))(17).abs2())/3.0
                       +((*(allcoeff_m3e[LO]))(5)*(*(allcoeff_m3e[LO]))(19).conjugate()+(*(allcoeff_m3e[LO]))(4)*(*(allcoeff_m3e[LO]))(16).conjugate()
                         +((*(allcoeff_m3e[LO]))(7)+(*(allcoeff_m3e[LO]))(13))*(*(allcoeff_m3e[LO]))(18).conjugate()+((*(allcoeff_m3e[LO]))(6)+(*(allcoeff_m3e[LO]))(12))*(*(allcoeff_m3e[LO]))(17).conjugate()).real()*2.0/3.0
                       +((*(allcoeff_m3e[LO]))(1)*(*(allcoeff_m3e[LO]))(19).conjugate()+(*(allcoeff_m3e[LO]))(0)*(*(allcoeff_m3e[LO]))(16).conjugate()).real()*4.0/3.0
                       +((*(allcoeff_m3e[LO]))(1)*(*(allcoeff_m3e[LO]))(18).conjugate()+(*(allcoeff_m3e[LO]))(0)*(*(allcoeff_m3e[LO]))(17).conjugate()).real()*2.0/3.0
                       -((*(allcoeff_m3e[LO]))(2)*(*(allcoeff_m3e[LO]))(19).conjugate()+(*(allcoeff_m3e[LO]))(3)*(*(allcoeff_m3e[LO]))(16).conjugate()).real()*8.0/3.0
                       -((*(allcoeff_m3e[LO]))(3)*(*(allcoeff_m3e[LO]))(17).conjugate()+(*(allcoeff_m3e[LO]))(2)*(*(allcoeff_m3e[LO]))(18).conjugate()).real()*4.0/3.0) / GammaMU;

//    double Mu3Erate = (alph*alph/(32.0*M_PI))*(pow(mMU,5.0))*
//                       (VLMuE**2.0 + VRMuE**2.0 - 4.0*(VLMuE*ARMuE + VRMuE*ALMuE)
//                        + (ALMuE**2.0 + ARMuE**2.0)*((16.0/3.0)*log(mMU/(2.0*mE))-14.0/9.0)
//                        + (1.d0/6.0)*(B1LMu3E**2.0 + B1RMu3E**2.0)
//                        + (1.d0/3.0)*(B2LMu3E**2.0 + B2RMu3E**2.0)
//                        + (1.d0/24.0)*(B3LMu3E**2.0 + B3RMu3E**2.0)
//                        + 6.0*(B4LMu3E**2.0 + B4RMu3E**2.0)
//                        - (B3LMu3E*B4LMu3E + B3RMu3E*B4RMu3E)
//                        + (2.0/3.0)*(VLMuE*B1LMu3E + VRMuE*B1RMu3E + VLMuE*B2LMu3E + VRMuE*B2RMu3E)
//                        - (4.0/3.0)*(ARMuE*B1LMu3E + ALMuE*B1RMu3E + ARMuE*B2LMu3E + ALMuE*B2RMu3E)
//                        + (1.0/3.0)*(2.0*(FLLMu3E**2.0 + FRRMu3E**2.0)
//                         + FLRMu3E**2.0 + FRLMu3E**2.0
//                         + 2.0*(B1LMu3E*FLLMu3E + B1RMu3E*FRRMu3E + B2LMu3E*FLRMu3E + B2RMu3E*FRLMu3E)
//                         + 4.0*(VLMuE*FLLMu3E + VRMuE*FRRMu3E)
//                         + 2.0*(VLMuE*FLRMu3E + VRMuE*FRLMu3E) 
//                         - 8.0*(ARMuE*FLLMu3E + ALMuE*FRRMu3E) 
//                         - 4.0*(ALMuE*FRLMu3E + ARMuE*FLRMu3E)))
//
//    double Brmu3e = Mu3Erate/(3.0*(pow(10.0,-19.0)))

                
//
//    C_1 = ((*(allcoeff[LO]))(0) + (*(allcoeff[NLO]))(0));
//    C_2 = ((*(allcoeff[LO]))(1) + (*(allcoeff[NLO]))(1));
//    C_2L = (*(allcoeff[LO]))(1);
//    C_3 = ((*(allcoeff[LO]))(2) + (*(allcoeff[NLO]))(2));
//    C_4 = ((*(allcoeff[LO]))(3) + (*(allcoeff[NLO]))(3));
//    C_5 = ((*(allcoeff[LO]))(4) + (*(allcoeff[NLO]))(4));
//    C_6 = ((*(allcoeff[LO]))(5) + (*(allcoeff[NLO]))(5));
//    C_7 = ((*(allcoeff[LO]))(6) + (*(allcoeff[NLO]))(6));
//    C_8L = (*(allcoeff[LO]))(7);
//    C_9 = ((*(allcoeff[LO]))(8) + (*(allcoeff[NLO]))(8));
//    C_10 = ((*(allcoeff[LO]))(9) + (*(allcoeff[NLO]))(9));
//    C_S = ((*(allcoeff[LO]))(10) + (*(allcoeff[NLO]))(10));
//    C_P = ((*(allcoeff[LO]))(11) + (*(allcoeff[NLO]))(11));
//
    return Brmu3e;
}

tau_3mu::tau_3mu(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

//double tau_3mu::computeThValue()
//{
////    double alph = mySM.getAle();
////    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
////    gslpp::vector<gslpp::complex> ** allcoeff_tm = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(2);
////    double BR_tau_mu_gamma = (alph*pow(mTAU,5.0) * ((*(allcoeff_tm[LO])) * (*(allcoeff_tm[LO])).conjugate()).abs() / (2.26e-12) );
//    return 0.0;
//}

double tau_3mu::computeThValue()
{
    double alph = mySM.getAle();
    double mMU = mySM.getLeptons(StandardModel::MU).getMass();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    double GammaTAU = 2.26735e-12;

    gslpp::vector<gslpp::complex> ** allcoeff_t3m = mySM.getMyLeptonFlavour()->ComputeCoeffli_3lj(2);

    double Brtau3mu = alph*alph/(32.0*M_PI) * pow(mTAU,5.0)
                    * ((*(allcoeff_t3m[LO]))(1).abs2()+(*(allcoeff_t3m[LO]))(0).abs2()
                       -4.0*((*(allcoeff_t3m[LO]))(1)*(*(allcoeff_t3m[LO]))(2).conjugate()+(*(allcoeff_t3m[LO]))(3)*(*(allcoeff_t3m[LO]))(0).conjugate()).real()
                       +((*(allcoeff_t3m[LO]))(3).abs2()+(*(allcoeff_t3m[LO]))(2).abs2())*(16.0*log(mTAU/mMU)-22.0)/3.0
                       +((*(allcoeff_t3m[LO]))(5).abs2()+(*(allcoeff_t3m[LO]))(4).abs2())/6.0
                       +(((*(allcoeff_t3m[LO]))(7)+(*(allcoeff_t3m[LO]))(13)).abs2()+((*(allcoeff_t3m[LO]))(6)+(*(allcoeff_t3m[LO]))(12)).abs2())/3.0
                       +(((*(allcoeff_t3m[LO]))(9)+(*(allcoeff_t3m[LO]))(15)).abs2()+((*(allcoeff_t3m[LO]))(8)+(*(allcoeff_t3m[LO]))(14)).abs2())/24.0
                       +((*(allcoeff_t3m[LO]))(11).abs2()+(*(allcoeff_t3m[LO]))(10).abs2())*6.0
                       -(((*(allcoeff_t3m[LO]))(9)+(*(allcoeff_t3m[LO]))(15))*(*(allcoeff_t3m[LO]))(11).conjugate()+((*(allcoeff_t3m[LO]))(8)+(*(allcoeff_t3m[LO]))(14))*(*(allcoeff_t3m[LO]))(10).conjugate()).real()
                       +((*(allcoeff_t3m[LO]))(1)*(*(allcoeff_t3m[LO]))(5).conjugate()+(*(allcoeff_t3m[LO]))(0)*(*(allcoeff_t3m[LO]))(4).conjugate()
                         +(*(allcoeff_t3m[LO]))(1)*((*(allcoeff_t3m[LO]))(7)+(*(allcoeff_t3m[LO]))(13)).conjugate()+(*(allcoeff_t3m[LO]))(0)*((*(allcoeff_t3m[LO]))(6)+(*(allcoeff_t3m[LO]))(12)).conjugate()).real()*2.0/3.0
                       -4.0*((*(allcoeff_t3m[LO]))(2)*(*(allcoeff_t3m[LO]))(5).conjugate()+(*(allcoeff_t3m[LO]))(3)*(*(allcoeff_t3m[LO]))(4).conjugate()
                             +(*(allcoeff_t3m[LO]))(3)*((*(allcoeff_t3m[LO]))(6)+(*(allcoeff_t3m[LO]))(12)).conjugate()+(*(allcoeff_t3m[LO]))(2)*((*(allcoeff_t3m[LO]))(7)+(*(allcoeff_t3m[LO]))(13)).conjugate()).real()/3.0
                       +2.0*((*(allcoeff_t3m[LO]))(19).abs2()+(*(allcoeff_t3m[LO]))(16).abs2())/3.0
                       +((*(allcoeff_t3m[LO]))(18).abs2()+(*(allcoeff_t3m[LO]))(17).abs2())/3.0
                       +((*(allcoeff_t3m[LO]))(5)*(*(allcoeff_t3m[LO]))(19).conjugate()+(*(allcoeff_t3m[LO]))(4)*(*(allcoeff_t3m[LO]))(16).conjugate()
                         +((*(allcoeff_t3m[LO]))(7)+(*(allcoeff_t3m[LO]))(13))*(*(allcoeff_t3m[LO]))(18).conjugate()+((*(allcoeff_t3m[LO]))(6)+(*(allcoeff_t3m[LO]))(12))*(*(allcoeff_t3m[LO]))(17).conjugate()).real()*2.0/3.0
                       +((*(allcoeff_t3m[LO]))(1)*(*(allcoeff_t3m[LO]))(19).conjugate()+(*(allcoeff_t3m[LO]))(0)*(*(allcoeff_t3m[LO]))(16).conjugate()).real()*4.0/3.0
                       +((*(allcoeff_t3m[LO]))(1)*(*(allcoeff_t3m[LO]))(18).conjugate()+(*(allcoeff_t3m[LO]))(0)*(*(allcoeff_t3m[LO]))(17).conjugate()).real()*2.0/3.0
                       -((*(allcoeff_t3m[LO]))(2)*(*(allcoeff_t3m[LO]))(19).conjugate()+(*(allcoeff_t3m[LO]))(3)*(*(allcoeff_t3m[LO]))(16).conjugate()).real()*8.0/3.0
                       -((*(allcoeff_t3m[LO]))(3)*(*(allcoeff_t3m[LO]))(17).conjugate()+(*(allcoeff_t3m[LO]))(2)*(*(allcoeff_t3m[LO]))(18).conjugate()).real()*4.0/3.0) / GammaTAU;

    return Brtau3mu;
}

tau_3e::tau_3e(const StandardModel& SM_i)
: ThObservable(SM_i), mySM(SM_i)
{}

//double tau_3e::computeThValue()
//{
////    double alph = mySM.getAle();
////    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
////    gslpp::vector<gslpp::complex> ** allcoeff_te = mySM.getMyLeptonFlavour()->ComputeCoeffli_lj_gamma(3);
////    double BR_tau_e_gamma = (alph*pow(mTAU,5.0) * ((*(allcoeff_te[LO])) * (*(allcoeff_te[LO])).conjugate()).abs() / (2.26e-12) );
//    return 0.0;
//}

double tau_3e::computeThValue()
{
    double alph = mySM.getAle();
    double mE = mySM.getLeptons(StandardModel::ELECTRON).getMass();
    double mTAU = mySM.getLeptons(StandardModel::TAU).getMass();
    double GammaTAU = 2.26735e-12;

    gslpp::vector<gslpp::complex> ** allcoeff_t3e = mySM.getMyLeptonFlavour()->ComputeCoeffli_3lj(3);

    double Brtau3e = alph*alph/(32.0*M_PI) * pow(mTAU,5.0)
                    * ((*(allcoeff_t3e[LO]))(1).abs2()+(*(allcoeff_t3e[LO]))(0).abs2()
                       -4.0*((*(allcoeff_t3e[LO]))(1)*(*(allcoeff_t3e[LO]))(2).conjugate()+(*(allcoeff_t3e[LO]))(3)*(*(allcoeff_t3e[LO]))(0).conjugate()).real()
                       +((*(allcoeff_t3e[LO]))(3).abs2()+(*(allcoeff_t3e[LO]))(2).abs2())*(16.0*log(mTAU/mE)-22.0)/3.0
                       +((*(allcoeff_t3e[LO]))(5).abs2()+(*(allcoeff_t3e[LO]))(4).abs2())/6.0
                       +(((*(allcoeff_t3e[LO]))(7)+(*(allcoeff_t3e[LO]))(13)).abs2()+((*(allcoeff_t3e[LO]))(6)+(*(allcoeff_t3e[LO]))(12)).abs2())/3.0
                       +(((*(allcoeff_t3e[LO]))(9)+(*(allcoeff_t3e[LO]))(15)).abs2()+((*(allcoeff_t3e[LO]))(8)+(*(allcoeff_t3e[LO]))(14)).abs2())/24.0
                       +((*(allcoeff_t3e[LO]))(11).abs2()+(*(allcoeff_t3e[LO]))(10).abs2())*6.0
                       -(((*(allcoeff_t3e[LO]))(9)+(*(allcoeff_t3e[LO]))(15))*(*(allcoeff_t3e[LO]))(11).conjugate()+((*(allcoeff_t3e[LO]))(8)+(*(allcoeff_t3e[LO]))(14))*(*(allcoeff_t3e[LO]))(10).conjugate()).real()
                       +((*(allcoeff_t3e[LO]))(1)*(*(allcoeff_t3e[LO]))(5).conjugate()+(*(allcoeff_t3e[LO]))(0)*(*(allcoeff_t3e[LO]))(4).conjugate()
                         +(*(allcoeff_t3e[LO]))(1)*((*(allcoeff_t3e[LO]))(7)+(*(allcoeff_t3e[LO]))(13)).conjugate()+(*(allcoeff_t3e[LO]))(0)*((*(allcoeff_t3e[LO]))(6)+(*(allcoeff_t3e[LO]))(12)).conjugate()).real()*2.0/3.0
                       -4.0*((*(allcoeff_t3e[LO]))(2)*(*(allcoeff_t3e[LO]))(5).conjugate()+(*(allcoeff_t3e[LO]))(3)*(*(allcoeff_t3e[LO]))(4).conjugate()
                             +(*(allcoeff_t3e[LO]))(3)*((*(allcoeff_t3e[LO]))(6)+(*(allcoeff_t3e[LO]))(12)).conjugate()+(*(allcoeff_t3e[LO]))(2)*((*(allcoeff_t3e[LO]))(7)+(*(allcoeff_t3e[LO]))(13)).conjugate()).real()/3.0
                       +2.0*((*(allcoeff_t3e[LO]))(19).abs2()+(*(allcoeff_t3e[LO]))(16).abs2())/3.0
                       +((*(allcoeff_t3e[LO]))(18).abs2()+(*(allcoeff_t3e[LO]))(17).abs2())/3.0
                       +((*(allcoeff_t3e[LO]))(5)*(*(allcoeff_t3e[LO]))(19).conjugate()+(*(allcoeff_t3e[LO]))(4)*(*(allcoeff_t3e[LO]))(16).conjugate()
                         +((*(allcoeff_t3e[LO]))(7)+(*(allcoeff_t3e[LO]))(13))*(*(allcoeff_t3e[LO]))(18).conjugate()+((*(allcoeff_t3e[LO]))(6)+(*(allcoeff_t3e[LO]))(12))*(*(allcoeff_t3e[LO]))(17).conjugate()).real()*2.0/3.0
                       +((*(allcoeff_t3e[LO]))(1)*(*(allcoeff_t3e[LO]))(19).conjugate()+(*(allcoeff_t3e[LO]))(0)*(*(allcoeff_t3e[LO]))(16).conjugate()).real()*4.0/3.0
                       +((*(allcoeff_t3e[LO]))(1)*(*(allcoeff_t3e[LO]))(18).conjugate()+(*(allcoeff_t3e[LO]))(0)*(*(allcoeff_t3e[LO]))(17).conjugate()).real()*2.0/3.0
                       -((*(allcoeff_t3e[LO]))(2)*(*(allcoeff_t3e[LO]))(19).conjugate()+(*(allcoeff_t3e[LO]))(3)*(*(allcoeff_t3e[LO]))(16).conjugate()).real()*8.0/3.0
                       -((*(allcoeff_t3e[LO]))(3)*(*(allcoeff_t3e[LO]))(17).conjugate()+(*(allcoeff_t3e[LO]))(2)*(*(allcoeff_t3e[LO]))(18).conjugate()).real()*4.0/3.0) / GammaTAU;

    return Brtau3e;
}

