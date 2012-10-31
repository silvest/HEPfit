/* 
 * Copyright (C) 2012 SUSYfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDM.h"
//#include "THDMcache.h"
#include <stdexcept>

const std::string THDM::THDMvars[NTHDMvars] = {"mHp","sin_ba","lambda6","lambda7","mA","m12_2","tanb","mH"};

THDM::THDM() : StandardModel(), mycache() {   
    //mycache = new THDMcache(*this);
}

bool THDM::Update(const std::map<std::string, double>& DPars) {
    computeCKM = false;
    computeYe = false;
    computeYn = false;
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    QCD::Update(DPars);
    if (computeCKM) {
        myCKM.setWolfenstein(lambda, A, rhob, etab);
        myCKM.getCKM(VCKM);
    }
    //In THDM U couple with v2, D with v1 and L with v1
    if (computeYu || computeCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v2() * sqrt(2.));
        Yu = VCKM.transpose()*Yu;
    }
    if (computeYd) {
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / v1() * sqrt(2.));
    }
    if (computeYe) {
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v1() * sqrt(2.));
    }
    if (computeYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v1() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }
}

void THDM::SetParameter(const std::string name, const double& value){    
    if(name.compare("mHp") == 0)
        mHp = value;
    else if(name.compare("tanb") == 0) {
        tanb = value;
        if(tanb > 0.){
        sinb = tanb / sqrt(1. + tanb*tanb);
        cosb = 1. / sqrt(1. + tanb*tanb);}
        else {
            throw std::runtime_error("error in THDM::SetParameter, tanb < 0!"); 
          }
        } 
    else if(name.compare("sin_ba") == 0)
        sin_ba = value;
    else if(name.compare("lambda6") == 0)
        lambda6 = value;
    else if(name.compare("lambda7") == 0)
        lambda7 = value;
    else if(name.compare("mA") == 0)
        mA = value;
    else if(name.compare("m12_2") == 0)
        m12_2 = value;
    else if(name.compare("mH") == 0)
        mH = value;
    else
        StandardModel::SetParameter(name,value);
}

bool THDM::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars));
}

bool THDM::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NTHDMvars; i++) {
        if (DPars.find(THDMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory THDM parameter " << THDMvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}



double THDM::v1() {
    return v() * cosb;
}

double THDM::v2() {
    return v() * sinb;
}

///////////////////////////////////////////////////////////////////////////////

double THDM::obliqueS(){
  
    double Mz2;
    complex B22prime_Mz_Mz2_mH_mA;
    complex B22prime_Mz_Mz2_mHp_mHp;
    complex B22prime_Mz_Mz2_mh_mA;
    complex B22prime_Mz_Mz2_Mz_mH;
    complex B22prime_Mz_Mz2_Mz_mh;
    
    complex B0prime_Mz_Mz2_Mz_mH;
    complex B0prime_Mz_Mz2_Mz_mh;
    
    double mh = mHl;
    Mz2 = Mz*Mz;
    DeltaS = 0.0;
    sin2_ba = sin_ba*sin_ba;
    cos2_ba = 1. - sin2_ba;
    
    B22prime_Mz_Mz2_mH_mA = - mycache.B22_Mz_Mz2_mH_mA(Mz,mH,mA) + mycache.B22_Mz_0_mH_mA(Mz,mH,mA);
    B22prime_Mz_Mz2_mHp_mHp = - mycache.B22_Mz_Mz2_mHp_mHp(Mz,mHp) + mycache.B22_Mz_0_mHp_mHp(Mz,mHp);
    B22prime_Mz_Mz2_mh_mA = - mycache.B22_Mz_Mz2_mh_mA(Mz,mh,mA) + mycache.B22_Mz_0_mh_mA(Mz,mh,mA);
    B22prime_Mz_Mz2_Mz_mH = - mycache.B22_Mz_Mz2_Mz_mH(Mz,mH) + mycache.B22_Mz_0_Mz_mH(Mz,mH);
    B22prime_Mz_Mz2_Mz_mh = - mycache.B22_Mz_0_Mz_mh(Mz,mh) + mycache.B22_Mz_0_Mz_mh(Mz,mh);
    B0prime_Mz_Mz2_Mz_mH = mycache.B0_Mz_Mz2_Mz_mH(Mz,mH) - mycache.B0_Mz_0_Mz_mH(Mz,mH);
    B0prime_Mz_Mz2_Mz_mh = mycache.B0_Mz_Mz2_Mz_mh(Mz,mh) - mycache.B0_Mz_0_Mz_mh(Mz,mh);
    
    DeltaS = 1./Mz2/M_PI*(sin2_ba * B22prime_Mz_Mz2_mH_mA.real() - B22prime_Mz_Mz2_mHp_mHp.real()
           + cos2_ba * (B22prime_Mz_Mz2_mh_mA.real() + B22prime_Mz_Mz2_Mz_mH.real() 
           - B22prime_Mz_Mz2_Mz_mh.real() - Mz2 * B0prime_Mz_Mz2_Mz_mH.real() 
           + Mz2 * B0prime_Mz_Mz2_Mz_mh.real()));
    
    return DeltaS;
   
}

double THDM::obliqueT(){
    
    double Mz2;
    complex B0_Mz_0_Mz_mH;
    complex B0_Mz_0_Mz_mh;
    complex B0_Mz_0_Mw_mH;
    complex B0_Mz_0_Mw_mh;    
    
    
    double M_w = Mw();
    double mh = mHl;
    Mz2 = Mz*Mz;
    Mw2 = M_w*M_w;
    sin2_ba = sin_ba*sin_ba;
    cos2_ba = 1. - sin2_ba;
    s_W2 = sW2(); 
    
    B0_Mz_0_Mw_mH = mycache.B0_Mz_0_Mw_mH(Mz,M_w,mH);
    B0_Mz_0_Mz_mH = mycache.B0_Mz_0_Mz_mH(Mz,mH);
    B0_Mz_0_Mz_mh = mycache.B0_Mz_0_Mw_mh(Mz,M_w,mh);
    B0_Mz_0_Mw_mh = mycache.B0_Mz_0_Mw_mh(Mz,M_w,mh); 
    
    DeltaT = 0.0;
  
    DeltaT = 1. / 16. / M_PI / Mw2 / s_W2 * (PV.F(mHp,mA)
           + sin2_ba * (PV.F(mHp,mH) - PV.F(mA,mH)) + cos2_ba * (PV.F(mHp,mh) 
           - PV.F(mA,mh) + PV.F(M_w,mH) - PV.F(M_w,mh) - PV.F(Mz,mH) 
           + PV.F(Mz,mh) + 4. * Mz2 * (B0_Mz_0_Mz_mH.real() - B0_Mz_0_Mz_mh.real()) 
           - 4. * Mw2 * (B0_Mz_0_Mw_mH.real() - B0_Mz_0_Mw_mh.real()))); 
     
    return DeltaT;
}

double THDM::obliqueU(){
    
    double Mz2;
    complex B22prime_Mz_Mw2_mA_mHp;
    complex B22prime_Mz_Mw2_mHp_mHp;
    complex B22prime_Mz_Mw2_mH_mHp;
    complex B22prime_Mz_Mw2_mh_mHp;
    complex B22prime_Mz_Mw2_Mw_mH;
    complex B22prime_Mz_Mw2_Mw_mh;
    
    complex B0prime_Mz_Mw2_Mw_mH;
    complex B0prime_Mz_Mw2_Mw_mh;
    
    double M_w = Mw();
    double mh = mHl;
    Mz2 = Mz*Mz;
    Mw2 = M_w*M_w;//
    sin2_ba = sin_ba*sin_ba;
    cos2_ba = 1. - sin2_ba;
      
    B22prime_Mz_Mw2_mA_mHp = - mycache.B22_Mz_Mw2_mA_mHp(Mz,M_w,mA,mHp) + mycache.B22_Mz_0_mA_mHp(Mz,mA,mHp);
    B22prime_Mz_Mw2_mHp_mHp = - mycache.B22_Mz_Mw2_mHp_mHp(Mz,M_w,mHp) + mycache.B22_Mz_0_mHp_mHp(Mz,mHp);
    B22prime_Mz_Mw2_mh_mHp = - mycache.B22_Mz_Mw2_mh_mHp(Mz,M_w,mh,mHp) + mycache.B22_Mz_0_mh_mHp(Mz,mh,mHp);
    B22prime_Mz_Mw2_Mw_mH = - mycache.B22_Mz_Mw2_Mw_mH(Mz,M_w,mH) + mycache.B22_Mz_0_Mw_mH(Mz,M_w,mH);
    B22prime_Mz_Mw2_Mw_mh = - mycache.B22_Mz_Mw2_Mw_mh(Mz,M_w,mh) + mycache.B22_Mz_0_Mw_mh(Mz,M_w,mh);
    B0prime_Mz_Mw2_Mw_mH = mycache.B0_Mz_Mw2_Mw_mH(Mz,M_w,mH) - mycache.B0_Mz_0_Mw_mH(Mz,M_w,mH);
    B0prime_Mz_Mw2_Mw_mh = mycache.B0_Mz_Mw2_Mw_mh(Mz,M_w,mh) - mycache.B0_Mz_0_Mw_mh(Mz,M_w,mh);
    B22prime_Mz_Mw2_mH_mHp = - mycache.B22_Mz_Mw2_mH_mHp(Mz,M_w,mH,mHp) + mycache.B22_Mz_0_mH_mHp(Mz,mH,mHp);
    
    DeltaU = 0.0;
    
    DeltaU = - obliqueS() + 1. / M_PI / Mz2 * (B22prime_Mz_Mw2_mA_mHp.real() 
           - 2. * B22prime_Mz_Mw2_mHp_mHp.real() + sin2_ba * B22prime_Mz_Mw2_mH_mHp.real() 
           + cos2_ba * (B22prime_Mz_Mw2_mh_mHp.real() + B22prime_Mz_Mw2_Mw_mH.real() 
           - B22prime_Mz_Mw2_Mw_mh.real() - Mw2 * B0prime_Mz_Mw2_Mw_mH.real() 
           + Mw2 * B0prime_Mz_Mw2_Mw_mh.real()));    
    
    return DeltaU;
 
}
