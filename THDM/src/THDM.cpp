/* 
 * File:   THDM.cpp
 * Author: giovannigrilli
 * 
 * Created on 30 aprile 2012, 15.51
 */

#include "THDM.h"
#include "THDMcache.h"

const std::string THDM::THDMvars[NTHDMvars] = {"mHp","sin_ba","lambda6","lambda7","mA","m12_2","tanb","mH"};

THDM::THDM() : StandardModel() {   
    mycache = new THDMcache(*this);
}

void THDM::Update(const std::map<std::string, double>& DPars) {
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
    else if(name.compare("tanb") == 0)
        tanb = value;
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

void THDM::SetSinb(double sinb){
    sinb = tanb / sqrt(1. + tanb*tanb);
    this->sinb = sinb;
}

void THDM::SetCosb(double sinb){
    cosb = 1. / sqrt(1. + tanb*tanb);
    this->cosb = cosb;
}

double THDM::v1() {
    return v() * cosb;
}

double THDM::v2() {
    return v() * sinb;
}

///////////////////////////////////////////////////////////////////////////////

double THDM::obliqueS(){
  
    mh = mHl;
    Mz2 = Mz*Mz;
    DeltaS = 0.0;
    sin2_ba = sin_ba*sin_ba;
    cos2_ba = 1. - sin2_ba;
    
    B22prime_Mz_Mz2_mH_mA = mycache->B22_Mz_Mz2_mH_mA(Mz,mH,mA) - mycache->B22_Mz_0_mH_mA(Mz,mH,mA);
    B22prime_Mz_Mz2_mHp_mHp = mycache->B22_Mz_Mz2_mHp_mHp(Mz,mHp) - mycache->B22_Mz_0_mHp_mHp(Mz,mHp);
    B22prime_Mz_Mz2_mh_mA = mycache->B22_Mz_Mz2_mh_mA(Mz,mh,mA) - mycache->B22_Mz_0_mh_mA(Mz,mh,mA);
    B22prime_Mz_Mz2_Mz_mH = mycache->B22_Mz_Mz2_Mz_mH(Mz,mH) - mycache->B22_Mz_0_Mz_mH(Mz,mH);
    B22prime_Mz_Mz2_Mz_mh = mycache->B22_Mz_0_Mz_mh(Mz,mh) - mycache->B22_Mz_0_Mz_mh(Mz,mh);
    B0prime_Mz_Mz2_Mz_mH = mycache->B0_Mz_Mz2_Mz_mH(Mz,mH) - mycache->B0_Mz_0_Mz_mH(Mz,mH);
    B0prime_Mz_Mz2_Mz_mh = mycache->B0_Mz_Mz2_Mz_mh(Mz,mh) - mycache->B0_Mz_0_Mz_mh(Mz,mh);
    
    DeltaS = 1./Mz2/M_PI*(sin2_ba * B22prime_Mz_Mz2_mH_mA.real() - B22prime_Mz_Mz2_mHp_mHp.real()
           + cos2_ba * (B22prime_Mz_Mz2_mh_mA.real() + B22prime_Mz_Mz2_Mz_mH.real() 
           - B22prime_Mz_Mz2_Mz_mh.real() - Mz2 * B0prime_Mz_Mz2_Mz_mH.real() 
           + Mz2 * B0prime_Mz_Mz2_Mz_mh.real()));
    
    return DeltaS;
   
}

double THDM::obliqueT(){
    
    mh = mHl;
    Mz2 = Mz*Mz;
    Mw2 = Mw_tree()*Mw_tree();
    sin2_ba = sin_ba*sin_ba;
    cos2_ba = 1. - sin2_ba;
    s_02 = s02(); 
    
    B0_Mz_0_Mw_mH = mycache->B0_Mz_0_Mw_mH(Mz,Mw_tree(),mH);
    B0_Mz_0_Mz_mH = mycache->B0_Mz_0_Mz_mH(Mz,mH);
    B0_Mz_0_Mz_mh = mycache->B0_Mz_0_Mw_mh(Mz,Mw_tree(),mh);
    B0_Mz_0_Mw_mh = mycache->B0_Mz_0_Mw_mh(Mz,Mw_tree(),mh); 
    
    DeltaT = 0.0;
  
    DeltaT = 1. / 16. / M_PI / Mw2 / s_02 * (PV.F(mHp,mA)
           + sin2_ba * (PV.F(mHp,mH) - PV.F(mA,mH)) + cos2_ba * (PV.F(mHp,mh) 
           - PV.F(mA,mh) + PV.F(Mw_tree(),mH) - PV.F(Mw_tree(),mh) - PV.F(Mz,mH) 
           + PV.F(Mz,mh) + 4. * Mz2 * (B0_Mz_0_Mz_mH.real() - B0_Mz_0_Mz_mh.real()) 
           - 4. * Mw2 * (B0_Mz_0_Mw_mH.real() - B0_Mz_0_Mw_mh.real()))); 
     
    return DeltaT;
}

double THDM::obliqueU(){
    
    mh = mHl;
    Mz2 = Mz*Mz;
    Mw2 = Mw_tree()*Mw_tree();//tree
    //s02 = s02(); 
    sin2_ba = sin_ba*sin_ba;
    cos2_ba = 1. - sin2_ba;
      
    B22prime_Mz_Mw2_mA_mHp = mycache->B22_Mz_Mw2_mA_mHp(Mz,Mw_tree(),mA,mHp) - mycache->B22_Mz_0_mA_mHp(Mz,mA,mHp);
    B22prime_Mz_Mw2_mHp_mHp = mycache->B22_Mz_Mw2_mHp_mHp(Mz,Mw_tree(),mHp) - mycache->B22_Mz_0_mHp_mHp(Mz,mHp);
    B22prime_Mz_Mw2_mh_mHp = mycache->B22_Mz_Mw2_mh_mHp(Mz,Mw_tree(),mh,mHp) - mycache->B22_Mz_0_mh_mHp(Mz,mh,mHp);
    B22prime_Mz_Mw2_Mw_mH = mycache->B22_Mz_Mw2_Mw_mH(Mz,Mw_tree(),mH) - mycache->B22_Mz_0_Mw_mH(Mz,Mw_tree(),mH);
    B22prime_Mz_Mw2_Mw_mh = mycache->B22_Mz_Mw2_Mw_mh(Mz,Mw_tree(),mh) - mycache->B22_Mz_0_Mw_mh(Mz,Mw_tree(),mh);
    B0prime_Mz_Mw2_Mw_mH = mycache->B0_Mz_Mw2_Mw_mH(Mz,Mw_tree(),mH) - mycache->B0_Mz_0_Mw_mH(Mz,Mw_tree(),mH);
    B0prime_Mz_Mw2_Mw_mh = mycache->B0_Mz_Mw2_Mw_mh(Mz,Mw_tree(),mh) - mycache->B0_Mz_0_Mw_mh(Mz,Mw_tree(),mh);
    B22prime_Mz_Mw2_mH_mHp = mycache->B22_Mz_Mw2_mH_mHp(Mz,Mw_tree(),mH,mHp) - mycache->B22_Mz_0_mH_mHp(Mz,mH,mHp);
    
    DeltaU = 0.0;
    
    DeltaU = - obliqueS() + 1. / M_PI / Mz2 * (B22prime_Mz_Mw2_mA_mHp.real() 
           - 2. * B22prime_Mz_Mw2_mHp_mHp.real() + sin2_ba * B22prime_Mz_Mw2_mH_mHp.real() 
           + cos2_ba * (B22prime_Mz_Mw2_mh_mHp.real() + B22prime_Mz_Mw2_Mw_mH.real() 
           - B22prime_Mz_Mw2_Mw_mh.real() - Mw2 * B0prime_Mz_Mw2_Mw_mH.real() 
           + Mw2 * B0prime_Mz_Mw2_Mw_mh.real()));    
    
    return DeltaU;
 
}
