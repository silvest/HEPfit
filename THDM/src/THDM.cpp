/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <StandardModelMatching.h>
#include "THDM.h"

const std::string THDM::THDMvars[NTHDMvars] = {"logtb","bma","mHh","mA","mHp","m12_2","lambda6","lambda7"};

THDM::THDM() : StandardModel() {   
    mycache = new THDMcache();
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bma", boost::cref(bma)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh", boost::cref(mHh)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mA", boost::cref(mA)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHp", boost::cref(mHp)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m12_2", boost::cref(m12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambda6", boost::cref(lambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambda7", boost::cref(lambda7)));
}

THDM::~THDM(){
    if (IsModelInitialized()) {
            if (myTHDMMatching != NULL) delete(myTHDMMatching);
            if (mycache != NULL) delete(mycache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool THDM::InitializeModel()
{
    myTHDMMatching = new THDMMatching(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelTHDM();
    return(true);
}
    
bool THDM::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool THDM::PreUpdate()
{    
    requireCKM = false;
    requireYe = false;
    requireYn = false;
    
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool THDM::Update(const std::map<std::string, double>& DPars) {
    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);

    return (true);
}

bool THDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);
    
    //In THDM U couple with v2, D with v1 and L with v1
    if (requireYu || requireCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v2() * sqrt(2.));
        Yu = VCKM.transpose()*Yu;
    }
    if (requireYd) {
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / v1() * sqrt(2.));
    }
    if (requireYe) {
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v1() * sqrt(2.));
    }
    if (requireYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v1() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and THDM and THDM-derived parameters in THDMMatching */
    myTHDMMatching->StandardModelMatching::updateSMParameters();
    myTHDMMatching->updateTHDMParameters();

    return (true);
}

double THDM::computeCosa() const
{
    return cos(atan(pow(10.,logtb))-bma);
}

double THDM::computeSina() const
{
    return sin(atan(pow(10.,logtb))-bma);
}

void THDM::setParameter(const std::string name, const double& value){    
    if(name.compare("logtb") == 0) {
        logtb = value;
//        if(tanb > 0.){
        tanb = pow(10.,logtb);
        sinb = tanb / sqrt(1. + tanb*tanb);
        cosb = 1. / sqrt(1. + tanb*tanb);}
//        else {
//            throw std::runtime_error("error in THDM::SetParameter, tanb < 0!"); 
//          }
//        } 
    else if(name.compare("bma") == 0) {
        bma = value;
        sin_ba = sin(bma);
    }
    else if(name.compare("mHh") == 0)
        mHh = value;
    else if(name.compare("mA") == 0)
        mA = value;
    else if(name.compare("mHp") == 0)
        mHp = value;
    else if(name.compare("m12_2") == 0)
        m12_2 = value;
    else if(name.compare("lambda6") == 0)
        lambda6 = value;
    else if(name.compare("lambda7") == 0)
        lambda7 = value;
    else
        StandardModel::setParameter(name,value);
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

///////////////////////////////////////////////////////////////////////////
// Flags

bool THDM::setFlag(const std::string name, const bool value)
{
    bool res = false;
    
    res = StandardModel::setFlag(name,value);

    return(res);
}

double THDM::v1() const {
    return v() * cosb;
}

double THDM::v2() const {
    return v() * sinb;
}

///////////////////////////////////////////////////////////////////////////////

double THDM::obliqueS() const {
  
    complex B00prime_Mz_Mz2_mH_mA;
    complex B00prime_Mz_Mz2_mHp_mHp;
    complex B00prime_Mz_Mz2_mh_mA;
    complex B00prime_Mz_Mz2_Mz_mH;
    complex B00prime_Mz_Mz2_Mz_mh;
    
    complex B0prime_Mz_Mz2_Mz_mH;
    complex B0prime_Mz_Mz2_Mz_mh;
    
    double mh = mHl;
    double Mz2 = Mz*Mz;
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;
    
    B00prime_Mz_Mz2_mH_mA = - mycache->B00_Mz_Mz2_mH_mA(Mz,mHh,mA) + mycache->B00_Mz_0_mH_mA(Mz,mHh,mA);
    B00prime_Mz_Mz2_mHp_mHp = - mycache->B00_Mz_Mz2_mHp_mHp(Mz,mHp) + mycache->B00_Mz_0_mHp_mHp(Mz,mHp);
    B00prime_Mz_Mz2_mh_mA = - mycache->B00_Mz_Mz2_mh_mA(Mz,mh,mA) + mycache->B00_Mz_0_mh_mA(Mz,mh,mA);
    B00prime_Mz_Mz2_Mz_mH = - mycache->B00_Mz_Mz2_Mz_mH(Mz,mHh) + mycache->B00_Mz_0_Mz_mH(Mz,mHh);
    B00prime_Mz_Mz2_Mz_mh = - mycache->B00_Mz_0_Mz_mh(Mz,mh) + mycache->B00_Mz_0_Mz_mh(Mz,mh);
    B0prime_Mz_Mz2_Mz_mH = mycache->B0_Mz_Mz2_Mz_mH(Mz,mHh) - mycache->B0_Mz_0_Mz_mH(Mz,mHh);
    B0prime_Mz_Mz2_Mz_mh = mycache->B0_Mz_Mz2_Mz_mh(Mz,mh) - mycache->B0_Mz_0_Mz_mh(Mz,mh);
    
    double DeltaS = 1./Mz2/M_PI*(sin2_ba * B00prime_Mz_Mz2_mH_mA.real() - B00prime_Mz_Mz2_mHp_mHp.real()
           + cos2_ba * (B00prime_Mz_Mz2_mh_mA.real() + B00prime_Mz_Mz2_Mz_mH.real()
           - B00prime_Mz_Mz2_Mz_mh.real() - Mz2 * B0prime_Mz_Mz2_Mz_mH.real()
           + Mz2 * B0prime_Mz_Mz2_Mz_mh.real()));
    
    return DeltaS;
   
}

double THDM::obliqueT() const {
    
    complex B0_Mz_0_Mz_mH;
    complex B0_Mz_0_Mz_mh;
    complex B0_Mz_0_Mw_mH;
    complex B0_Mz_0_Mw_mh;    
    
    
    double M_w = Mw();
    double mh = mHl;
    double Mz2 = Mz*Mz;
    double Mw2 = M_w*M_w;
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;
    double s_W2 = sW2(); 
    
    B0_Mz_0_Mw_mH = mycache->B0_Mz_0_Mw_mH(Mz,M_w,mHh);
    B0_Mz_0_Mz_mH = mycache->B0_Mz_0_Mz_mH(Mz,mHh);
    B0_Mz_0_Mz_mh = mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh);
    B0_Mz_0_Mw_mh = mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh); 
    
    double DeltaT = 1. / 16. / M_PI / Mw2 / s_W2 * (F(mHp,mA)
           + sin2_ba * (F(mHp,mHh) - F(mA,mHh)) + cos2_ba * (F(mHp,mh) 
           - F(mA,mh) + F(M_w,mHh) - F(M_w,mh) - F(Mz,mHh) 
           + F(Mz,mh) + 4. * Mz2 * (B0_Mz_0_Mz_mH.real() - B0_Mz_0_Mz_mh.real()) 
           - 4. * Mw2 * (B0_Mz_0_Mw_mH.real() - B0_Mz_0_Mw_mh.real()))); 
     
    return DeltaT;
}

double THDM::obliqueU() const {
    
    complex B00prime_Mz_Mw2_mA_mHp;
    complex B00prime_Mz_Mw2_mHp_mHp;
    complex B00prime_Mz_Mw2_mH_mHp;
    complex B00prime_Mz_Mw2_mh_mHp;
    complex B00prime_Mz_Mw2_Mw_mH;
    complex B00prime_Mz_Mw2_Mw_mh;
    
    complex B0prime_Mz_Mw2_Mw_mH;
    complex B0prime_Mz_Mw2_Mw_mh;
    
    double M_w = Mw();
    double mh = mHl;
    double Mz2 = Mz*Mz;
    double Mw2 = M_w*M_w;//
    double sin2_ba = sin_ba*sin_ba;
    double cos2_ba = 1. - sin2_ba;
      
    B00prime_Mz_Mw2_mA_mHp = - mycache->B00_Mz_Mw2_mA_mHp(Mz,M_w,mA,mHp) + mycache->B00_Mz_0_mA_mHp(Mz,mA,mHp);
    B00prime_Mz_Mw2_mHp_mHp = - mycache->B00_Mz_Mw2_mHp_mHp(Mz,M_w,mHp) + mycache->B00_Mz_0_mHp_mHp(Mz,mHp);
    B00prime_Mz_Mw2_mh_mHp = - mycache->B00_Mz_Mw2_mh_mHp(Mz,M_w,mh,mHp) + mycache->B00_Mz_0_mh_mHp(Mz,mh,mHp);
    B00prime_Mz_Mw2_Mw_mH = - mycache->B00_Mz_Mw2_Mw_mH(Mz,M_w,mHh) + mycache->B00_Mz_0_Mw_mH(Mz,M_w,mHh);
    B00prime_Mz_Mw2_Mw_mh = - mycache->B00_Mz_Mw2_Mw_mh(Mz,M_w,mh) + mycache->B00_Mz_0_Mw_mh(Mz,M_w,mh);
    B0prime_Mz_Mw2_Mw_mH = mycache->B0_Mz_Mw2_Mw_mH(Mz,M_w,mHh) - mycache->B0_Mz_0_Mw_mH(Mz,M_w,mHh);
    B0prime_Mz_Mw2_Mw_mh = mycache->B0_Mz_Mw2_Mw_mh(Mz,M_w,mh) - mycache->B0_Mz_0_Mw_mh(Mz,M_w,mh);
    B00prime_Mz_Mw2_mH_mHp = - mycache->B00_Mz_Mw2_mH_mHp(Mz,M_w,mHh,mHp) + mycache->B00_Mz_0_mH_mHp(Mz,mHh,mHp);
    
    double DeltaU = - obliqueS() + 1. / M_PI / Mz2 * (B00prime_Mz_Mw2_mA_mHp.real()
           - 2. * B00prime_Mz_Mw2_mHp_mHp.real() + sin2_ba * B00prime_Mz_Mw2_mH_mHp.real()
           + cos2_ba * (B00prime_Mz_Mw2_mh_mHp.real() + B00prime_Mz_Mw2_Mw_mH.real()
           - B00prime_Mz_Mw2_Mw_mh.real() - Mw2 * B0prime_Mz_Mw2_Mw_mH.real()
           + Mw2 * B0prime_Mz_Mw2_Mw_mh.real()));    
    
    return DeltaU;
 
}

///////////////////////////////////////////////////////////////////////////////

double THDM::F(const double m0, const double m1) const {
    double m12 = m1 * m1;
    double m02 = m0 * m0;
    double F;
    
    if ( m0<=0.0 || m1<0.0 )
        throw std::runtime_error("Invalid argument for THDM::F()\n"); 
    
    if(m0 == 0. && m1 != 0.) {
        F=0.5 * m12;
    } else if(m0 != 0. && m1 == 0.){
        F=0.5 * m02;
    } else if((m0 == 0. && m1 == 0.) || (fabs(m0-m1) < LEPS)){
        F=0.;
    } else if (m0 != 0 && m1 != 0){
        F=0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12);
    } else
        throw std::runtime_error("Error in THDM::F()");
    return (F);
}
