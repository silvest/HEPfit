/* 
 * File:   StandardModel.cpp
 * Author: silvest
 * 
 * Created on November 30, 2010, 1:27 PM
 */

//#include <bt/assign/list_oost/assign/list_of.hpp> // for 'map_list_of()'
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include "StandardModel.h"
#include "CKM.h"
#include "EWSM.h"
#include "StandardModelMatching.h"


const std::string StandardModel::SMvars[NSMvars] = {"GF", "mneutrino_1", "mneutrino_2",
    "mneutrino_3", "melectron", "mmu", "mtau", "lambda", "A", "rhob", "etab", "ale",
    "dAle5Mz", "mHl", "muw", "phiEpsK","DeltaMK", "KbarEpsK", "Dmk", "SM_M12D" };

StandardModel::StandardModel(const bool bDebug_i) : QCD(), VCKM(3, 3, 0.), UPMNS(3, 3, 0.), Yu(3, 3, 0.),
Yd(3, 3, 0.), Yn(3, 3, 0.), Ye(3, 3, 0.) {
    leptons[NEUTRINO_1].setCharge(0.);
    leptons[NEUTRINO_2].setCharge(0.);    
    leptons[NEUTRINO_3].setCharge(0.);    
    leptons[ELECTRON].setCharge(-1.);
    leptons[MU].setCharge(-1.);    
    leptons[TAU].setCharge(-1.); 
    leptons[NEUTRINO_1].setIsospin(1./2.);
    leptons[NEUTRINO_2].setIsospin(1./2.);    
    leptons[NEUTRINO_3].setIsospin(1./2.);    
    leptons[ELECTRON].setIsospin(-1./2.);
    leptons[MU].setIsospin(-1./2.);   
    leptons[TAU].setIsospin(-1./2.);
}

bool StandardModel::SetFlag(const std::string name , const bool& value){  
    return (false);
}

bool StandardModel::Init(const std::map<std::string, double>& DPars) {
    Update(DPars);
    return(CheckParameters(DPars));
}

bool StandardModel::PreUpdate(){
    computeCKM = false;
    computeYe = false;
    computeYn = false;
    
    if(!QCD::PreUpdate())  return (false);
    
    return (true);
}

bool StandardModel::Update(const std::map<std::string, double>& DPars) {
    if(!PreUpdate()) return (false);
    
    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate())  return (false);
    
    return (true);
}

bool StandardModel::PostUpdate(){
    if(!QCD::PostUpdate()) return (false);
    
    if (computeCKM) {
        myCKM.setWolfenstein(lambda, A, rhob, etab);
        myCKM.getCKM(VCKM);
    }
    if (computeYu || computeCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v() * sqrt(2.));
        Yu = VCKM.transpose()*Yu;
    }
    if (computeYd) {
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (computeYe) {
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (computeYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }
    
     return (true);
}

void StandardModel::SetParameter(const std::string name, const double& value) {
    if (name.compare("GF") == 0)
        GF = value;
    else if (name.compare("ale") == 0)
        ale = value;
    else if (name.compare("dAle5Mz") == 0)
        dAle5Mz = value;
    else if (name.compare("mHl") == 0)
        mHl = value;
    else if (name.compare("muw") == 0)
        muw = value;
    else if (name.compare("SM_M12D") == 0)
        SM_M12D = value;
    else if (name.compare("phiEpsK") == 0)
        phiEpsK = value;
    else if (name.compare("KbarEpsK") == 0)
        KbarEpsK = value;
    else if (name.compare("Dmk") == 0)
        Dmk = value;
    else if (name.compare("DeltaMK") == 0)
        DeltaMK = value;
    else if (name.compare("mneutrino_1") == 0) {
        leptons[NEUTRINO_1].setMass(value);
        computeYn = true;
    } else if (name.compare("mneutrino_2") == 0) {
        leptons[NEUTRINO_2].setMass(value);
        computeYn = true;
    } else if (name.compare("mneutrino_3") == 0) {
        leptons[NEUTRINO_3].setMass(value);
        computeYn = true;
    } else if (name.compare("melectron") == 0) {
        leptons[ELECTRON].setMass(value);
        computeYe = true;
    } else if (name.compare("mmu") == 0) {
        leptons[MU].setMass(value);
        computeYe = true;
    } else if (name.compare("mtau") == 0) {
        leptons[TAU].setMass(value);
        computeYe = true;
    } else if (name.compare("lambda") == 0) {
        lambda = value;
        computeCKM = true;
    } else if (name.compare("A") == 0) {
        A = value;
        computeCKM = true;
    } else if (name.compare("rhob") == 0) {
        rhob = value;
        computeCKM = true;
    } else if (name.compare("etab") == 0) {
        etab = value;
        computeCKM = true;
    } else
        QCD::SetParameter(name, value);
}

bool StandardModel::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NSMvars; i++) {
        if (DPars.find(SMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory SM parameter " << SMvars[i] << std::endl;
            return false;
        }
    }
    return(QCD::CheckParameters(DPars));
}


///////////////////////////////////////////////////////////////////////////
// Matching

bool StandardModel::InitializeModel(){
    
    myStandardModelMatching = new StandardModelMatching(*this);
    SetModelInitialized(true);
    myEWSM = new EWSM(*this);
    return(true);
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const {
    return ( 1. / sqrt(sqrt(2.) * GF) );
}

double StandardModel::Mw_tree() const {
    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz;
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp)) );
}

double StandardModel::Mw0() const {
    return ( sqrt(c02())*Mz );
}

double StandardModel::s02() const {
    return ( ( 1.0 - sqrt(1.0 - 4.0*M_PI*alphaMz()/sqrt(2.0)/GF/Mz/Mz) )/2.0 );
}

double StandardModel::c02() const {
    return ( 1.0 - s02() );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::DeltaAlphaLepton() const {
    return myEWSM->DeltaAlphaLepton();
}

double StandardModel::DeltaAlphaL5q() const {
    return myEWSM->DeltaAlphaL5q();
}
    
double StandardModel::DeltaAlpha() const {
    return myEWSM->DeltaAlpha();
}
    
double StandardModel::alphaMz() const {
    return myEWSM->alphaMz();
}
    
double StandardModel::Mw() const {
    return myEWSM->Mw_SM();
}

double StandardModel::cW2() const {
    return myEWSM->cW2_SM();
}

double StandardModel::sW2() const {
    return myEWSM->sW2_SM();
}
    
complex StandardModel::rhoZ_l(const lepton l) const {
    return myEWSM->rhoZ_l_SM(l);
}
    
complex StandardModel::rhoZ_q(const quark q) const {
    return myEWSM->rhoZ_q_SM(q);
}
    
complex StandardModel::kappaZ_l(const lepton l) const {
    return myEWSM->kappaZ_l_SM(l);
}
    
complex StandardModel::kappaZ_q(const quark q) const {
    return myEWSM->kappaZ_q_SM(q);
}

complex StandardModel::gVl(const lepton l) const {
    return myEWSM->gVl_SM(l);
}

complex StandardModel::gVq(const quark q) const {
    return myEWSM->gVq_SM(q);
}
    
complex StandardModel::gAl(const lepton l) const {
    return myEWSM->gAl_SM(l);
}

complex StandardModel::gAq(const quark q) const {
    return myEWSM->gAq_SM(q);
}

double StandardModel::GammaW() const {
    return myEWSM->GammaW_SM();
}

double StandardModel::sigma_l_LEP2(const StandardModel::lepton l, const double s,
                                   const double Mw, const double GammaZ,
                                   const bool bDP, const bool bWEAK, const bool bQED) const {
    return (myEWSM->sigma_l(l, s, Mw, GammaZ, bDP, bWEAK, bQED));
}

double StandardModel::sigma_q_LEP2(const StandardModel::quark q, const double s,
                                   const double Mw, const double GammaZ,
                                   const bool bDP, const bool bWEAK, const bool bQED) const {
    return (myEWSM->sigma_q(q, s, Mw, GammaZ, bDP, bWEAK, bQED));
}

double StandardModel::AFB_l_LEP2(const StandardModel::lepton l, const double s,
                                 const double Mw, const double GammaZ,
                                 const bool bDP, const bool bWEAK, const bool bQED) const {
    return (myEWSM->AFB_l(l, s, Mw, GammaZ, bDP, bWEAK, bQED));
}

double StandardModel::AFB_q_LEP2(const StandardModel::quark q, const double s,
                                 const double Mw, const double GammaZ,
                                 const bool bDP, const bool bWEAK, const bool bQED) const {
    return (myEWSM->AFB_q(q, s, Mw, GammaZ, bDP, bWEAK, bQED));
}

double StandardModel::DsigmaLEP2_l(const StandardModel::lepton l, const double s, const double cos_theta,  
                                   const double W, const double X, const double Y, const double GammaZ) const{
    return (myEWSM->dsigmaLEP2_l(l, s, Mw(), cos_theta, W, X, Y, GammaZ));
}

double StandardModel::DsigmaLEP2_q(const StandardModel::quark q, const double s, 
                                   const double cos_theta, const double W, 
                                   const double X, const double Y, const double GammaZ) const{
    return (myEWSM->dsigmaLEP2_q(q, s, Mw(), cos_theta, W, X, Y, GammaZ));
}
    

////////////////////////////////////////////////////////////////////////
// CKM parameters

// Angles

double StandardModel::getBeta() const {
    return (-VCKM(1, 0) * VCKM(1, 2).conjugate() / (VCKM(2, 0) * VCKM(2, 2).conjugate())).arg();
}

double StandardModel::getGamma() const {
    return (-VCKM(0, 0) * VCKM(0, 2).conjugate() / (VCKM(1, 0) * VCKM(1, 2).conjugate())).arg();
}

double StandardModel::getAlpha() const {
    return (-VCKM(2, 0) * VCKM(2, 2).conjugate() / (VCKM(0, 0) * VCKM(0, 2).conjugate())).arg();
}

double StandardModel::getBetas() const {
    return (-VCKM(2, 1) * VCKM(2, 2).conjugate() / (VCKM(1, 1) * VCKM(1, 2).conjugate())).arg();
}

// Lambda_q

complex StandardModel::getlamt() const {
    return VCKM(2, 0) * VCKM(2, 1).conjugate();
}

complex StandardModel::getlamc() const {
    return VCKM(1, 0) * VCKM(1, 1).conjugate();
}

complex StandardModel::getlamu() const {
    return VCKM(0, 0) * VCKM(0, 1).conjugate();
}

complex StandardModel::getlamt_d() const {
    return VCKM(2, 0) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_d() const {
    return VCKM(1, 0) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_d() const {
    return VCKM(0, 0) * VCKM(0, 2).conjugate();
}

complex StandardModel::getlamt_s() const {
    return VCKM(2, 1) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_s() const {
    return VCKM(1, 1) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_s() const {
    return VCKM(0, 1) * VCKM(0, 2).conjugate();
}

