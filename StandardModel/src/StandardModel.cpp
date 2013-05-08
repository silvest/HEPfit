/* 
 * Copyright (C) 2012-2013 SusyFit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <TF1.h>
#include <Math/WrappedTF1.h>
#include <Math/BrentRootFinder.h>
#include <gsl/gsl_sf_zeta.h>
#include "StandardModel.h"
#include "CKM.h"
#include "EWSM.h"
#include "StandardModelMatching.h"


const std::string StandardModel::SMvars[NSMvars] = {
    "ale", "dAle5Mz", "GF", "mHl", "delMw", "delSin2th_l", 
    "delRhoZ_nu", "delRhoZ_e", "delRhoZ_b", "muw",
    "mneutrino_1", "mneutrino_2", "mneutrino_3", "melectron", "mmu", "mtau", 
    "lambda", "A", "rhob", "etab", 
    "EpsK", "phiEpsK", "DeltaMK", "KbarEpsK", "Dmk", "SM_M12D"
};

const std::string StandardModel::SMflags[NSMflags] = {
    "FixedAllSMparams", "EWCHMN", "EWABC", "EWABC2", "EWBURGESS", 
    "withoutNonUniversalVCinEpsilons", "NPZbbbarLinearize", 
    "ApproximateGuOverGb", "ApproximateGdOverGb", "RhoZbFromGqOverGb", 
    "TestSubleadingTwoLoopEW"
};

StandardModel::StandardModel(const bool bDebug_i) 
: QCD(), VCKM(3, 3, 0.), UPMNS(3, 3, 0.), Yu(3, 3, 0.), Yd(3, 3, 0.), Yn(3, 3, 0.), 
        Ye(3, 3, 0.), bDebug(bDebug_i) 
{
    FlagFixedAllSMparams = false;
    FlagEWCHMN = false;
    FlagEWABC = false;
    FlagEWABC2 = false;
    FlagEWBURGESS = false;
    FlagWithoutNonUniversalVC = false;
    FlagApproximateGqOverGb = false;
    FlagRhoZbFromGuOverGb = false;
    FlagRhoZbFromGdOverGb = false;
    FlagTestSubleadingTwoLoopEW = false;
    FlagNPZbbbarLinearize = false;
    
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

bool StandardModel::Init(const std::map<std::string, double>& DPars)
{
    Update(DPars);
    return(CheckParameters(DPars));
}

bool StandardModel::PreUpdate()
{
    computeCKM = false;
    computeYe = false;
    computeYn = false;
    
    if(!QCD::PreUpdate())  return (false);
    
    return (true);
}

bool StandardModel::Update(const std::map<std::string, double>& DPars)
{
    if(!PreUpdate()) return (false);

    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        SetParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate())  return (false);
    
    return (true);
}

bool StandardModel::PostUpdate()
{
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

void StandardModel::SetParameter(const std::string name, const double& value)
{
    if (name.compare("ale") == 0)
        ale = value;
    else if (name.compare("dAle5Mz") == 0)
        dAle5Mz = value;    
    else if (name.compare("GF") == 0)
        GF = value;
    else if (name.compare("mHl") == 0)
        mHl = value;
    else if (name.compare("delMw") == 0)
        delMw = value;
    else if (name.compare("delSin2th_l") == 0)
        delSin2th_l = value;
    else if (name.compare("delRhoZ_nu") == 0)
        delRhoZ_nu = value;
    else if (name.compare("delRhoZ_e") == 0)
        delRhoZ_e = value;
    else if (name.compare("delRhoZ_b") == 0)
        delRhoZ_b = value;
    else if (name.compare("muw") == 0)
        muw = value;
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
    } else if (name.compare("EpsK") == 0)
        EpsK = value;
    else if (name.compare("phiEpsK") == 0)
        phiEpsK = value;
    else if (name.compare("DeltaMK") == 0)
        DeltaMK = value;
    else if (name.compare("KbarEpsK") == 0)
        KbarEpsK = value;
    else if (name.compare("Dmk") == 0)
        Dmk = value;
    else if (name.compare("SM_M12D") == 0)
        SM_M12D = value;
    else
        QCD::SetParameter(name, value);
}

bool StandardModel::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NSMvars; i++) {
        if (DPars.find(SMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory SM parameter " << SMvars[i] << std::endl;
            return false;
        }
    }
    return(QCD::CheckParameters(DPars));
}


///////////////////////////////////////////////////////////////////////////
// Flags

bool StandardModel::SetFlag(const std::string name, const bool& value) 
{  
    bool res = false;
    if (name.compare("FixedAllSMparams") == 0) {
        FlagFixedAllSMparams = value;
        res = true;
    } else if (name.compare("EWCHMN") == 0) {
        FlagEWCHMN = value;
        res = true;
    } else if (name.compare("EWABC") == 0) {
        FlagEWABC = value;
        res = true;
    } else if (name.compare("EWABC2") == 0) {
        FlagEWABC2 = value;
        res = true;
    } else if (name.compare("EWBURGESS") == 0) {
        FlagEWBURGESS = value;
        res = true;
    } else if (name.compare("withoutNonUniversalVCinEpsilons") == 0) {
        FlagWithoutNonUniversalVC = value;
        res = true;
    } else if (name.compare("NPZbbbarLinearize") == 0) {
        FlagNPZbbbarLinearize = value;
        res = true;
    } else if (name.compare("ApproximateGqOverGb") == 0) {
        FlagApproximateGqOverGb = value;
        res = true;
    } else if (name.compare("RhoZbFromGuOverGb") == 0) {
        FlagRhoZbFromGuOverGb = value;
        res = true;
    } else if (name.compare("RhoZbFromGdOverGb") == 0) {
        FlagRhoZbFromGdOverGb = value;
        res = true;
    } else if (name.compare("TestSubleadingTwoLoopEW") == 0) {
        FlagTestSubleadingTwoLoopEW = value;
        res = true;
    }
    return(res);
}


///////////////////////////////////////////////////////////////////////////
// Initialization and Matching

bool StandardModel::InitializeModel() 
{
    myStandardModelMatching = new StandardModelMatching(*this);
    SetModelInitialized(true);
    myEWSM = new EWSM(*this);
    this->SetEWSMflags(*myEWSM);
    return(true);
}

void StandardModel::SetEWSMflags(EWSM& myEWSM) 
{
    myEWSM.setSchemeMw(EWSM::APPROXIMATEFORMULA);
    //myEWSM.setSchemeRhoZ(EWSM::OMSI);
    myEWSM.setSchemeRhoZ(EWSM::NORESUM);
    myEWSM.setSchemeKappaZ(EWSM::APPROXIMATEFORMULA);
    
    /* Should be modified by hand! */
    std::cout << "Schemes for EWPOs:" << std::endl
              << "  Mw: APPROXIMATEFORMULA" << std::endl
              << "  rhoZf: NORESUM" << std::endl
              << "  kappaZf: APPROXIMATEFORMULA" << std::endl;
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const 
{
    return ( 1. / sqrt(sqrt(2.) * GF) );
}

double StandardModel::Mw_tree() const 
{
    double tmp = 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz;
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp)) );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::ale_OS(const double mu, orders order) const
{
    if (mu < 50.0) 
        throw std::runtime_error("out of range in StandardModel::Als_OS()"); 
    
    double N = 20.0/3.0;
    double beta1 = N/3.0;
    double beta2 = N/4.0;
    double alpha_ini = alphaMz();
    double v = 1.0 + 2.0*beta1*alpha_ini/M_PI*log(Mz/mu);

    switch (order) {
        case LO:
            return ( alpha_ini/v );
        case FULLNLO:
            return ( alpha_ini/v*(1.0 - beta2/beta1*alpha_ini/M_PI*log(v)/v) );
        default:
            throw std::runtime_error("Error in StandardModel::Als_OS()"); 
    }
}


////////////////////////////////////////////////////////////////////////

double StandardModel::Mw0() const 
{
    return ( sqrt(c02())*Mz );
}

double StandardModel::s02() const
{
    double tmp = 1.0 - 4.0*M_PI*alphaMz()/sqrt(2.0)/GF/Mz/Mz;
    if (tmp < 0.0)
        throw std::runtime_error("Error in StandardModel::s02()");
    
    return ( ( 1.0 - sqrt(tmp) )/2.0 );
}

double StandardModel::c02() const 
{
    return ( 1.0 - s02() );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::DeltaAlphaLepton(const double s) const 
{
    return myEWSM->DeltaAlphaLepton(s);
}

double StandardModel::DeltaAlphaL5q() const 
{
    return myEWSM->DeltaAlphaL5q();
}
    
double StandardModel::DeltaAlpha() const 
{
    return myEWSM->DeltaAlpha();
}
    
double StandardModel::alphaMz() const 
{
    return myEWSM->alphaMz();
}
    
double StandardModel::Mw() const 
{
    return myEWSM->Mw_SM();
}

double StandardModel::cW2() const 
{
    return myEWSM->cW2_SM();
}

double StandardModel::sW2() const 
{
    return myEWSM->sW2_SM();
}
    
complex StandardModel::rhoZ_l(const lepton l) const
{
    return ( myEWSM->rhoZ_l_SM(l) + myEWSM->delRhoZ_l(l) );
}
    
complex StandardModel::rhoZ_q(const quark q) const 
{
    return ( myEWSM->rhoZ_q_SM(q) + myEWSM->delRhoZ_q(q) );
}
    
complex StandardModel::kappaZ_l(const lepton l) const 
{
    return myEWSM->kappaZ_l_SM(l);
}
    
complex StandardModel::kappaZ_q(const quark q) const 
{
    return myEWSM->kappaZ_q_SM(q);
}

complex StandardModel::gVl(const lepton l) const 
{
    double Ql = getLeptons(l).getCharge();
    return ( gAl(l)
             *(1.0 - 4.0*fabs(Ql)*(myEWSM->kappaZ_l_SM(l))*myEWSM->sW2_SM()) );
}

complex StandardModel::gVq(const quark q) const 
{
    double Qq = getQuarks(q).getCharge();
    return ( gAq(q)
             *(1.0 - 4.0*fabs(Qq)*(myEWSM->kappaZ_q_SM(q))*myEWSM->sW2_SM()) );
}
    
complex StandardModel::gAl(const lepton l) const 
{
    double I3l = getLeptons(l).getIsospin();
    return ( sqrt(myEWSM->rhoZ_l_SM(l) + myEWSM->delRhoZ_l(l))*I3l );
}

complex StandardModel::gAq(const quark q) const 
{
    double I3q = getQuarks(q).getIsospin();
    return ( sqrt(myEWSM->rhoZ_q_SM(q) + myEWSM->delRhoZ_q(q))*I3q );
}

double StandardModel::GammaW() const 
{
    return myEWSM->GammaW_SM();
}

double StandardModel::epsilon1_SM() const 
{
    double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real() + myEWSM->delRhoZ_l(ELECTRON);
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );

    return DeltaRhoPrime;
}

double StandardModel::epsilon2_SM() const 
{
    double s_W2 = myEWSM->sW2_SM(), c_W2 = myEWSM->cW2_SM();
    double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real() + myEWSM->delRhoZ_l(ELECTRON);
    double sin2thetaEff = myEWSM->kappaZ_l_SM(ELECTRON).real()*s_W2;
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    double DeltaKappaPrime = sin2thetaEff/s02() - 1.0;
    double DeltaRW = 1.0 - M_PI*alphaMz()/(sqrt(2.0)*GF*Mz*Mz*s_W2*c_W2);
    
    return ( c02()*DeltaRhoPrime + s02()*DeltaRW/(c02() - s02()) 
             - 2.0*s02()*DeltaKappaPrime );
}

double StandardModel::epsilon3_SM() const 
{
    double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real() + myEWSM->delRhoZ_l(ELECTRON);
    double sin2thetaEff = myEWSM->kappaZ_l_SM(ELECTRON).real()*myEWSM->sW2_SM();
    double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    double DeltaKappaPrime = sin2thetaEff/s02() - 1.0;

    return ( c02()*DeltaRhoPrime + (c02() - s02())*DeltaKappaPrime );
}

double StandardModel::epsilonb_SM() const 
{
    /* epsilon_b from g_A^b
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    //double rhoZe = myEWSM->rhoZ_l_SM(ELECTRON).real() + myEWSM->delRhoZ_l(ELECTRON);
    //double rhoZb = myEWSM->rhoZ_q_SM(BOTTOM).real() + myEWSM->delRhoZ_q(BOTTOM);
    //double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    //double eps1 = DeltaRhoPrime;
    //return ( - 1.0 + sqrt(rhoZb)/(1.0 + eps1/2.0) );

    /* epsilon_b from Re(g_V^b/g_A^b), i.e. Re(kappaZ_b)
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    complex kappaZe = myEWSM->kappaZ_l_SM(ELECTRON);
    complex kappaZb = myEWSM->kappaZ_q_SM(BOTTOM);
    if (FlagWithoutNonUniversalVC) 
        return ( kappaZe.real()/kappaZb.real() - 1.0 ); 
    else 
        return ( (kappaZe.real() + myEWSM->kappaZ_q_SM_FlavorDep(BOTTOM).real())
                 /kappaZb.real() - 1.0 );     
    
    /* epsilon_b from Gamma_b via Eqs.(11), (12) and (16) of IJMP A7, 
     * 1031 (1998) by Altarelli et al. 
     * Note: mb has to be mb=4.7, since Eq.(16) were derived with this value. 
     */
    //double als_Mz = Als(Mz, FULLNNLO);
    //double delta_als = (als_Mz - 0.119)/M_PI;
    //double delta_alpha = (alphaMz() - 1.0/128.90)/ale;
    //double Gamma_b_Born = 0.3798*( 1.0 + delta_als - 0.42*delta_alpha);
    //double a = als_Mz/M_PI;
    //double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    //double mb = Mrun(Mz, getQuarks(BOTTOM).getMass(), FULLNNLO);// This is wrong!
    //double mb = 4.7;
    //std::cout << "mb = " << mb << std::endl;
    //double beta = sqrt(1.0 - 4.0*mb*mb/Mz/Mz);
    //double Nc = 3.0; 
    //double factor = GF*Mz*Mz*Mz/6.0/M_PI/sqrt(2.0);
    //double Gamma_b = factor*beta*((3.0 - beta*beta)/2.0*myEWSM->gVq_SM(BOTTOM).abs2()
    //                              + beta*beta*myEWSM->gAq_SM(BOTTOM).abs2())
    //                 *Nc*RQCD*(1.0 + alphaMz()/12.0/M_PI);
    //return ( (Gamma_b/Gamma_b_Born - 1.0 - 1.42*epsilon1_SM() 
    //          + 0.54*epsilon3_SM() )/2.29 );
}

double StandardModel::epsilon1() const
{ 
    return epsilon1_SM();
}

double StandardModel::epsilon2() const 
{
    return epsilon2_SM();    
}
    
double StandardModel::epsilon3() const 
{
    return epsilon3_SM();
}

double StandardModel::epsilonb() const 
{
    return epsilonb_SM();    
}

double StandardModel::deltaGVb() const 
{
    return 0.0;
}

double StandardModel::deltaGAb() const
{
    return 0.0;
}


////////////////////////////////////////////////////////////////////////
// CKM parameters

// Angles

double StandardModel::getBeta() const 
{
    return (-VCKM(1, 0) * VCKM(1, 2).conjugate() / (VCKM(2, 0) * VCKM(2, 2).conjugate())).arg();
}

double StandardModel::getGamma() const 
{
    return (-VCKM(0, 0) * VCKM(0, 2).conjugate() / (VCKM(1, 0) * VCKM(1, 2).conjugate())).arg();
}

double StandardModel::getAlpha() const 
{
    return (-VCKM(2, 0) * VCKM(2, 2).conjugate() / (VCKM(0, 0) * VCKM(0, 2).conjugate())).arg();
}

double StandardModel::getBetas() const 
{
    return (-VCKM(2, 1) * VCKM(2, 2).conjugate() / (VCKM(1, 1) * VCKM(1, 2).conjugate())).arg();
}

// Lambda_q

complex StandardModel::getlamt() const 
{
    return VCKM(2, 0) * VCKM(2, 1).conjugate();
}

complex StandardModel::getlamc() const 
{
    return VCKM(1, 0) * VCKM(1, 1).conjugate();
}

complex StandardModel::getlamu() const 
{
    return VCKM(0, 0) * VCKM(0, 1).conjugate();
}

complex StandardModel::getlamt_d() const
{
    return VCKM(2, 0) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_d() const 
{
    return VCKM(1, 0) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_d() const 
{
    return VCKM(0, 0) * VCKM(0, 2).conjugate();
}

complex StandardModel::getlamt_s() const 
{
    return VCKM(2, 1) * VCKM(2, 2).conjugate();
}

complex StandardModel::getlamc_s() const 
{
    return VCKM(1, 1) * VCKM(1, 2).conjugate();
}

complex StandardModel::getlamu_s() const 
{
    return VCKM(0, 1) * VCKM(0, 2).conjugate();
}

