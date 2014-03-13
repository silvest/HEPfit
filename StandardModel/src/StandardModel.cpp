/* 
 * Copyright (C) 2012 SusyFit Collaboration
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
    "Mz", "AlsMz", "GF", "ale", "dAle5Mz", "mHl", "delMw", "delSin2th_l", "delGammaZ",
    "mneutrino_1", "mneutrino_2", "mneutrino_3", "melectron", "mmu", "mtau", 
    "lambda", "A", "rhob", "etab", "muw",
    "EpsK", "phiEpsK", "DeltaMK", "KbarEpsK", "Dmk", "SM_M12D"
};

StandardModel::StandardModel() 
: QCD(), VCKM(3, 3, 0.), UPMNS(3, 3, 0.), Yu(3, 3, 0.), Yd(3, 3, 0.), Yn(3, 3, 0.), 
        Ye(3, 3, 0.)
{
    FlagWithoutNonUniversalVC = false;
    FlagNoApproximateGammaZ = false;
    FlagMw = "APPROXIMATEFORMULA";
    FlagRhoZ = "NORESUM";
    FlagKappaZ = "APPROXIMATEFORMULA";
    
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


///////////////////////////////////////////////////////////////////////////
// Initialization

bool StandardModel::InitializeModel()
{
    std::cout << "Model: " << ModelName() << std::endl;
    myEWSM = new EWSM(*this);
    myStandardModelMatching = new StandardModelMatching(*this);
    setModelInitialized(true);
    return(true);
}


///////////////////////////////////////////////////////////////////////////
// Parameters

bool StandardModel::Init(const std::map<std::string, double>& DPars)
{
    std::map<std::string, double> myDPars(DPars);
    myDPars["AlsM"] = myDPars.at("AlsMz");
    myDPars["MAls"] = myDPars.at("Mz");
    Update(myDPars);
    return(CheckParameters(myDPars));
}

bool StandardModel::PreUpdate()
{
    requireCKM = false;
    requireYe = false;
    requireYn = false;
    
    if(!QCD::PreUpdate())  return (false);
    
    return (true);
}

bool StandardModel::Update(const std::map<std::string, double>& DPars)
{
    if(!PreUpdate()) return (false);

    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate())  return (false);
    
    return (true);
}

bool StandardModel::PostUpdate()
{
    if(!QCD::PostUpdate()) return (false);

    /* Set the CKM and PMNS matrices */
    computeCKM();

    /* Set the Yukawa matrices */
    computeYukawas();

    /* Necessary for updating StandardModel parameters in StandardModelMatching */
    if (ModelName()=="StandardModel")
        myStandardModelMatching->updateSMParameters();
    
    return (true);
}

void StandardModel::setParameter(const std::string name, const double& value)
{
    if (name.compare("GF") == 0)
        GF = value;
    else if(name.compare("Mz")==0)
        Mz = value;
    else if(name.compare("AlsMz")==0) {
        AlsMz = value;
        AlsM = value;
    }
    else if (name.compare("ale") == 0)
        ale = value;
    else if (name.compare("dAle5Mz") == 0)
        dAle5Mz = value;    
    else if (name.compare("mHl") == 0)
        mHl = value;
    else if (name.compare("delMw") == 0)
        delMw = value;
    else if (name.compare("delSin2th_l") == 0)
        delSin2th_l = value;
    else if (name.compare("delGammaZ") == 0)
        delGammaZ = value;
    else if (name.compare("mneutrino_1") == 0) {
        leptons[NEUTRINO_1].setMass(value);
        requireYn = true;
    } else if (name.compare("mneutrino_2") == 0) {
        leptons[NEUTRINO_2].setMass(value);
        requireYn = true;
    } else if (name.compare("mneutrino_3") == 0) {
        leptons[NEUTRINO_3].setMass(value);
        requireYn = true;
    } else if (name.compare("melectron") == 0) {
        leptons[ELECTRON].setMass(value);
        requireYe = true;
    } else if (name.compare("mmu") == 0) {
        leptons[MU].setMass(value);
        requireYe = true;
    } else if (name.compare("mtau") == 0) {
        leptons[TAU].setMass(value);
        requireYe = true;
    } else if (name.compare("lambda") == 0) {
        lambda = value;
        requireCKM = true;
    } else if (name.compare("A") == 0) {
        A = value;
        requireCKM = true;
    } else if (name.compare("rhob") == 0) {
        rhob = value;
        requireCKM = true;
    } else if (name.compare("etab") == 0) {
        etab = value;
        requireCKM = true;
    } else if (name.compare("muw") == 0)
        muw = value;
    else if (name.compare("EpsK") == 0)
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
        QCD::setParameter(name, value);
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

void StandardModel::computeCKM()
{
    if (requireCKM) {
        myCKM.setWolfenstein(lambda, A, rhob, etab);
        myCKM.getCKM(VCKM);
    }
    UPMNS = matrix<complex>::Id(3);
}

void StandardModel::computeYukawas()
{
    /* THE FOLLOWING CODES HAVE TO BE MODIFIED!!
     *   The Yukawa matrices have to be computed at a common scale
     *   for all the fermions!!! */
    if (requireYu || requireCKM) {
        Yu = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / v() * sqrt(2.));
        Yu = VCKM.transpose()*Yu;
    }
    if (requireYd) {
        Yd = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yd.assign(i, i, this->quarks[DOWN + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (requireYe) {
        Ye = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (requireYn) {
        Yn = matrix<complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v() * sqrt(2.));
        Yn = Yn * UPMNS.hconjugate();
    }
}


///////////////////////////////////////////////////////////////////////////
// Flags

bool StandardModel::setFlag(const std::string name, const bool value)
{  
    bool res = false;
    if (name.compare("CacheInEWSM") == 0) {
        myEWSM->setFlagCacheInEWSM(value);
        res = true;
    } else if (name.compare("CacheInEWSMcache") == 0) {
        myEWSM->getMyCache()->setFlagCacheInEWSMcache(value);
        res = true;
    } else if (name.compare("WithoutNonUniversalVC") == 0) {
        FlagWithoutNonUniversalVC = value;
        res = true;
    } else if (name.compare("NoApproximateGammaZ") == 0) {
        FlagNoApproximateGammaZ = value;
        res = true;
    } else
        res = QCD::setFlag(name,value);

    return(res);
}

bool StandardModel::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if (name.compare("Mw") == 0) {
        if (myEWSM->checkEWPOscheme(value)) {
            FlagMw = value;
            res = true;
        } else
            throw std::runtime_error("StandardModel::setFlagStr(): Invalid flag "
                                     + name + "=" + value);

    } else if (name.compare("RhoZ") == 0) {
        if (myEWSM->checkEWPOscheme(value)) {
            FlagRhoZ = value;
            res = true;
        } else
            throw std::runtime_error("StandardModel::setFlagStr(): Invalid flag "
                                     + name + "=" + value);
    } else if (name.compare("KappaZ") == 0) {
        if (myEWSM->checkEWPOscheme(value)) {
            FlagKappaZ = value;
            res = true;
        } else
            throw std::runtime_error("StandardModel::setFlagStr(): Invalid flag "
                                     + name + "=" + value);
    } else
        res = QCD::setFlagStr(name,value);

    return(res);
}

bool StandardModel::CheckFlags() const
{
    return(QCD::CheckFlags());
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const 
{
    return ( 1. / sqrt(sqrt(2.) * GF) );
}

double StandardModel::Mw_tree() const 
{
    return ( Mz/sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0*M_PI*ale/sqrt(2.0)/GF/Mz/Mz)) );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::ale_OS(const double mu, orders order) const
{
    if (mu < 50.0) 
        throw std::runtime_error("out of range in StandardModel::ale_OS()"); 
    
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
            throw std::runtime_error("Error in StandardModel::ale_OS()"); 
    }
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

double StandardModel::GammaW() const
{
    return myEWSM->GammaW_SM();
}


////////////////////////////////////////////////////////////////////////
// CKM parameters

// Angles

double StandardModel::computeBeta() const 
{
    return (-VCKM(1, 0) * VCKM(1, 2).conjugate() / (VCKM(2, 0) * VCKM(2, 2).conjugate())).arg();
}

double StandardModel::computeGamma() const 
{
    return (-VCKM(0, 0) * VCKM(0, 2).conjugate() / (VCKM(1, 0) * VCKM(1, 2).conjugate())).arg();
}

double StandardModel::computeAlpha() const 
{
    return (-VCKM(2, 0) * VCKM(2, 2).conjugate() / (VCKM(0, 0) * VCKM(0, 2).conjugate())).arg();
}

double StandardModel::computeBetas() const 
{
    return (-VCKM(2, 1) * VCKM(2, 2).conjugate() / (VCKM(1, 1) * VCKM(1, 2).conjugate())).arg();
}

// Lambda_q

complex StandardModel::computelamt() const 
{
    return VCKM(2, 0) * VCKM(2, 1).conjugate();
}

complex StandardModel::computelamc() const 
{
    return VCKM(1, 0) * VCKM(1, 1).conjugate();
}

complex StandardModel::computelamu() const 
{
    return VCKM(0, 0) * VCKM(0, 1).conjugate();
}

complex StandardModel::computelamt_d() const
{
    return VCKM(2, 0) * VCKM(2, 2).conjugate();
}

complex StandardModel::computelamc_d() const 
{
    return VCKM(1, 0) * VCKM(1, 2).conjugate();
}

complex StandardModel::computelamu_d() const 
{
    return VCKM(0, 0) * VCKM(0, 2).conjugate();
}

complex StandardModel::computelamt_s() const 
{
    return VCKM(2, 1) * VCKM(2, 2).conjugate();
}

complex StandardModel::computelamc_s() const 
{
    return VCKM(1, 1) * VCKM(1, 2).conjugate();
}

complex StandardModel::computelamu_s() const 
{
    return VCKM(0, 1) * VCKM(0, 2).conjugate();
}

double StandardModel::computeRt() const
{
    return (VCKM(2, 0) * VCKM(2, 2).conjugate()
            /(VCKM(1, 0) * VCKM(1, 2).conjugate())).abs();
}

double StandardModel::computeRts() const
{
    return (VCKM(2, 1) * VCKM(2, 2).conjugate()
            /(VCKM(1, 1) * VCKM(1, 2).conjugate())).abs();
}

double StandardModel::computeRb() const
{
    return (VCKM(0, 0) * VCKM(0, 2).conjugate()
            /(VCKM(1, 0) * VCKM(1, 2).conjugate())).abs();
}


