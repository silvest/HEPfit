/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <gsl/gsl_sf_zeta.h>
#include <algorithm>
#include "StandardModel.h"
#include "EWSMcache.h"
#include "EWSMOneLoopEW.h"
#include "EWSMTwoLoopQCD.h"
#include "EWSMThreeLoopQCD.h"
#include "EWSMTwoLoopEW.h"
#include "EWSMThreeLoopEW2QCD.h"
#include "EWSMThreeLoopEW.h"
#include "EWSMApproximateFormulae.h"
#include "LeptonFlavour.h"
#include "gslpp_function_adapter.h"
/* BEGIN: REMOVE FROM THE PACKAGE */
#include "EWSMTwoFermionsLEP2.h"
#include <functional>
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;
/* END: REMOVE FROM THE PACKAGE */
  
std::string StandardModel::SMvars[NSMvars] = {
    "lambda", "A", "rhob", "etab", "Mz", "AlsMz", "GF", "ale", "dAle5Mz", "mHl", 
    "delMw", "delSin2th_l", "delSin2th_q", "delSin2th_b", "delGammaZ", "delsigma0H", "delR0l", "delR0c", "delR0b",
    "mneutrino_1", "mneutrino_2", "mneutrino_3", "melectron", "mmu", "mtau", "muw"
};

const double StandardModel::GeVminus2_to_nb = 389379.338;
const double StandardModel::Mw_error = 0.00001; /* 0.01 MeV */

StandardModel::StandardModel()
: QCD(), Yu(3, 3, 0.), Yd(3, 3, 0.), Yn(3, 3, 0.),
Ye(3, 3, 0.), SMM(*this), SMFlavour(*this)
{
    FlagWithoutNonUniversalVC = false;
    FlagNoApproximateGammaZ = false;
    FlagMw = "APPROXIMATEFORMULA";
    FlagRhoZ = "NORESUM";
    FlagKappaZ = "APPROXIMATEFORMULA";
    FlagWolfenstein = true;
    FlagUseVud = false;
    FlagFixMuwMut = false;

    FlagMWinput = false;
    
    FlagSMAux = false;

    /* Internal flags for EWPO (for debugging) */
    flag_order[EW1] = true;
    flag_order[EW1QCD1] = true;
    flag_order[EW1QCD2] = true;
    flag_order[EW2] = true;
    flag_order[EW2QCD1] = true;
    flag_order[EW3] = true;

    //Flags for LEP2 observables
    flagLEP2[Weak] = true;
    flagLEP2[WeakBox] = true;
    flagLEP2[ISR] = true;
    flagLEP2[QEDFSR] = true;
    flagLEP2[QCDFSR] = true;
    
    bSigmaForAFB = false;
    bSigmaForR = false;
    
    // Caches for EWPO
    FlagCacheInStandardModel = true; // use caches in the current class
    useDeltaAlphaLepton_cache = false;
    useDeltaAlpha_cache = false;
    useMw_cache = false;
    useGammaW_cache = false;
    DeltaAlphaLepton_cache = 0.0;
    DeltaAlpha_cache = 0.0;
    Mw_cache = 0.0;
    GammaW_cache = 0.0;
    for (int i = 0; i < 12; ++i) {
        useRhoZ_f_cache[i] = false;
        useKappaZ_f_cache[i] = false;
        rhoZ_f_cache[i] = gslpp::complex(0.0, 0.0, false);
        kappaZ_f_cache[i] = gslpp::complex(0.0, 0.0, false);
    }

    myEWSMcache = NULL;
    myOneLoopEW = NULL;
    myTwoLoopQCD = NULL;
    myThreeLoopQCD = NULL;
    myTwoLoopEW = NULL;
    myThreeLoopEW2QCD = NULL;
    myThreeLoopEW = NULL;
    myApproximateFormulae = NULL;
    /* BEGIN: REMOVE FROM THE PACKAGE */
    myTwoFermionsLEP2 = NULL;
    /* END: REMOVE FROM THE PACKAGE */

    // Particle(std::string name, double mass, double mass_scale = 0., double width = 0., double charge = 0.,double isospin = 0.);
    leptons[NEUTRINO_1] = Particle("NEUTRINO_1", 0., 0., 0., 0., .5);
    leptons[NEUTRINO_2] = Particle("NEUTRINO_2", 0., 0., 0., 0., .5);
    leptons[NEUTRINO_3] = Particle("NEUTRINO_3", 0., 0., 0., 0., .5);
    leptons[ELECTRON] = Particle("ELECTRON", 0., 0., 0., -1., -.5);
    leptons[MU] = Particle("MU", 0., 0., 0., -1., -.5);
    leptons[TAU] = Particle("TAU", 0., 0., 0., -1., -.5);

    ModelParamMap.insert(std::make_pair("Mz", std::cref(Mz)));
    ModelParamMap.insert(std::make_pair("AlsMz", std::cref(AlsMz)));
    ModelParamMap.insert(std::make_pair("GF", std::cref(GF)));
    ModelParamMap.insert(std::make_pair("ale", std::cref(ale)));
    ModelParamMap.insert(std::make_pair("dAle5Mz", std::cref(dAle5Mz)));
//    ModelParamMap.insert(std::make_pair("Mw_inp", std::cref(Mw_inp)));
    ModelParamMap.insert(std::make_pair("mHl", std::cref(mHl)));
    ModelParamMap.insert(std::make_pair("delMw", std::cref(delMw)));
    ModelParamMap.insert(std::make_pair("delSin2th_l", std::cref(delSin2th_l)));
    ModelParamMap.insert(std::make_pair("delSin2th_q", std::cref(delSin2th_q)));
    ModelParamMap.insert(std::make_pair("delSin2th_b", std::cref(delSin2th_b)));
    ModelParamMap.insert(std::make_pair("delGammaZ", std::cref(delGammaZ)));
    ModelParamMap.insert(std::make_pair("delsigma0H", std::cref(delsigma0H)));
    ModelParamMap.insert(std::make_pair("delR0l", std::cref(delR0l)));
    ModelParamMap.insert(std::make_pair("delR0c", std::cref(delR0c)));
    ModelParamMap.insert(std::make_pair("delR0b", std::cref(delR0b)));
    ModelParamMap.insert(std::make_pair("mneutrino_1", std::cref(leptons[NEUTRINO_1].getMass())));
    ModelParamMap.insert(std::make_pair("mneutrino_2", std::cref(leptons[NEUTRINO_2].getMass())));
    ModelParamMap.insert(std::make_pair("mneutrino_3", std::cref(leptons[NEUTRINO_3].getMass())));
    ModelParamMap.insert(std::make_pair("melectron", std::cref(leptons[ELECTRON].getMass())));
    ModelParamMap.insert(std::make_pair("mmu", std::cref(leptons[MU].getMass())));
    ModelParamMap.insert(std::make_pair("mtau", std::cref(leptons[TAU].getMass())));
    ModelParamMap.insert(std::make_pair("lambda", std::cref(lambda)));
    ModelParamMap.insert(std::make_pair("A", std::cref(A)));
    ModelParamMap.insert(std::make_pair("rhob", std::cref(rhob)));
    ModelParamMap.insert(std::make_pair("etab", std::cref(etab)));
    ModelParamMap.insert(std::make_pair("muw", std::cref(muw)));
    
    iterationNo = 0;
    realorder = LO;

    w_GSL1 = gsl_integration_workspace_alloc (200);
}

StandardModel::~StandardModel()
{
    if (IsModelInitialized()) {
        if (myEWSMcache != NULL) delete(myEWSMcache);
        if (myOneLoopEW != NULL) delete(myOneLoopEW);
        if (myTwoLoopQCD != NULL) delete(myTwoLoopQCD);
        if (myThreeLoopQCD != NULL) delete(myThreeLoopQCD);
        if (myTwoLoopEW != NULL) delete(myTwoLoopEW);
        if (myThreeLoopEW2QCD != NULL) delete(myThreeLoopEW2QCD);
        if (myThreeLoopEW != NULL) delete(myThreeLoopEW);
        if (myApproximateFormulae != NULL) delete(myApproximateFormulae);
        if (myLeptonFlavour != NULL) delete(myLeptonFlavour);
        /* BEGIN: REMOVE FROM THE PACKAGE */
        if (myTwoFermionsLEP2 != NULL) delete(myTwoFermionsLEP2);
        /* END: REMOVE FROM THE PACKAGE */
    }
}


///////////////////////////////////////////////////////////////////////////
// Initialization

bool StandardModel::InitializeModel()
{
    myEWSMcache = new EWSMcache(*this); ///< A pointer to an object of type EWSMcache.
    myOneLoopEW = new EWSMOneLoopEW(*myEWSMcache); ///< A pointer to an object of type EWSMOneLoopEW.
    myTwoLoopQCD = new EWSMTwoLoopQCD(*myEWSMcache); ///< A pointer to an object of type EWSMTwoLoopQCD.
    myThreeLoopQCD = new EWSMThreeLoopQCD(*myEWSMcache); ///< A pointer to an object of type EWSMThreeLoopQCD.
    myTwoLoopEW = new EWSMTwoLoopEW(*myEWSMcache); ///< A pointer to an object of type EWSMTwoLoopEW.
    myThreeLoopEW2QCD = new EWSMThreeLoopEW2QCD(*myEWSMcache); ///< A pointer to an object of type EWSMThreeLoopEW2QCD.
    myThreeLoopEW = new EWSMThreeLoopEW(*myEWSMcache); ///< A pointer to an object of type EWSMThreeLoopEW.
    myApproximateFormulae = new EWSMApproximateFormulae(*myEWSMcache); ///< A pointer to an object of type EWSMApproximateFormulae.
    myLeptonFlavour = new LeptonFlavour(*this);
    /* BEGIN: REMOVE FROM THE PACKAGE */
    myTwoFermionsLEP2 = new EWSMTwoFermionsLEP2(*myEWSMcache); ///< A pointer to an object of type EWSMTwoFermionsLEP2.
    /* END: REMOVE FROM THE PACKAGE */
    setModelInitialized(true);
    return (true);
}


///////////////////////////////////////////////////////////////////////////
// Parameters

bool StandardModel::Init(const std::map<std::string, double>& DPars)
{
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        if (it->first.compare("AlsM") == 0 || it->first.compare("MAls") == 0)
            throw std::runtime_error("ERROR: inappropriate parameter " + it->first
                + " in model initialization");
        else if (FlagFixMuwMut && it->first.compare("mut") == 0)
            throw std::runtime_error("ERROR: cannot use " + it->first
                + " when FlagFixMuwMut is true: use only muw");

    std::map<std::string, double> myDPars(DPars);
    myDPars["AlsM"] = myDPars.at("AlsMz"); // do not change!
    myDPars["MAls"] = myDPars.at("Mz");
    if (FlagFixMuwMut)
        myDPars["mut"] = myDPars.at("muw") * 163. / 80.4 ;
    return (QCD::Init(myDPars));
}

bool StandardModel::PreUpdate()
{
    requireCKM = false;
    requireYe = false;
    requireYn = false;

    if (!QCD::PreUpdate()) return (false);

    return (true);
}

bool StandardModel::Update(const std::map<std::string, double>& DPars)
{
    if (!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if (!PostUpdate()) return (false);

    return (true);
}

bool StandardModel::PostUpdate()
{
    if (!QCD::PostUpdate()) return (false);

    /* Set the CKM and PMNS matrices */
    computeCKM();
    
    /* Compute the 5-quark contribution to the running of alpha*/
    dAl5hMz = Dalpha5hMz();   

    /* Set the Yukawa matrices */
    if (!isModelSUSY()) {
        computeYukawas();
    }

    /* Check whether the parameters for the EWPO are updated or not */
    if (!checkSMparamsForEWPO()) {
        useDeltaAlphaLepton_cache = false;
        useDeltaAlpha_cache = false;
        useMw_cache = false;
        useGammaW_cache = false;
        for (int i = 0; i < 12; ++i) {
            useRhoZ_f_cache[i] = false;
            useKappaZ_f_cache[i] = false;
        }
    }
    SMFlavour.setSMupdated();
    /* Necessary for updating StandardModel parameters in StandardModelMatching */
    if (!isModelSUSY()) SMM.getObj().updateSMParameters();

    iterationNo++;

    return (true);
}

void StandardModel::setParameter(const std::string name, const double& value)
{
    if (name.compare("Mz") == 0) {
        Mz = value;
        QCD::setParameter("MAls", value);
    } else if (name.compare("AlsMz") == 0) {
        AlsMz = value;
        QCD::setParameter("AlsM", value);
    } else if (name.compare("GF") == 0)
        GF = value;
    else if (name.compare("ale") == 0)
        ale = value;
    else if (name.compare("dAle5Mz") == 0 && !FlagMWinput)
        dAle5Mz = value;
    else if (name.compare("Mw_inp") == 0 && FlagMWinput)
        Mw_inp = value;
    else if (name.compare("mHl") == 0)
        mHl = value;
    else if (name.compare("delMw") == 0)
        delMw = value;
    else if (name.compare("delSin2th_l") == 0)
        delSin2th_l = value;
    else if (name.compare("delSin2th_q") == 0)
        delSin2th_q = value;
    else if (name.compare("delSin2th_b") == 0)
        delSin2th_b = value;
    else if (name.compare("delGammaZ") == 0)
        delGammaZ = value;
    else if (name.compare("delsigma0H") == 0)
        delsigma0H = value;
    else if (name.compare("delR0l") == 0)
        delR0l = value;
    else if (name.compare("delR0c") == 0)
        delR0c = value;
    else if (name.compare("delR0b") == 0)
        delR0b = value;
    else if (name.compare("mneutrino_1") == 0) 
        leptons[NEUTRINO_1].setMass(value);
    else if (name.compare("mneutrino_2") == 0) 
        leptons[NEUTRINO_2].setMass(value);
    else if (name.compare("mneutrino_3") == 0) 
        leptons[NEUTRINO_3].setMass(value);
    else if (name.compare("melectron") == 0) 
        leptons[ELECTRON].setMass(value);
    else if (name.compare("mmu") == 0) 
        leptons[MU].setMass(value);
    else if (name.compare("mtau") == 0) 
        leptons[TAU].setMass(value);
    else if (name.compare("lambda") == 0 && FlagWolfenstein) {
        lambda = value;
        requireCKM = true;
    } else if (name.compare("A") == 0 && FlagWolfenstein) {
        A = value;
        requireCKM = true;
    } else if (name.compare("rhob") == 0 && FlagWolfenstein) {
        rhob = value;
        requireCKM = true;
    } else if (name.compare("etab") == 0 && FlagWolfenstein) {
        etab = value;
        requireCKM = true;
    } else if (name.compare("V_us") == 0 && !FlagWolfenstein && !FlagUseVud) {
        Vus = value;
        requireCKM = true;
    } else if (name.compare("V_ud") == 0 && !FlagWolfenstein && FlagUseVud) {
        Vud = value;
        requireCKM = true;
    } else if (name.compare("V_cb") == 0 && !FlagWolfenstein) {
        Vcb = value;
        requireCKM = true;
    } else if (name.compare("V_ub") == 0 && !FlagWolfenstein) {
        Vub = value;
        requireCKM = true;
    } else if (name.compare("gamma") == 0 && !FlagWolfenstein) {
        gamma = value;
        requireCKM = true;
    } else if (name.compare("muw") == 0) {
            /* Update mut if FlagFixMuwMut is activated */
        muw = value;            
        if (FlagFixMuwMut)  {
            mut = muw / 80.4 * 163.;
        }
    else
        QCD::setParameter(name, value);
}

bool StandardModel::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NSMvars; i++) {
        if (DPars.find(SMvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory SM parameter " << SMvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(SMvars[i]);
        }
    }
    return (QCD::CheckParameters(DPars));
}

void StandardModel::computeCKM()
{
    if (requireCKM) {
        if (FlagWolfenstein) {
            myCKM.computeCKMwithWolfenstein(lambda, A, rhob, etab);
            Vus = myCKM.getV_us().abs();
            Vcb = myCKM.getV_cb().abs();
            Vub = myCKM.getV_ub().abs();
            gamma = myCKM.computeGamma();
        } else if (FlagUseVud) {
            myCKM.computeCKM(Vud, Vcb, Vub, gamma, FlagUseVud);
            lambda = myCKM.getLambda();
            A = myCKM.getA();
            rhob = myCKM.getRhoBar();
            etab = myCKM.getEtaBar();
            Vus = myCKM.getV_us().abs();
        } else {
            myCKM.computeCKM(Vus, Vcb, Vub, gamma);
            lambda = myCKM.getLambda();
            A = myCKM.getA();
            rhob = myCKM.getRhoBar();
            etab = myCKM.getEtaBar();
            Vud = myCKM.getV_ud().abs();
        }
    }
    myPMNS.computePMNS(s12, s13, s23, delta, alpha21, alpha31); // WARNING: This does not do anything since the input values are not set.
}

///////////////////////////////////////////////////////////////////////////
// Flagsvoid StandardModel::computeYukawas()
{
    if (requireYu || requireCKM) {
        Yu.reset();
        for (int i = 0; i < 3; i++) {
            Yu.assign(i, i, this->getmq(quark(UP + 2 * i),  v()/ sqrt(2.))/ v() * sqrt(2.));
//            std::cout << quarks[UP + 2 * i].getName() << " mass at EW scale is " << this->getmq(quark(UP + 2 * i),  v() / sqrt(2.)) << std::endl;
        }
//        std::cout << "(top MSbar mass is " << this->Mp2Mbar(this->getMtpole()) << ")" << std::endl;
        Yu = Yu * myCKM.getCKM();
    }
    if (requireYd) {
        Yd.reset();
        for (int i = 0; i < 3; i++) {
            Yd.assign(i, i, this->getmq(quark(DOWN + 2 * i),  v() / sqrt(2.)) / v() * sqrt(2.));
//            std::cout << quarks[DOWN + 2 * i].getName() << " mass at " << v() / sqrt(2) << " is " << this->getmq(quark(DOWN + 2 * i),  v() / sqrt(2.)) << std::endl;            
            }
    }
    if (requireYe) {
        Ye = gslpp::matrix<gslpp::complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / v() * sqrt(2.));
    }
    if (requireYn) {
        Yn = gslpp::matrix<gslpp::complex>::Id(3);
        for (int i = 0; i < 3; i++)
            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / v() * sqrt(2.));
        Yn = Yn * myPMNS.getPMNS().hconjugate();
    }
}



bool StandardModel::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if (name.compare("CacheInStandardModel") == 0) {
        setFlagCacheInStandardModel(value);
        res = true;
    } else if (name.compare("CacheInEWSMcache") == 0) {
        getMyEWSMcache()->setFlagCacheInEWSMcache(value);
        res = true;
    } else if (name.compare("Wolfenstein") == 0) {
        FlagWolfenstein = value;
        if(!FlagWolfenstein) {
            SMvars[std::distance(SMvars,std::find(SMvars,SMvars+NSMvars,"lambda"))] = "V_us";
            SMvars[std::distance(SMvars,std::find(SMvars,SMvars+NSMvars,"A"))] = "V_cb";
            SMvars[std::distance(SMvars,std::find(SMvars,SMvars+NSMvars,"rhob"))] = "V_ub";
            SMvars[std::distance(SMvars,std::find(SMvars,SMvars+NSMvars,"etab"))] = "gamma";
            
            ModelParamMap.insert(std::make_pair("V_us", std::cref(Vus)));
            ModelParamMap.insert(std::make_pair("V_cb", std::cref(Vcb)));
            ModelParamMap.insert(std::make_pair("V_ub", std::cref(Vub)));
            ModelParamMap.insert(std::make_pair("gamma", std::cref(gamma)));
        }
        res = true;
    } else if (name.compare("WithoutNonUniversalVC") == 0) {
        FlagWithoutNonUniversalVC = value;
        res = true;
    } else if (name.compare("NoApproximateGammaZ") == 0) {
        FlagNoApproximateGammaZ = value;
        res = true;
    } else if (name.compare("MWinput") == 0) {
        FlagMWinput = value;
        if (FlagMWinput) {
            SMvars[std::distance(SMvars,std::find(SMvars,SMvars+NSMvars,"dAle5Mz"))] = "Mw_inp";
            ModelParamMap.insert(std::make_pair("Mw_inp", std::cref(Mw_inp)));
            // Point the different flags towards the approximate formulae, when available
            FlagNoApproximateGammaZ = false;
            FlagMw = "APPROXIMATEFORMULA";
            FlagRhoZ = "NORESUM";
            FlagKappaZ = "APPROXIMATEFORMULA";
        }
        res = true;
    } else if (name.compare("SMAux") == 0) {
        FlagSMAux = value;
        res = true;
    } else if (name.compare("FixMuwMut") == 0) {
        FlagFixMuwMut = value;
        res = true;
    } else if (name.compare("UseVud") == 0) {
        FlagUseVud = value;
        if (FlagUseVud && FlagWolfenstein)
            throw std::runtime_error("UseVud can only be used when Wolfenstein is false");
        else if(FlagUseVud) {
            SMvars[std::distance(SMvars,std::find(SMvars,SMvars+NSMvars,"V_us"))] = "V_ud";
            ModelParamMap.erase("V_us");
            ModelParamMap.insert(std::make_pair("V_ud", std::cref(Vud)));
        }
        res = true;
    } else
        res = QCD::setFlag(name, value);
    
    if (!res) res = SMFlavour.setFlag(name, value);

    return (res);
}

bool StandardModel::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if (name.compare("Mw") == 0) {
        if (checkEWPOscheme(value)) {
            FlagMw = value;
            res = true;
        } else
            throw std::runtime_error("StandardModel::setFlagStr(): Invalid flag "
                + name + "=" + value);

    } else if (name.compare("RhoZ") == 0) {
        if (checkEWPOscheme(value)) {
            FlagRhoZ = value;
            res = true;
        } else
            throw std::runtime_error("StandardModel::setFlagStr(): Invalid flag "
                + name + "=" + value);
    } else if (name.compare("KappaZ") == 0) {
        if (checkEWPOscheme(value)) {
            FlagKappaZ = value;
            res = true;
        } else
            throw std::runtime_error("StandardModel::setFlagStr(): Invalid flag "
                + name + "=" + value);
    } else
        res = QCD::setFlagStr(name, value);
    
    if (FlagMWinput) {
        // Point the different flags towards the approximate formulae, when available
        FlagNoApproximateGammaZ = false;
        FlagMw = "APPROXIMATEFORMULA";
        FlagRhoZ = "NORESUM";
        FlagKappaZ = "APPROXIMATEFORMULA";
    }

    return (res);
}

bool StandardModel::CheckFlags() const
{
    return (QCD::CheckFlags());
}


///////////////////////////////////////////////////////////////////////////
// For EWPO caches

bool StandardModel::checkSMparamsForEWPO()
{
    // 11 parameters in QCD:
    // AlsMz, Mz, mup, mdown, mcharm, mstrange, mtop, mbottom,
    // mut, mub, muc
    // 19 parameters in StandardModel
    // GF, ale, dAle5Mz, mHl,
    // mneutrino_1, mneutrino_2, mneutrino_3, melectron, mmu, mtau,
    // delMw, delSin2th_l, delSin2th_q, delSin2th_b, delGammaZ, delsigma0H, delR0l, delR0c, delR0b,
    // 3 flags in StandardModel
    // FlagMw_cache, FlagRhoZ_cache, FlagKappaZ_cache

    // Note: When modifying the array below, the constant NumSMParams has to
    // be modified accordingly.
    double SMparams[NumSMParamsForEWPO] = {
        AlsMz, Mz, GF, ale, FlagMWinput? Mw_inp: dAle5Mz,
        mHl, mtpole,
        leptons[NEUTRINO_1].getMass(),
        leptons[NEUTRINO_2].getMass(),
        leptons[NEUTRINO_3].getMass(),
        leptons[ELECTRON].getMass(),
        leptons[MU].getMass(),
        leptons[TAU].getMass(),
        quarks[UP].getMass(),
        quarks[DOWN].getMass(),
        quarks[CHARM].getMass(),
        quarks[STRANGE].getMass(),
        quarks[BOTTOM].getMass(),
        mut, mub, muc,
        delMw, delSin2th_l, delSin2th_q, delSin2th_b, delGammaZ, delsigma0H, delR0l, delR0c, delR0b,
        SchemeToDouble(FlagMw),
        SchemeToDouble(FlagRhoZ),
        SchemeToDouble(FlagKappaZ)
    };

    // check updated parameters
    bool bNotUpdated = true;
    for (int i = 0; i < NumSMParamsForEWPO; ++i) {
        if (SMparamsForEWPO_cache[i] != SMparams[i]) {
            SMparamsForEWPO_cache[i] = SMparams[i];
            bNotUpdated &= false;
        }
    }

    return bNotUpdated;
}

////////////////////////////////////////////////////////////////////////

double StandardModel::ale_OS(const double mu, orders order) const
{
    if (mu < 50.0)
        throw std::runtime_error("out of range in StandardModel::ale_OS()");

    double N = 20.0 / 3.0;
    double beta1 = N / 3.0;
    double beta2 = N / 4.0;
    double alpha_ini = alphaMz();
    double v = 1.0 + 2.0 * beta1 * alpha_ini / M_PI * log(Mz / mu);

    switch (order) {
        case LO:
            return ( alpha_ini / v);
        case FULLNLO:
            return ( alpha_ini / v * (1.0 - beta2 / beta1 * alpha_ini / M_PI * log(v) / v));
        default:
            throw std::runtime_error("Error in StandardModel::ale_OS()");
    }
}

double StandardModel::Beta_s(int nm, unsigned int nf) const
{
    unsigned int nu = nf % 2 == 0 ? nf / 2 : nf / 2;
    unsigned int nd = nf % 2 == 0 ? nf / 2 : 1 + nf / 2;
    double Qu = 2. / 3., Qd = -1. / 3., Qbar2 = nu * Qu * Qu + nd * Qd * Qd,
            Qbar4 = nu * Qu * Qu * Qu * Qu + nd * Qd * Qd * Qd * Qd;

    switch(nm)
    {
        case 00:
            return(Beta0((double) nf));
        case 10:
            return(Beta1((double) nf));
        case 20:
            return(Beta2((double) nf));
        case 30:
            return(Beta3((double) nf));
        case 01:
            return(-4. * TF * Qbar2 );
        case 11:
            return((4. * CF - 8. * CA) * TF * Qbar2 );
        case 02:
            return(11./3. * TF * Qbar2 * Beta_e(00, nf) + 2. * TF * Qbar4);
        default:
            throw std::runtime_error("StandardModel::Beta_s(): case not implemented");
    }
}

double StandardModel::Beta_e(int nm, unsigned int nf) const
{
    unsigned int nu = nf % 2 == 0 ? nf / 2 : nf / 2;
    unsigned int nd = nf % 2 == 0 ? nf / 2 : 1 + nf / 2;
    double Qu = 2. / 3., Qd = -1. / 3., Qbar2 = nu * Qu * Qu + nd * Qd * Qd,
            Qbar4 = nu * Qu * Qu * Qu * Qu + nd * Qd * Qd * Qd * Qd;

    switch(nm)
    {
        case 00:
            return(4./3. * (Qbar2 * Nc + 3.)); // QL^2 = 1
        case 10:
            return(4. * (Qbar4 * Nc + 3.));
        case 01:
            return(4. * CF * Nc * Qbar2);
        default:
            throw std::runtime_error("StandardModel::Beta_e(): case not implemented");
    }
}

double StandardModel::Als(double mu, orders order, bool qed_flag, bool Nf_thr) const
{
    switch (order)
    {
        case LO:
            realorder = order;
            return AlsByOrder(mu, LO, qed_flag, Nf_thr);
        case FULLNLO:
            realorder = order;
            return (AlsByOrder(mu, LO, qed_flag, Nf_thr) + AlsByOrder(mu, NLO, qed_flag, Nf_thr));
        case FULLNNLO:
            realorder = order;
            return (AlsByOrder(mu, LO, qed_flag, Nf_thr) + AlsByOrder(mu, NLO, qed_flag, Nf_thr) + AlsByOrder(mu, NNLO, qed_flag, Nf_thr));
        case FULLNNNLO:
            realorder = order;
            return (AlsByOrder(mu, LO, qed_flag, Nf_thr) + AlsByOrder(mu, NLO, qed_flag, Nf_thr) + AlsByOrder(mu, NNLO, qed_flag, Nf_thr) + AlsByOrder(mu, NNNLO, qed_flag, Nf_thr));
        default:
            throw std::runtime_error("StandardModel::Als(): " + orderToString(order) + " is not implemented.");
    }
}

double StandardModel::AlsByOrder(double mu, orders order, bool qed_flag, bool Nf_thr) const
{
    int i, nfAls = (int) Nf(Mz), nfmu = Nf_thr ? (int) Nf(mu) : nfAls;
    double als, alstmp, mutmp;
    orders fullord;
    
    for (i = 0; i < CacheSize; ++i)
        if ((mu == als_cache[0][i]) && ((double) order == als_cache[1][i]) &&
                (AlsMz == als_cache[2][i]) && (Mz == als_cache[3][i]) &&
                (mut == als_cache[4][i]) && (mub == als_cache[5][i]) &&
                (muc == als_cache[6][i]) && (double) qed_flag == als_cache[7][i]
                && (double) Nf_thr == als_cache[8][i] && alphaMz() == als_cache[9][i])
            return als_cache[10][i];

    switch (order)
    {
        case LO:
        case NLO:
        case NNLO:
        case NNNLO:
            if (nfAls == nfmu) 
                als = AlsWithInit(mu, AlsMz, Mz, order, qed_flag);
            fullord = FullOrder(order);
            if (nfAls > nfmu) {
                mutmp = BelowTh(Mz);
                alstmp = AlsWithInit(mutmp, AlsMz, Mz, realorder, qed_flag);
                alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(nfAls), alstmp, nfAls, fullord)); // WARNING: QED threshold corrections not implemented yet
                for (i = nfAls - 1; i > nfmu; i--) {
                    mutmp = BelowTh(mutmp - MEPS);
                    alstmp = AlsWithInit(mutmp, alstmp, AboveTh(mutmp) - MEPS, realorder, qed_flag);
                    alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(i), alstmp, i, fullord)); // WARNING: QED threshold corrections not implemented yet
                }
                als = AlsWithInit(mu, alstmp, AboveTh(mu) - MEPS, order, qed_flag);
            }

            if (nfAls < nfmu) {
                mutmp = AboveTh(Mz) - MEPS;
                alstmp = AlsWithInit(mutmp, AlsMz, Mz, realorder, qed_flag);
                alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(nfAls + 1), alstmp, nfAls + 1, fullord)); // WARNING: QED threshold corrections not implemented yet
                for (i = nfAls + 1; i < nfmu; i++) {
                    mutmp = AboveTh(mutmp) - MEPS;
                    alstmp = AlsWithInit(mutmp, alstmp, BelowTh(mutmp) + MEPS, realorder, qed_flag); 
                    alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(i + 1), alstmp, i + 1, fullord)); // WARNING: QED threshold corrections not implemented yet
                }
                als = AlsWithInit(mu, alstmp, BelowTh(mu) + MEPS, order, qed_flag);
            }

            CacheShift(als_cache, 11);
            als_cache[0][0] = mu;
            als_cache[1][0] = (double) order;       
            als_cache[2][0] = AlsMz;
            als_cache[3][0] = Mz;
            als_cache[4][0] = mut;
            als_cache[5][0] = mub;
            als_cache[6][0] = muc;
            als_cache[7][0] = (double) qed_flag;
            als_cache[8][0] = (double) Nf_thr;
            als_cache[9][0] = alphaMz();
            als_cache[10][0] = als;

             return als;
        default:
            throw std::runtime_error("StandardModel::Als(): " + orderToString(order) + " is not implemented.");
    }
}

double StandardModel::AlsWithInit(double mu, double alsi, double mu_i, orders order, bool qed_flag) const
{
    double nf = Nf(mu), alei = Ale(mu_i, FULLNLO); // CHANGE ME!
    double b00s = Beta_s(00, nf), b00e = Beta_e(00, nf);
    double v = 1. + b00s * alsi / 2. / M_PI * log(mu / mu_i);
    double ve = 1. - b00e * alei / 2. / M_PI * log(mu / mu_i);
    double logv = log(v), logve = log(ve);
    double rho = 1. / (1. + b00e * alei / b00s / alsi);
    double als = QCD::AlsWithInit(mu, alsi, mu_i, order);
    double b01s = Beta_s(01,nf), b01s00e = b01s / b00e;

    if (qed_flag)
        switch (order)
        {
            case LO:
                break;
            case NLO:
                als += alsi * alsi / 4. / M_PI / v / v * b01s00e * logve;
                break;
            case NNLO:
                als += alsi * alsi * alsi / 4. / 4. / M_PI / M_PI / v / v / v * (
                        b01s00e * b01s00e * logve * logve + b01s00e * Beta_s(10, nf) / b00s *
                        (-2. * logv * logve + rho * ve * logve));
                break;
            case NNNLO:
                als += alsi * alsi * alei / 4. / 4. / M_PI / M_PI / v / v / ve * (Beta_s(02, nf) / b00e *
                        (ve - 1.) + Beta_s(11, nf) / b00s * rho * ve * (logve - logv) + b01s00e * Beta_e(10, nf) /
                        b00e * (logve - ve + 1.) + b01s * Beta_s(10, nf) / b00s / b00s * rho * logv +
                        b01s00e * Beta_e(01, nf) / b00s * (rho * ve * (logv - logve) - logv));
                break;
            case FULLNLO:
                return (AlsWithInit(mu, alsi, mu_i, LO, true) + AlsWithInit(mu, alsi, mu_i, NLO, true));
            case FULLNNLO:
                return (AlsWithInit(mu, alsi, mu_i, LO, true) + AlsWithInit(mu, alsi, mu_i, NLO, true)+ AlsWithInit(mu, alsi, mu_i, NNLO, true));
            case FULLNNNLO:
                return (AlsWithInit(mu, alsi, mu_i, LO, true) + AlsWithInit(mu, alsi, mu_i, NLO, true)+ AlsWithInit(mu, alsi, mu_i, NNLO, true) + AlsWithInit(mu, alsi, mu_i, NNNLO, true));
            default:
                throw std::runtime_error("StandardModel::AlsWithInit(): " + orderToString(order) + " is not implemented.");
        }

    return (als);
}

double StandardModel::Ale(const double mu, orders order, bool Nf_thr) const
{
    int i, nfAle = (int) Nf(Mz), nfmu = Nf_thr ? (int) Nf(mu) : nfAle;
    double ale, aletmp, mutmp, aleMz = alphaMz();
    orders fullord;

    for (i = 0; i < CacheSize; ++i)
        if ((mu == ale_cache[0][i]) && ((double) order == ale_cache[1][i]) &&
                (AlsMz == ale_cache[2][i]) && (Mz == ale_cache[3][i]) &&
                (mut == ale_cache[4][i]) && (mub == ale_cache[5][i]) &&
                (muc == ale_cache[6][i])
                && (double) Nf_thr == ale_cache[7][i] && aleMz == ale_cache[8][i])
            return ale_cache[9][i];

    switch (order)
    {
        case FULLNLO:
            return (Ale(mu, LO, Nf_thr) + Ale(mu, NLO, Nf_thr));
        case FULLNNLO:
            return (Ale(mu, LO, Nf_thr) + Ale(mu, NLO, Nf_thr) + Ale(mu, NNLO, Nf_thr));
        case FULLNNNLO:
            return (Ale(mu, LO, Nf_thr) + Ale(mu, NLO, Nf_thr) + Ale(mu, NNLO, Nf_thr) + Ale(mu, NNNLO, Nf_thr));
        case LO:
            if (nfAle == nfmu)
                return(AleWithInit(mu, aleMz, Mz, order));
        case NLO:
        case NNLO:
        case NNNLO:
            if (nfAle == nfmu)
                return(0.);
            fullord = FullOrder(order);
            if (nfAle > nfmu) {
                mutmp = BelowTh(Mz);
                aletmp = AleWithInit(mutmp, aleMz, Mz, fullord);
//                aletmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(nfAle), alstmp, nfAls, fullord)); // WARNING: QED threshold corrections not implemented yet
                for (i = nfAle - 1; i > nfmu; i--) {
                    mutmp = BelowTh(mutmp - MEPS);
                    aletmp = AleWithInit(mutmp, aletmp, AboveTh(mutmp) - MEPS, fullord);
//                    aletmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(i), aletmp, i, fullord)); // WARNING: QED threshold corrections not implemented yet
                }
                ale = AleWithInit(mu, aletmp, AboveTh(mu) - MEPS, order);
            }

            if (nfAle < nfmu) {
                mutmp = AboveTh(Mz) - MEPS;
                aletmp = AleWithInit(mutmp, aleMz, Mz, fullord);
//                alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(nfAls + 1), alstmp, nfAls + 1, fullord)); // WARNING: QED threshold corrections not implemented yet
                for (i = nfAle + 1; i < nfmu; i++) {
                    mutmp = AboveTh(mutmp) - MEPS;
                    aletmp = AleWithInit(mutmp, aletmp, BelowTh(mutmp) + MEPS, fullord); 
//                    alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(i + 1), alstmp, i + 1, fullord)); // WARNING: QED threshold corrections not implemented yet
                }
                ale = AleWithInit(mu, aletmp, BelowTh(mu) + MEPS, order);
            }

            CacheShift(ale_cache, 10);
            ale_cache[0][0] = mu;
            ale_cache[1][0] = (double) order;    
            ale_cache[2][0] = AlsMz;
            ale_cache[3][0] = Mz;
            ale_cache[4][0] = mut;
            ale_cache[5][0] = mub;
            ale_cache[6][0] = muc;
            ale_cache[7][0] = (double) Nf_thr;
            ale_cache[8][0] = aleMz;
            ale_cache[9][0] = ale;

            return ale;
        default:
            throw std::runtime_error("StandardModel::Ale(): " + orderToString(order) + " is not implemented.");
    }
}

double StandardModel::AleWithInit(double mu, double alei, double mu_i, orders order) const
{
    if (fabs(mu - mu_i) < MEPS) return(alei);

    double nf = Nf(mu), alsi = (mu_i == Mz ? AlsMz : Als(mu_i, FULLNNNLO, true));
    double b00e = Beta_e(00, nf), b00s = Beta_s(00, nf);
    double ve = 1. - b00e * alei / 2. / M_PI * log(mu / mu_i);
    double logv = log(1. + b00s * alsi / 2. / M_PI * log(mu / mu_i)), logve = log(ve);

    switch (order)
    {
        case LO:
            return (alei / ve);
        case NLO:
            return (- alei * alei / 4. / M_PI / ve / ve * (Beta_e(10, nf) / b00e * logve - Beta_e(01, nf) / b00s * logv) );
            // Higher order terms ? Need to understand eq. (35)
        case FULLNLO:
            return (AleWithInit(mu, alei, mu_i, LO) + AleWithInit(mu, alei, mu_i, NLO));
        default:
            throw std::runtime_error("StandardModel::AleWithInit(): " + orderToString(order) + " is not implemented.");
    }
}

double StandardModel::DeltaAlphaLepton(const double s) const
{
    if (s == Mz * Mz)
        if (FlagCacheInStandardModel)
            if (useDeltaAlphaLepton_cache)
                return DeltaAlphaLepton_cache;

    double DeltaAlphaL = 0.0;
    if (flag_order[EW1])
        DeltaAlphaL += myOneLoopEW->DeltaAlpha_l(s);
    if (flag_order[EW1QCD1])
        DeltaAlphaL += myTwoLoopQCD->DeltaAlpha_l(s);
    if (flag_order[EW1QCD2])
        DeltaAlphaL += myThreeLoopQCD->DeltaAlpha_l(s);
    if (flag_order[EW2])
        DeltaAlphaL += myTwoLoopEW->DeltaAlpha_l(s);
    if (flag_order[EW2QCD1])
        DeltaAlphaL += myThreeLoopEW2QCD->DeltaAlpha_l(s);
    if (flag_order[EW3])
        DeltaAlphaL += myThreeLoopEW->DeltaAlpha_l(s);

    if (s == Mz * Mz) {
        DeltaAlphaLepton_cache = DeltaAlphaL;
        useDeltaAlphaLepton_cache = true;
    }
    return DeltaAlphaL;
}

double StandardModel::DeltaAlphaL5q() const
{
    double Mz2 = Mz*Mz;
    return (DeltaAlphaLepton(Mz2) + dAl5hMz);
}

double StandardModel::DeltaAlphaTop(const double s) const
{
    double DeltaAlpha = 0.0;
    if (flag_order[EW1])
        DeltaAlpha += myOneLoopEW->DeltaAlpha_t(s);
    if (flag_order[EW1QCD1])
        DeltaAlpha += myTwoLoopQCD->DeltaAlpha_t(s);
    if (flag_order[EW1QCD2])
        DeltaAlpha += myThreeLoopQCD->DeltaAlpha_t(s);
    if (flag_order[EW2])
        DeltaAlpha += myTwoLoopEW->DeltaAlpha_t(s);
    if (flag_order[EW2QCD1])
        DeltaAlpha += myThreeLoopEW2QCD->DeltaAlpha_t(s);
    if (flag_order[EW3])
        DeltaAlpha += myThreeLoopEW->DeltaAlpha_t(s);

    return DeltaAlpha;
}

double StandardModel::DeltaAlpha() const
{
    if (FlagCacheInStandardModel)
        if (useDeltaAlpha_cache)
            return DeltaAlpha_cache;

    double Mz2 = Mz*Mz;
    DeltaAlpha_cache = DeltaAlphaL5q() + DeltaAlphaTop(Mz2);
    useDeltaAlpha_cache = true;
    return DeltaAlpha_cache;
}

double StandardModel::alphaMz() const
{
    return (ale / (1.0 - DeltaAlpha()));
//    return(1./127.918); // FOR HEFFDF1 TEST: VALUE IN hep-ph/0512066
//    return(1./127.955); // FOR HEFFDF1 TEST: VALUE IN 2007.04191
}

double StandardModel::Alstilde5(const double mu) const
{
    double mu_0 = Mz;
    double alphatilde_e = alphaMz()/4./M_PI;
    double alphatilde_s = AlsMz/4./M_PI;
    unsigned int nf = 5;

    double B00S = Beta0(nf), B10S = Beta1(nf), B20S = Beta2(nf), B30S = gsl_sf_zeta_int(3) * 352864./81. - 598391./1458,
            B01S = -22./9., B11S = -308./27., B02S = 4945./243.; 

    double B00E = 80./9., B01E = 176./9., B10E = 464./27.; 

    double B10soB00s = B10S / B00S;
    double B01soB00e = B01S/B00E;

    double vs= 1. + 2. * B00S * alphatilde_s * log(mu/ mu_0);
    double ve= 1. - 2. * B00E * alphatilde_e * log(mu/ mu_0);
    double ps= B00S * alphatilde_s /(B00S * alphatilde_s + B00E * alphatilde_e);

    double logve = log(ve);
    double logvs = log(vs);
    double logeos = log(ve/vs);
    double logsoe = log(vs/ve);
    double asovs = alphatilde_s/vs;
    double aeove = alphatilde_e/ve;

    double result = 0;

    result = asovs - pow(asovs, 2) * (logvs * B10soB00s - logve * B01soB00e) 
            +  pow(asovs, 3) * ((1. - vs) * B20S / B00S + B10soB00s * B10soB00s * (logvs * logvs - logvs
            + vs - 1.) + B01soB00e * B01soB00e * logve * logve + (-2. * logvs * logve 
            + ps * ve * logve) * B01S * B10S/(B00E * B00S)) 
            +  pow(asovs, 4) * (0.5 * B30S *(1. - vs * vs)/ B00S + ((2. * vs - 3.) * logvs + vs * vs 
            - vs) * B20S * B10soB00s /(B00S) + B10soB00s * B10soB00s * B10soB00s * (- pow(logvs,3) 
            + 5. * pow(logvs,2) / 2. + 2. * (1. - vs) * logvs - (vs - 1.) * (vs - 1.)* 0.5))
            + pow(asovs, 2) * (aeove) * ((ve - 1.) * B02S / B00E 
            + ps * ve * logeos * B11S /B00S +(logve - ve + 1.) * B01soB00e * B10E/(B00E) 
            + logvs * ps * B01S * B10soB00s/(B00S) +(logsoe * ve * ps - logvs) * B01soB00e * B01E/( B00S));
    return (result);
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::v() const
{
    return ( 1. / sqrt(sqrt(2.) * GF));
}


///////////////////////////////////////////////////////////////////////////

double StandardModel::Mw_tree() const
{
    if (FlagMWinput){
        return Mw_inp;
    } else
        return ( Mz / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0 * M_PI * ale / sqrt(2.0) / GF / Mz / Mz)));
}

double StandardModel::s02() const
{
    double tmp = 1.0 - 4.0 * M_PI * alphaMz() / sqrt(2.0) / GF / Mz / Mz;
    if (tmp < 0.0)
        throw std::runtime_error("Error in s02()");

    return ( (1.0 - sqrt(tmp)) / 2.0);
}

double StandardModel::c02() const
{
    return ( 1.0 - s02());
}

double StandardModel::Mw() const
{
    /* Debug */
    //std::cout << std::boolalpha
    //          << checkScheme(schemeMw_cache,schemeMw,false)
    //          << " [cache:" << schemeMw_cache
    //          << " current:" << schemeMw << "]" << std::endl;

    if (FlagMWinput)
        return Mw_inp;

    if (FlagCacheInStandardModel)
        if (useMw_cache)
            return Mw_cache;

    double Mw;
    if (FlagMw.compare("APPROXIMATEFORMULA") == 0)
        Mw = myApproximateFormulae->Mw();
    else {
        //std::cout << std::setprecision(12)
        //          << "TEST: Mw_tree = " << Mw_tree() << std::endl;

        double DeltaRho[orders_EW_size], DeltaR_rem[orders_EW_size];
        ComputeDeltaRho(Mw_tree(), DeltaRho);
        ComputeDeltaR_rem(Mw_tree(), DeltaR_rem);
        Mw = resumMw(Mw_tree(), DeltaRho, DeltaR_rem);

        /* Mw from iterations */
        double Mw_org = Mw_tree();
        while (fabs(Mw - Mw_org) > Mw_error) {
            Mw_org = Mw;
            ComputeDeltaRho(Mw, DeltaRho);
            ComputeDeltaR_rem(Mw, DeltaR_rem);
            Mw = resumMw(Mw, DeltaRho, DeltaR_rem);
            /* TEST */
            //int prec_def = std::cout.precision();
            //std::cout << std::setprecision(12) << "TEST: Mw_org = " << Mw_org
            //        << "  Mw_new = " << Mw << std::endl;
            //std::cout.precision(prec_def);
        }
    }

//    Mw = 80.426; // FOR HEFFDF1 TEST: VALUE IN hep-ph/0512066
//    Mw = 80.379; // FOR HEFFDF1 TEST: VALUE IN 2007.04191
    Mw_cache = Mw;
    useMw_cache = true;
    return Mw;
}

double StandardModel::Dalpha5hMz() const
{
    if (FlagMWinput){
        return (myApproximateFormulae->dAlpha5hMw());
    } else
        return dAle5Mz;
}

double StandardModel::cW2(double Mw_i) const
{
    return ( Mw_i * Mw_i / Mz / Mz);
}

double StandardModel::cW2() const
{
    return ( cW2(Mw()));
//    return (1.0 - 0.2312); // FOR HEFFDF1 TEST
}

double StandardModel::sW2(double Mw_i) const
{
    return ( 1.0 - cW2(Mw_i));
}

double StandardModel::sW2() const
{
    return ( 1.0 - cW2());
}

double StandardModel::DeltaR() const
{
    /* in the experimental/running-width scheme */
    double myMw = Mw();
    double sW2 = 1.0 - myMw * myMw / Mz / Mz;
    double tmp = sqrt(2.0) * GF * sW2 * myMw * myMw / M_PI / ale;
    if (FlagMw.compare("NORESUM") == 0
            || FlagMw.compare("APPROXIMATEFORMULA") == 0) {
        return (tmp - 1.0);
    } else {
        return (1.0 - 1.0 / tmp);
    }
}

void StandardModel::ComputeDeltaRho(const double Mw_i,
        double DeltaRho[orders_EW_size]) const
{
    if (flag_order[EW1])
        DeltaRho[EW1] = myOneLoopEW->DeltaRho(Mw_i);
    else
        DeltaRho[EW1] = 0.0;
    if (flag_order[EW1QCD1])
        DeltaRho[EW1QCD1] = myTwoLoopQCD->DeltaRho(Mw_i);
    else
        DeltaRho[EW1QCD1] = 0.0;
    if (flag_order[EW1QCD2])
        DeltaRho[EW1QCD2] = myThreeLoopQCD->DeltaRho(Mw_i);
    else
        DeltaRho[EW1QCD2] = 0.0;
    if (flag_order[EW2])
        DeltaRho[EW2] = myTwoLoopEW->DeltaRho(Mw_i);
    else
        DeltaRho[EW2] = 0.0;
    if (flag_order[EW2QCD1])
        DeltaRho[EW2QCD1] = myThreeLoopEW2QCD->DeltaRho(Mw_i);
    else
        DeltaRho[EW2QCD1] = 0.0;
    if (flag_order[EW3])
        DeltaRho[EW3] = myThreeLoopEW->DeltaRho(Mw_i);
    else
        DeltaRho[EW3] = 0.0;
}

void StandardModel::ComputeDeltaR_rem(const double Mw_i,
        double DeltaR_rem[orders_EW_size]) const
{
    if (flag_order[EW1])
        DeltaR_rem[EW1] = myOneLoopEW->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW1] = 0.0;
    if (flag_order[EW1QCD1])
        DeltaR_rem[EW1QCD1] = myTwoLoopQCD->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW1QCD1] = 0.0;
    if (flag_order[EW1QCD2])
        DeltaR_rem[EW1QCD2] = myThreeLoopQCD->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW1QCD2] = 0.0;
    if (flag_order[EW2])
        DeltaR_rem[EW2] = myTwoLoopEW->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW2] = 0.0;
    if (flag_order[EW2QCD1])
        DeltaR_rem[EW2QCD1] = myThreeLoopEW2QCD->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW2QCD1] = 0.0;
    if (flag_order[EW3])
        DeltaR_rem[EW3] = myThreeLoopEW->DeltaR_rem(Mw_i);
    else
        DeltaR_rem[EW3] = 0.0;
}


////////////////////////////////////////////////////////////////////////

double StandardModel::Mzbar() const
{
    double G0 = GF * pow(Mz, 3.0) / 24.0 / sqrt(2.0) / M_PI;
    double sW2tree = 1.0 - Mw_tree() * Mw_tree() / Mz / Mz;
    double Gz = 6.0 * G0; // neutrinos
    Gz += 3.0 * G0 * (pow(1.0 - 4.0 * sW2tree, 2.0) + 1.0); // e, mu and tau
    Gz += 6.0 * G0 * (pow(1.0 - 8.0 / 3.0 * sW2tree, 2.0) + 1.0)
            * (1.0 + AlsMz / M_PI); // u and c
    Gz += 9.0 * G0 * (pow(1.0 - 4.0 / 3.0 * sW2tree, 2.0) + 1.0)
            * (1.0 + AlsMz / M_PI); // d, s and b

    //Gz = 2.4952; // experimental data
    //std::cout << "Gz=" << Gz << std::endl; // for test

    return ( Mz - Gz * Gz / 2.0 / Mz);
}

double StandardModel::MwbarFromMw(const double Mw) const
{
    double AlsMw = Als(Mw, FULLNLO);
    double Gw_SM = 3.0 * GF * pow(Mw, 3.0) / 2.0 / sqrt(2.0) / M_PI
            * (1.0 + 2.0 * AlsMw / 3.0 / M_PI);

    return ( Mw - Gw_SM * Gw_SM / 2.0 / Mw);
}

double StandardModel::MwFromMwbar(const double Mwbar) const
{
    double AlsMw = Als(Mwbar, FULLNNLO);
    double Gw_SM = 3.0 * GF * pow(Mwbar, 3.0) / 2.0 / sqrt(2.0) / M_PI
            * (1.0 + 2.0 * AlsMw / 3.0 / M_PI);

    return (Mwbar + Gw_SM * Gw_SM / 2.0 / Mwbar);
}

double StandardModel::DeltaRbar() const
{
    double Mwbar_SM = MwbarFromMw(Mw());
    double sW2bar = 1.0 - Mwbar_SM * Mwbar_SM / Mzbar() / Mzbar();
    double tmp = sqrt(2.0) * GF * sW2bar * Mwbar_SM * Mwbar_SM / M_PI / ale;

    return (tmp - 1.0);
}


////////////////////////////////////////////////////////////////////////

double StandardModel::rho_GammaW(const Particle fi, const Particle fj) const
{
    double rhoW = 0.0;
    if (flag_order[EW1])
        rhoW = myOneLoopEW->rho_GammaW(fi, fj, Mw());
    return rhoW;
}

double StandardModel::GammaW(const Particle fi, const Particle fj) const
{
    if ((fi.getIndex()) % 2 || (fj.getIndex() + 1) % 2)
        throw std::runtime_error("Error in StandardModel::GammaW()");

    double G0 = GF * pow(Mw(), 3.0) / 6.0 / sqrt(2.0) / M_PI;
    gslpp::complex V(0.0, 0.0, false);

    if (fi.is("TOP"))
        return (0.0);

    if (fj.getIndex() - fi.getIndex() == 1)
        V = gslpp::complex(1.0, 0.0, false);
    else
        V = gslpp::complex(0.0, 0.0, false);

    if (fi.is("LEPTON"))
        return ( V.abs2() * G0 * rho_GammaW(fi, fj));
    else {
        double AlsMw = AlsWithInit(Mw(), AlsMz, Mz, FULLNLO, false);
        return ( 3.0 * V.abs2() * G0 * rho_GammaW(fi, fj)*(1.0 + AlsMw / M_PI));
    }
}

double StandardModel::GammaW() const
{
    if (FlagCacheInStandardModel)
        if (useGammaW_cache)
            return GammaW_cache;

    double GammaWtmp = 0.;

    for (int i = 0; i < 6; i += 2)
        GammaWtmp += GammaW(leptons[i], leptons[i + 1]) + GammaW(quarks[i], quarks[i + 1]);

    GammaW_cache = GammaWtmp;
    useGammaW_cache = true;
    return GammaWtmp;
}


double StandardModel::BrW(const Particle fi, const Particle fj) const
{
    double GammW = GammaW();
    double GammWij = GammaW(fi, fj);

    return GammWij/GammW;
}


double StandardModel::RWlilj(const Particle li, const Particle lj) const
{
    double GammWli, GammWlj;
    
    if (li.is("ELECTRON"))
        GammWli = GammaW(leptons[NEUTRINO_1],li);
    else if (li.is("MU"))
        GammWli = GammaW(leptons[NEUTRINO_2],li);        
    else if (li.is("TAU"))
        GammWli = GammaW(leptons[NEUTRINO_3],li);        
    else
        throw std::runtime_error("Error in StandardModel::RWlilj. li must be a charged lepton");
    
    if (lj.is("ELECTRON"))
        GammWlj = GammaW(leptons[NEUTRINO_1],lj);
    else if (lj.is("MU"))
        GammWlj = GammaW(leptons[NEUTRINO_2],lj);        
    else if (lj.is("TAU"))
        GammWlj = GammaW(leptons[NEUTRINO_3],lj);        
    else
        throw std::runtime_error("Error in StandardModel::RWlilj. lj must be a charged lepton");
    
    return GammWli/GammWlj;
}

double StandardModel::Ruc() const       //AG:added
{
    return 0.5 * ( R0_f(quarks[UP]) + R0_f(quarks[CHARM]) );
}

double StandardModel::RWc() const
{  
    double GammWcX, GammWhad;

//  Add all the  W-> cX decays
//  In GammaW fermion masses are ignored and CKM=1 but uses that SM CKM is unitary => I only need W->cs
    GammWcX = GammaW(quarks[CHARM], quarks[STRANGE]);
    
//  For the same reasons, I only need to add the W-> ud decays into the hadronic part
    GammWhad = GammWcX
            + GammaW(quarks[UP], quarks[DOWN]);

    return GammWcX/GammWhad;
}

////////////////////////////////////////////////////////////////////////

double StandardModel::A_f(const Particle f) const
{
    double Re_kappa = kappaZ_f(f).real();
    double Re_gV_over_gA = 1.0 - 4.0 * fabs(f.getCharge()) * Re_kappa * sW2();
    return ( 2.0 * Re_gV_over_gA / (1.0 + pow(Re_gV_over_gA, 2.0)));
}

double StandardModel::AFB(const Particle f) const
{
    return (3.0 / 4.0 * A_f(leptons[ELECTRON]) * A_f(f));
}

double StandardModel::sin2thetaEff(const Particle f) const
{
    double Re_kappa = kappaZ_f(f).real();
    return ( Re_kappa * sW2());
}

double StandardModel::GammaZ(const Particle f) const
{
    if (f.is("TOP"))
        return 0.0;
    double Gamma;
    if (!IsFlagNoApproximateGammaZ()) {
            
        /* SM contribution with the approximate formula */
        if (f.is("NEUTRINO_1") || f.is("NEUTRINO_2") || f.is("NEUTRINO_3"))
            Gamma = myApproximateFormulae->X_full("Gamma_nu");
        else if (f.is("ELECTRON") || f.is("MU"))
            Gamma = myApproximateFormulae->X_full("Gamma_e_mu");
        else if (f.is("TAU"))
            Gamma = myApproximateFormulae->X_full("Gamma_tau");
        else if (f.is("UP"))
            Gamma = myApproximateFormulae->X_full("Gamma_u");
        else if (f.is("CHARM"))
            Gamma = myApproximateFormulae->X_full("Gamma_c");
        else if (f.is("DOWN") || f.is("STRANGE"))
            Gamma = myApproximateFormulae->X_full("Gamma_d_s");
        else if (f.is("BOTTOM"))
            Gamma = myApproximateFormulae->X_full("Gamma_b");
        else
            throw std::runtime_error("Error in StandardModel::GammaZ()");
        
    } else {
        gslpp::complex myrhoZ_f = rhoZ_f(f);
        gslpp::complex gV_over_gA = gV_f(f) / gA_f(f);
        double G0 = GF * pow(Mz, 3.0) / 24.0 / sqrt(2.0) / M_PI;
        if (f.is("LEPTON")) {
            double myalphaMz = alphaMz();
            double Q = f.getCharge();
            double xl = pow(f.getMass() / Mz, 2.0);
            Gamma = G0 * myrhoZ_f.abs() * sqrt(1.0 - 4.0 * xl)
                    * ((1.0 + 2.0 * xl)*(gV_over_gA.abs2() + 1.0) - 6.0 * xl)
                    * (1.0 + 3.0 / 4.0 * myalphaMz / M_PI * pow(Q, 2.0));
        } else if (f.is("QUARK")) {
            Gamma = 3.0 * G0 * myrhoZ_f.abs()*(gV_over_gA.abs2() * RVq((QCD::quark) (f.getIndex() - 6)) + RAq((QCD::quark) (f.getIndex() - 6)));

            /* Nonfactorizable EW-QCD corrections */
            Gamma += Delta_EWQCD((QCD::quark) (f.getIndex() - 6));
        } else
            throw std::runtime_error("Error in StandardModel::GammaZ()");
    }

    return Gamma;
}

double StandardModel::Gamma_inv() const
{
    return ( GammaZ(leptons[NEUTRINO_1]) + GammaZ(leptons[NEUTRINO_2])
            + GammaZ(leptons[NEUTRINO_3]));
}

double StandardModel::Gamma_had() const
{
    double Gamma_had_tmp = 0.0;
    
    if (!IsFlagNoApproximateGammaZ()){
            
        /* SM contribution with the approximate formula */
        return myApproximateFormulae->X_full("Gamma_had");
    
    } else {
    
        Gamma_had_tmp = GammaZ(quarks[UP]) + GammaZ(quarks[DOWN]) + GammaZ(quarks[CHARM])
            + GammaZ(quarks[STRANGE]) + GammaZ(quarks[BOTTOM]);

    /* Singlet vector contribution (not included in the approximate formula) */
        double G0 = GF * pow(Mz, 3.0) / 24.0 / sqrt(2.0) / M_PI;
        Gamma_had_tmp += 4.0 * 3.0 * G0 * RVh();

        return Gamma_had_tmp;
    }
}

double StandardModel::Gamma_Z() const
{
    if (!IsFlagNoApproximateGammaZ()){
            
        /* SM contribution with the approximate formula */
        return myApproximateFormulae->X_full("GammaZ");

    } else {
        return ( GammaZ(leptons[ELECTRON]) + GammaZ(leptons[MU]) + GammaZ(leptons[TAU])
            + Gamma_inv() + Gamma_had());
    }
}


double StandardModel::RZlilj(const Particle li, const Particle lj) const
{
    double GammZli, GammZlj;
    
    if ( li.is("ELECTRON") || li.is("MU") || li.is("TAU") )
        GammZli = GammaZ(li);        
    else
        throw std::runtime_error("Error in StandardModel::RZlilj. li must be a charged lepton");
    
    if ( lj.is("ELECTRON") || lj.is("MU") || lj.is("TAU") )
        GammZlj = GammaZ(lj);        
    else
        throw std::runtime_error("Error in StandardModel::RZlilj. lj must be a charged lepton");
    
    return GammZli/GammZlj;
}


double StandardModel::sigma0_had() const
{
    if (!IsFlagNoApproximateGammaZ()){
            
        /* SM contribution with the approximate formula */
        return (myApproximateFormulae->X_full("sigmaHadron")
            / GeVminus2_to_nb);

    } else {
        return (12.0 * M_PI * GammaZ(leptons[ELECTRON]) * Gamma_had()
            / Mz / Mz / Gamma_Z() / Gamma_Z());
    }
}

double StandardModel::R0_f(const Particle f) const
{
                
    if (f.is("ELECTRON")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_electron"));
        else
            return (Gamma_had() / GammaZ(leptons[ELECTRON]));
    }  else if (f.is("MU")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_muon"));
        else
            return (Gamma_had() / GammaZ(leptons[MU]));
    }  else if (f.is("TAU")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_tau"));
        else
            return (Gamma_had() / GammaZ(leptons[TAU]));
    } else if (f.is("NEUTRINO_1")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_neutrino"));
        else
            return (GammaZ(leptons[NEUTRINO_1]) / Gamma_had());
    } else if (f.is("NEUTRINO_2")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_neutrino"));
        else
            return (GammaZ(leptons[NEUTRINO_2]) / Gamma_had());
    } else if (f.is("NEUTRINO_3")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_neutrino"));
        else
            return (GammaZ(leptons[NEUTRINO_3]) / Gamma_had());
    }  else if (f.is("UP")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_up"));
        else
            return (GammaZ(quarks[UP]) / Gamma_had());

    }  else if (f.is("STRANGE")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_strange"));
        else
            return (GammaZ(quarks[STRANGE]) / Gamma_had());

    }  else if (f.is("CHARM")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_charm"));
        else
            return (GammaZ(quarks[CHARM]) / Gamma_had());

    } else if (f.is("BOTTOM")) {
        if (!IsFlagNoApproximateGammaZ())
            /* SM contribution with the approximate formula */
            return (myApproximateFormulae->X_full("R0_bottom"));
        else
            return (GammaZ(quarks[BOTTOM]) / Gamma_had());

    } else throw std::runtime_error("StandardModel::R0_f called with wrong argument");     
    
}

double StandardModel::R_inv() const
{
    return (Gamma_inv() / GammaZ(leptons[ELECTRON]));

}

double StandardModel::N_nu() const
{
    double Nnu = 0.0;
    double Gl = 0.0;    
    double Rl = 0.0;
    
    // Don't assume lepton universality: average over lepton flavours
    Gl = GammaZ(leptons[ELECTRON]) + GammaZ(leptons[MU]) + GammaZ(leptons[TAU]);
    Rl = (1.0/3.0) * ( R0_f(leptons[ELECTRON]) + R0_f(leptons[MU]) + R0_f(leptons[TAU]) );
    
    Nnu = sqrt( 12.0 * M_PI * Rl / Mz / Mz / sigma0_had() ) - Rl -3.0;
    
    Nnu = (Gl/Gamma_inv()) * Nnu;
    
    return Nnu;

}


////////////////////////////////////////////////////////////////////////

gslpp::complex StandardModel::gV_f(const Particle f) const
{
    return ( gA_f(f)
            *(1.0 - 4.0 * fabs(f.getCharge())*(kappaZ_f(f)) * sW2()));
}

gslpp::complex StandardModel::gA_f(const Particle f) const
{
    return ( sqrt(rhoZ_f(f)) * f.getIsospin());
}

gslpp::complex StandardModel::rhoZ_f(const Particle f) const
{
    if (f.getName().compare("TOP") == 0) return (gslpp::complex(0.0, 0.0, false));
    if (FlagRhoZ.compare("APPROXIMATEFORMULA") == 0)
        throw std::runtime_error("No approximate formula is available for rhoZ^f");
    else {

        if (FlagCacheInStandardModel)
            if (useRhoZ_f_cache[f.getIndex()])
                return rhoZ_f_cache[f.getIndex()];

        double myMw = Mw();

        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(myMw, DeltaRho);

        /* compute delta rho_rem^f */
        gslpp::complex deltaRho_remf[orders_EW_size];
        deltaRho_remf[EW1] = gslpp::complex(0.0, 0.0, false);
        deltaRho_remf[EW1QCD1] = gslpp::complex(0.0, 0.0, false);
        deltaRho_remf[EW1QCD2] = gslpp::complex(0.0, 0.0, false);
        deltaRho_remf[EW2] = gslpp::complex(0.0, 0.0, false);
        deltaRho_remf[EW2QCD1] = gslpp::complex(0.0, 0.0, false);
        deltaRho_remf[EW3] = gslpp::complex(0.0, 0.0, false);
        if (flag_order[EW1])
            deltaRho_remf[EW1] = myOneLoopEW->deltaRho_rem_f(f, myMw);
        if (flag_order[EW1QCD1])
#ifdef WITHIMTWOLOOPQCD
            deltaRho_remf[EW1QCD1] = gslpp::complex(myTwoLoopQCD->deltaRho_rem_f(f, myMw).real(),
                myTwoLoopQCD->deltaRho_rem_f(f, myMw).imag(), false);
#else
            deltaRho_remf[EW1QCD1] = gslpp::complex(myTwoLoopQCD->deltaRho_rem_f(f, myMw).real(), 0.0, false);
#endif
        if (flag_order[EW1QCD2])
            deltaRho_remf[EW1QCD2] = gslpp::complex(myThreeLoopQCD->deltaRho_rem_f(f, myMw).real(), 0.0, false);
        if (flag_order[EW2])
            deltaRho_remf[EW2] = gslpp::complex(myTwoLoopEW->deltaRho_rem_f(f, myMw).real(), 0.0, false);
        if (flag_order[EW2QCD1])
            deltaRho_remf[EW2QCD1] = gslpp::complex(myThreeLoopEW2QCD->deltaRho_rem_f(f, myMw).real(), 0.0, false);
        if (flag_order[EW3])
            deltaRho_remf[EW3] = gslpp::complex(myThreeLoopEW->deltaRho_rem_f(f, myMw).real(), 0.0, false);

        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(myMw);

        /* Re[rho_Z^f] with or without resummation */
        double deltaRho_rem_f_real[orders_EW_size];
        for (int j = 0; j < orders_EW_size; ++j)
            deltaRho_rem_f_real[j] = deltaRho_remf[j].real();
        double ReRhoZf = resumRhoZ(DeltaRho, deltaRho_rem_f_real, DeltaRbar_rem, f.is("BOTTOM"));

        /* Im[rho_Z^f] without resummation */
        double ImRhoZf = 0.0;
        for (int j = 0; j < orders_EW_size; ++j)
            ImRhoZf += deltaRho_remf[j].imag();

        rhoZ_f_cache[f.getIndex()] = gslpp::complex(ReRhoZf, ImRhoZf, false);
        useRhoZ_f_cache[f.getIndex()] = true;
        return (gslpp::complex(ReRhoZf, ImRhoZf, false));
    }
}

gslpp::complex StandardModel::kappaZ_f(const Particle f) const
{
    if (f.is("TOP")) return (gslpp::complex(0.0, 0.0, false));

    if (FlagCacheInStandardModel)
        if (useKappaZ_f_cache[f.getIndex()])
            return kappaZ_f_cache[f.getIndex()];

    double myMw = Mw();

    double ReKappaZf = 0.0, ImKappaZf = 0.0;
    if (FlagKappaZ.compare("APPROXIMATEFORMULA") == 0) {

//  Choose the correct formulae for the effective angle        
        if ( f.is("BOTTOM") ){
            ReKappaZf = myApproximateFormulae->sin2thetaEff_b_full() / sW2();            
        } else if ( f.is("ELECTRON") || f.is("MU") || f.is("TAU") ) {
            ReKappaZf = myApproximateFormulae->sin2thetaEff_l_full() / sW2();             
        } else {
            ReKappaZf = myApproximateFormulae->sin2thetaEff(f) / sW2();
        }
        
        ImKappaZf = myOneLoopEW->deltaKappa_rem_f(f, myMw).imag();
#ifdef WITHIMTWOLOOPQCD
        ImKappaZf += myTwoLoopQCD->deltaKappa_rem_f(f, myMw).imag();

        /* TEST */
        //ImKappaZf -= myCache->ale()*myCache->alsMz()/24.0/M_PI*(cW2() - sW2())/sW2()/sW2();
#endif
    } else {
        /* compute Delta rho */
        double DeltaRho[orders_EW_size];
        ComputeDeltaRho(myMw, DeltaRho);

        /* compute delta kappa_rem^f */
        gslpp::complex deltaKappa_remf[orders_EW_size];
        deltaKappa_remf[EW1] = gslpp::complex(0.0, 0.0, false);
        deltaKappa_remf[EW1QCD1] = gslpp::complex(0.0, 0.0, false);
        deltaKappa_remf[EW1QCD2] = gslpp::complex(0.0, 0.0, false);
        deltaKappa_remf[EW2] = gslpp::complex(0.0, 0.0, false);
        deltaKappa_remf[EW2QCD1] = gslpp::complex(0.0, 0.0, false);
        deltaKappa_remf[EW3] = gslpp::complex(0.0, 0.0, false);
        if (flag_order[EW1])
            deltaKappa_remf[EW1] = myOneLoopEW->deltaKappa_rem_f(f, myMw);
        if (flag_order[EW1QCD1])
#ifdef WITHIMTWOLOOPQCD
            deltaKappa_remf[EW1QCD1] = gslpp::complex(myTwoLoopQCD->deltaKappa_rem_f(f, myMw).real(),
                myTwoLoopQCD->deltaKappa_rem_f(f, myMw).imag(), false);
#else
            deltaKappa_remf[EW1QCD1] = gslpp::complex(myTwoLoopQCD->deltaKappa_rem_f(f, myMw).real(), 0.0, false);
#endif
        if (flag_order[EW1QCD2])
            deltaKappa_remf[EW1QCD2] = gslpp::complex(myThreeLoopQCD->deltaKappa_rem_f(f, myMw).real(), 0.0, false);
        if (flag_order[EW2])
            deltaKappa_remf[EW2] = gslpp::complex(myTwoLoopEW->deltaKappa_rem_f(f, myMw).real(), 0.0, false);
        if (flag_order[EW2QCD1])
            deltaKappa_remf[EW2QCD1] = gslpp::complex(myThreeLoopEW2QCD->deltaKappa_rem_f(f, myMw).real(), 0.0, false);
        if (flag_order[EW3])
            deltaKappa_remf[EW3] = gslpp::complex(myThreeLoopEW->deltaKappa_rem_f(f, myMw).real(), 0.0, false);

        /* compute Delta rbar_rem */
        double DeltaRbar_rem = 0.0;
        if (flag_order[EW1])
            DeltaRbar_rem = myOneLoopEW->DeltaRbar_rem(myMw);

        /* Re[kappa_Z^f] with or without resummation */
        double deltaKappa_rem_f_real[orders_EW_size];
        for (int j = 0; j < orders_EW_size; ++j)
            deltaKappa_rem_f_real[j] = deltaKappa_remf[j].real();

        ReKappaZf = resumKappaZ(DeltaRho, deltaKappa_rem_f_real, DeltaRbar_rem, f.is("BOTTOM"));

        /* O(alpha^2) correction to Re[kappa_Z^f] from the Z-gamma mixing */
        ReKappaZf += 35.0 * alphaMz() * alphaMz() / 18.0 / sW2()
                *(1.0 - 8.0 / 3.0 * ReKappaZf * sW2());

        /* Im[kappa_Z^f] without resummation */
        for (int j = 0; j < orders_EW_size; ++j)
            ImKappaZf += deltaKappa_remf[j].imag();
    }

    kappaZ_f_cache[f.getIndex()] = gslpp::complex(ReKappaZf, ImKappaZf, false);
    useKappaZ_f_cache[f.getIndex()] = true;
    return (gslpp::complex(ReKappaZf, ImKappaZf, false));
}

gslpp::complex StandardModel::deltaRhoZ_f(const Particle f) const
{
    Particle p1 = f, pe = leptons[ELECTRON];

    if (f.is("TOP") || f.is("ELECTRON")) return (gslpp::complex(0.0, 0.0, false));

    /* In the case of BOTTOM, the top contribution has to be subtracted.
     * The remaining contribution is the same as that for DOWN and STRANGE. */
    if (f.is("BOTTOM")) p1 = quarks[DOWN];

    double myMw = Mw();
    double cW2 = myMw * myMw / Mz / Mz, sW2 = 1.0 - cW2;

    gslpp::complex ul = (3.0 * myEWSMcache->v_f(pe, myMw) * myEWSMcache->v_f(pe, myMw)
            + myEWSMcache->a_f(pe) * myEWSMcache->a_f(pe)) / 4.0 / cW2 * myOneLoopEW->FZ(Mz*Mz, myMw)
            + myOneLoopEW->FW(Mz*Mz, pe, myMw);
    gslpp::complex uf = (3.0 * myEWSMcache->v_f(p1, myMw) * myEWSMcache->v_f(p1, myMw)
            + myEWSMcache->a_f(p1) * myEWSMcache->a_f(p1)) / 4.0 / cW2 * myOneLoopEW->FZ(Mz*Mz, myMw)
            + myOneLoopEW->FW(Mz*Mz, p1, myMw);

    gslpp::complex dRho = 2.0 * (uf - ul);
    dRho *= ale / 4.0 / M_PI / sW2;
    return dRho;
}

gslpp::complex StandardModel::deltaKappaZ_f(const Particle f) const
{
    Particle p1 = f, pe = leptons[ELECTRON];

    if (f.is("TOP") || f.is("ELECTRON")) return (gslpp::complex(0.0, 0.0, false));

    /* In the case of BOTTOM, the top contribution has to be subtracted.
     * The remaining contribution is the same as that for DOWN and STRANGE. */
    if (f.is("BOTTOM")) p1 = quarks[DOWN];

    double myMw = Mw();
    double cW2 = myMw * myMw / Mz / Mz, sW2 = 1.0 - cW2;
    gslpp::complex ul = (3.0 * myEWSMcache->v_f(pe, myMw) * myEWSMcache->v_f(pe, myMw)
            + myEWSMcache->a_f(pe) * myEWSMcache->a_f(pe)) / 4.0 / cW2 * myOneLoopEW->FZ(Mz*Mz, myMw)
            + myOneLoopEW->FW(Mz*Mz, pe, myMw);
    double deltal = myEWSMcache->delta_f(pe, myMw);
    gslpp::complex uf = (3.0 * myEWSMcache->v_f(p1, myMw) * myEWSMcache->v_f(p1, myMw)
            + myEWSMcache->a_f(p1) * myEWSMcache->a_f(p1)) / 4.0 / cW2 * myOneLoopEW->FZ(Mz*Mz, myMw)
            + myOneLoopEW->FW(Mz*Mz, p1, myMw);
    double deltaf = myEWSMcache->delta_f(p1, myMw);

    gslpp::complex dKappa = (deltaf * deltaf - deltal * deltal) / 4.0 / cW2 * myOneLoopEW->FZ(Mz*Mz, myMw)
            - uf + ul;
    dKappa *= ale / 4.0 / M_PI / sW2;
    return dKappa;
}


////////////////////////////////////////////////////////////////////////

double StandardModel::epsilon1() const
{
    double rhoZe = rhoZ_f(leptons[ELECTRON]).real();
    double DeltaRhoPrime = 2.0 * (sqrt(rhoZe) - 1.0);

    return DeltaRhoPrime;
}

double StandardModel::epsilon2() const
{
    double rhoZe = rhoZ_f(leptons[ELECTRON]).real();
    double sin2thetaEff = kappaZ_f(leptons[ELECTRON]).real() * sW2();
    double DeltaRhoPrime = 2.0 * (sqrt(rhoZe) - 1.0);
    double DeltaKappaPrime = sin2thetaEff / s02() - 1.0;
    double DeltaRW = 1.0 - M_PI * alphaMz() / (sqrt(2.0) * GF * Mz * Mz * sW2() * cW2());

    return ( c02() * DeltaRhoPrime + s02() * DeltaRW / (c02() - s02())
            - 2.0 * s02() * DeltaKappaPrime);
}

double StandardModel::epsilon3() const
{
    double rhoZe = rhoZ_f(leptons[ELECTRON]).real();
    double sin2thetaEff = kappaZ_f(leptons[ELECTRON]).real() * sW2();
    double DeltaRhoPrime = 2.0 * (sqrt(rhoZe) - 1.0);
    double DeltaKappaPrime = sin2thetaEff / s02() - 1.0;

    return ( c02() * DeltaRhoPrime + (c02() - s02()) * DeltaKappaPrime);
}

double StandardModel::epsilonb() const
{
    /* epsilon_b from g_A^b
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    //double rhoZe = rhoZ_l_SM(StandardModel::ELECTRON).real();
    //double rhoZb = rhoZ_q_SM(QCD::BOTTOM).real();
    //double DeltaRhoPrime = 2.0*( sqrt(rhoZe) - 1.0 );
    //double eps1 = DeltaRhoPrime;
    //return ( - 1.0 + sqrt(rhoZb)/(1.0 + eps1/2.0) );

    /* epsilon_b from Re(g_V^b/g_A^b), i.e. Re(kappaZ_b)
     * see Eq.(13) of IJMP A7, 1031 (1998) by Altarelli et al. */
    gslpp::complex kappaZe = kappaZ_f(leptons[ELECTRON]);
    gslpp::complex kappaZb = kappaZ_f(quarks[BOTTOM]);
    if (IsFlagWithoutNonUniversalVC())
        return ( kappaZe.real() / kappaZb.real() - 1.0);
    else
        return ( (kappaZe.real() + deltaKappaZ_f(quarks[BOTTOM]).real())
            / kappaZb.real() - 1.0);

    /* epsilon_b from Gamma_b via Eqs.(11), (12) and (16) of IJMP A7,
     * 1031 (1998) by Altarelli et al.
     * Note: mb has to be mb=4.7, since Eq.(16) were derived with this value.
     */
    //double als_Mz = Als(myCache->Mz(), FULLNNLO);
    //double delta_als = (als_Mz - 0.119)/M_PI;
    //double delta_alpha = (alphaMz() - 1.0/128.90)/myCache->ale();
    //double Gamma_b_Born = 0.3798*( 1.0 + delta_als - 0.42*delta_alpha);
    //double a = als_Mz/M_PI;
    //double RQCD = 1.0 + 1.2*a - 1.1*a*a - 13.0*a*a*a;
    //double mb = Mrun(myCache->Mz(), quarks[QCD::BOTTOM].getMass(), FULLNNLO);// This is wrong!
    //double mb = 4.7;
    //std::cout << "mb = " << mb << std::endl;
    //double beta = sqrt(1.0 - 4.0*mb*mb/myCache->Mz()/myCache->Mz());
    //double Nc = 3.0;
    //double factor = myCache->GF()*myCache->Mz()*myCache->Mz()*myCache->Mz()/6.0/M_PI/sqrt(2.0);
    //double Gamma_b = factor*beta*((3.0 - beta*beta)/2.0*gVq_SM(QCD::BOTTOM).abs2()
    //                              + beta*beta*gAq_SM(QCD::BOTTOM).abs2())
    //                 *Nc*RQCD*(1.0 + alphaMz()/12.0/M_PI);
    //return ( (Gamma_b/Gamma_b_Born - 1.0 - 1.42*epsilon1_SM()
    //          + 0.54*epsilon3_SM() )/2.29 );
}


////////////////////////////////////////////////////////////////////////

double StandardModel::resumMw(const double Mw_i, const double DeltaRho[orders_EW_size],
        const double DeltaR_rem[orders_EW_size]) const
{
    if ((FlagMw.compare("APPROXIMATEFORMULA") == 0)
            || (DeltaR_rem[EW2QCD1] != 0.0)
            || (DeltaR_rem[EW3] != 0.0))
        throw std::runtime_error("Error in StandardModel::resumMw()");

    if (!flag_order[EW2] && FlagMw.compare("NORESUM") != 0)
        throw std::runtime_error("Error in StandardModel::resumMw()");

    double cW2_TMP = Mw_i * Mw_i / Mz / Mz;
    double sW2_TMP = 1.0 - cW2_TMP;

    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G = 0.0;
    if (FlagMw.compare("NORESUM") == 0) {
        for (int j = 0; j < orders_EW_size; ++j) {
            DeltaRho_sum += DeltaRho[(orders_EW) j];
        }
    } else {
        // conversion: alpha(0) --> G_F
        f_AlphaToGF = sqrt(2.0) * GF * pow(Mz, 2.0) * sW2_TMP * cW2_TMP / M_PI / ale;
        DeltaRho_sum = f_AlphaToGF * DeltaRho[EW1]
                + f_AlphaToGF * DeltaRho[EW1QCD1]
                + f_AlphaToGF * DeltaRho[EW1QCD2]
                + pow(f_AlphaToGF, 2.0) * DeltaRho[EW2]
                + pow(f_AlphaToGF, 2.0) * DeltaRho[EW2QCD1]
                + pow(f_AlphaToGF, 3.0) * DeltaRho[EW3];
        DeltaRho_G = f_AlphaToGF * DeltaRho[EW1];
    }

    double R;
    double DeltaR_rem_sum = 0.0;
    double DeltaR_EW1 = 0.0, DeltaR_EW2_rem = 0.0;
    if (FlagMw.compare("NORESUM") == 0) {
        for (int j = 0; j < orders_EW_size; ++j)
            DeltaR_rem_sum += DeltaR_rem[(orders_EW) j];

        // Full EW one-loop contribution (without the full DeltaAlphaL5q)
        DeltaR_EW1 = -cW2_TMP / sW2_TMP * DeltaRho[EW1] + DeltaR_rem[EW1];

        // Full EW two-loop contribution without reducible corrections
        DeltaR_EW2_rem = myApproximateFormulae->DeltaR_TwoLoopEW_rem(Mw_i);

        // subtract the EW two-loop contributions from DeltaRho_sum and DeltaR_rem_sum
        DeltaRho_sum -= DeltaRho[EW2];
        DeltaR_rem_sum -= DeltaR_rem[EW2];

        // R = 1 + Delta r, including the full EW two-loop contribution
        R = 1.0 + DeltaAlphaL5q() - cW2_TMP / sW2_TMP * DeltaRho_sum
                + DeltaR_rem_sum;
        R += DeltaAlphaL5q() * DeltaAlphaL5q() + 2.0 * DeltaAlphaL5q() * DeltaR_EW1
                + DeltaR_EW2_rem;
    } else if (FlagMw.compare("OMSI") == 0) {
        // R = 1/(1 - Delta r)
        R = 1.0 / (1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum)
                / (1.0 - DeltaAlphaL5q()
                - DeltaR_rem[EW1] - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2]);
    } else if (FlagMw.compare("INTERMEDIATE") == 0) {
        // R = 1/(1 - Delta r)
        R = 1.0 / ((1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum)
                *(1.0 - DeltaAlphaL5q() - DeltaR_rem[EW1])
                - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2]);
    } else if (FlagMw.compare("OMSII") == 0) {
        // R = 1/(1 - Delta r)
        R = 1.0 / ((1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum)*(1.0 - DeltaAlphaL5q())
                - (1.0 + cW2_TMP / sW2_TMP * DeltaRho_G) * DeltaR_rem[EW1]
                - DeltaR_rem[EW1QCD1] - DeltaR_rem[EW2]);
    } else
        throw std::runtime_error("Error in StandardModel::resumMw()");

    if (FlagMw.compare("NORESUM") == 0) {
        /* Mzbar and Mwbar are defined in the complex-pole scheme. */

        double tmp = 4.0 * M_PI * ale / sqrt(2.0) / GF / Mzbar() / Mzbar();
        if (tmp * R > 1.0) throw std::runtime_error("StandardModel::resumMw(): Negative (1-tmp*R)");
        double Mwbar = Mzbar() / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp * R));

        return MwFromMwbar(Mwbar);
    } else {
        double tmp = 4.0 * M_PI * ale / sqrt(2.0) / GF / Mz / Mz;
        if (tmp * R > 1.0) throw std::runtime_error("StandardModel::resumMw(): Negative (1-tmp*R)");

        return (Mz / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - tmp * R)));
    }
}

double StandardModel::resumRhoZ(const double DeltaRho[orders_EW_size],
        const double deltaRho_rem[orders_EW_size],
        const double DeltaRbar_rem, const bool bool_Zbb) const
{
    if ((FlagRhoZ.compare("APPROXIMATEFORMULA") == 0)
            || (deltaRho_rem[EW1QCD2] != 0.0)
            || (deltaRho_rem[EW2QCD1] != 0.0)
            || (deltaRho_rem[EW3] != 0.0))
        throw std::runtime_error("Error in StandardModel::resumRhoZ()");

    if (!flag_order[EW2] && FlagRhoZ.compare("NORESUM") != 0)
        throw std::runtime_error("Error in StandardModel::resumRhoZ()");

    double Mw_TMP = Mw();
    double cW2_TMP = cW2();
    double sW2_TMP = sW2();

    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G;
    double DeltaRbar_rem_G, deltaRho_rem_G, deltaRho_rem_G2;
    // conversion: alpha(0) --> G_F
    f_AlphaToGF = sqrt(2.0) * GF * pow(Mz, 2.0)
            * sW2_TMP * cW2_TMP / M_PI / ale;
    DeltaRho_sum = f_AlphaToGF * DeltaRho[EW1]
            + f_AlphaToGF * DeltaRho[EW1QCD1]
            + f_AlphaToGF * DeltaRho[EW1QCD2]
            + pow(f_AlphaToGF, 2.0) * DeltaRho[EW2]
            + pow(f_AlphaToGF, 2.0) * DeltaRho[EW2QCD1]
            + pow(f_AlphaToGF, 3.0) * DeltaRho[EW3];
    DeltaRho_G = f_AlphaToGF * DeltaRho[EW1];
    DeltaRbar_rem_G = f_AlphaToGF*DeltaRbar_rem;
    deltaRho_rem_G = f_AlphaToGF * (deltaRho_rem[EW1]
            + deltaRho_rem[EW1QCD1]);
    deltaRho_rem_G2 = pow(f_AlphaToGF, 2.0) * deltaRho_rem[EW2];

    /* Real parts */
    double rhoZ;
    if (!bool_Zbb) {
        if (FlagRhoZ.compare("OMSI") == 0) {
            rhoZ = (1.0 + deltaRho_rem_G + deltaRho_rem_G2)
                    / (1.0 - DeltaRho_sum * (1.0 - DeltaRbar_rem_G));
        } else if (FlagRhoZ.compare("INTERMEDIATE") == 0) {
            rhoZ = (1.0 + deltaRho_rem_G)
                    / (1.0 - DeltaRho_sum * (1.0 - DeltaRbar_rem_G))
                    + deltaRho_rem_G2;
        } else if (FlagRhoZ.compare("NORESUM") == 0
                || FlagRhoZ.compare("OMSII") == 0) {
            rhoZ = 1.0 + DeltaRho_sum - DeltaRho_G * DeltaRbar_rem_G
                    + DeltaRho_G * DeltaRho_G
                    + deltaRho_rem_G * (1.0 + DeltaRho_G) + deltaRho_rem_G2;
        } else
            throw std::runtime_error("Error in StandardModel::resumRhoZ()");
    } else {
        /* Z to bb */
        double OnePlusTaub = 1.0 + taub();
        double OnePlusTaub2 = OnePlusTaub*OnePlusTaub;
        double rhoZbL;
        deltaRho_rem_G += f_AlphaToGF * ale / 4.0 / M_PI / sW2_TMP
                * pow(mtpole / Mw_TMP, 2.0);
        if (FlagRhoZ.compare("NORESUM") == 0) {
            rhoZ = (1.0 + DeltaRho_sum - DeltaRho_G * DeltaRbar_rem_G
                    + DeltaRho_G * DeltaRho_G
                    + deltaRho_rem_G * (1.0 + DeltaRho_G) + deltaRho_rem_G2)
                    * OnePlusTaub2;
        } else if (FlagRhoZ.compare("OMSI") == 0) {
            rhoZbL = OnePlusTaub2 / (1.0 - DeltaRho_sum);
            rhoZ = rhoZbL / (1.0 - rhoZbL * deltaRho_rem_G);
        } else if (FlagRhoZ.compare("INTERMEDIATE") == 0) {
            rhoZbL = OnePlusTaub2 / (1.0 - DeltaRho_sum);
            rhoZ = rhoZbL * (1.0 + rhoZbL * deltaRho_rem_G);
        } else if (FlagRhoZ.compare("OMSII") == 0) {
            rhoZbL = OnePlusTaub2 / (1.0 - DeltaRho_sum);
            rhoZ = rhoZbL * (1.0 + deltaRho_rem_G);
        } else
            throw std::runtime_error("Error in StandardModel::resumRhoZ()");
    }

    return rhoZ;
}

double StandardModel::resumKappaZ(const double DeltaRho[orders_EW_size],
        const double deltaKappa_rem[orders_EW_size],
        const double DeltaRbar_rem, const bool bool_Zbb) const
{
    if ((FlagKappaZ.compare("APPROXIMATEFORMULA") == 0)
            || (deltaKappa_rem[EW2QCD1] != 0.0)
            || (deltaKappa_rem[EW3] != 0.0))
        throw std::runtime_error("Error in StandardModel::resumKappaZ()");

    if (!flag_order[EW2] && FlagKappaZ.compare("NORESUM") != 0)
        throw std::runtime_error("Error in StandardModel::resumKappaZ()");

    double Mw_TMP = Mw();
    double cW2_TMP = cW2();
    double sW2_TMP = sW2();

    double f_AlphaToGF, DeltaRho_sum = 0.0, DeltaRho_G;
    double DeltaRbar_rem_G, deltaKappa_rem_G, deltaKappa_rem_G2;
    // conversion: alpha(0) --> G_F
    f_AlphaToGF = sqrt(2.0) * GF * pow(Mz, 2.0)
            * sW2_TMP * cW2_TMP / M_PI / ale;
    DeltaRho_sum = f_AlphaToGF * DeltaRho[EW1]
            + f_AlphaToGF * DeltaRho[EW1QCD1]
            + f_AlphaToGF * DeltaRho[EW1QCD2]
            + pow(f_AlphaToGF, 2.0) * DeltaRho[EW2]
            + pow(f_AlphaToGF, 2.0) * DeltaRho[EW2QCD1]
            + pow(f_AlphaToGF, 3.0) * DeltaRho[EW3];
    DeltaRho_G = f_AlphaToGF * DeltaRho[EW1];
    DeltaRbar_rem_G = f_AlphaToGF*DeltaRbar_rem;
    deltaKappa_rem_G = f_AlphaToGF * (deltaKappa_rem[EW1]
            + deltaKappa_rem[EW1QCD1]
            + deltaKappa_rem[EW1QCD2]);
    deltaKappa_rem_G2 = pow(f_AlphaToGF, 2.0) * deltaKappa_rem[EW2];

    /* Real parts */
    double kappaZ;
    if (!bool_Zbb) {
        if (FlagKappaZ.compare("OMSI") == 0) {
            kappaZ = (1.0 + deltaKappa_rem_G + deltaKappa_rem_G2)
                    *(1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum * (1.0 - DeltaRbar_rem_G));
        } else if (FlagKappaZ.compare("INTERMEDIATE") == 0) {
            kappaZ = (1.0 + deltaKappa_rem_G)
                    *(1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum * (1.0 - DeltaRbar_rem_G))
                    + deltaKappa_rem_G2;
        } else if (FlagKappaZ.compare("NORESUM") == 0
                || FlagKappaZ.compare("OMSII") == 0) {
            kappaZ = 1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum
                    - cW2_TMP / sW2_TMP * DeltaRho_G * DeltaRbar_rem_G
                    + deltaKappa_rem_G * (1.0 + cW2_TMP / sW2_TMP * DeltaRho_G)
                    + deltaKappa_rem_G2;
        } else
            throw std::runtime_error("Error in StandardModel::resumKappaZ()");
    } else {
        /* Z to bb */
        double OnePlusTaub = 1.0 + taub();
        double kappaZbL;
        deltaKappa_rem_G -= f_AlphaToGF * ale / 8.0 / M_PI / sW2_TMP
                * pow(mtpole / Mw_TMP, 2.0);
        if (FlagKappaZ.compare("NORESUM") == 0) {
            kappaZ = (1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum
                    - cW2_TMP / sW2_TMP * DeltaRho_G * DeltaRbar_rem_G
                    + deltaKappa_rem_G * (1.0 + cW2_TMP / sW2_TMP * DeltaRho_G)
                    + deltaKappa_rem_G2) / OnePlusTaub;
        } else if (FlagKappaZ.compare("OMSI") == 0) {
            kappaZbL = (1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum) / OnePlusTaub;
            kappaZ = kappaZbL * (1.0 + deltaKappa_rem_G);
        } else if (FlagKappaZ.compare("INTERMEDIATE") == 0
                || FlagKappaZ.compare("OMSII") == 0) {
            kappaZbL = (1.0 + cW2_TMP / sW2_TMP * DeltaRho_sum) / OnePlusTaub;
            kappaZ = kappaZbL + deltaKappa_rem_G;
        } else
            throw std::runtime_error("Error in StandardModel::resumKappaZ()");
    }

    return kappaZ;
}

double StandardModel::taub() const
{
    double taub_tmp = 0.0;
    double Xt = myEWSMcache->Xt_GF();
    if (flag_order[EW1])
        taub_tmp += -2.0 * Xt;
    if (flag_order[EW1QCD1])
        taub_tmp += 2.0 / 3.0 * M_PI * Xt * myEWSMcache->alsMt();
    if (flag_order[EW1QCD2])
        taub_tmp += 0.0;
    if (flag_order[EW2])
        taub_tmp += -2.0 * Xt * Xt * myTwoLoopEW->tau_2();
    if (flag_order[EW2QCD1])
        taub_tmp += 0.0;
    if (flag_order[EW3])
        taub_tmp += 0.0;

    return taub_tmp;
}

double StandardModel::Delta_EWQCD(const QCD::quark q) const
{
    switch (q) {
        case QCD::UP:
        case QCD::CHARM:
            return ( -0.000113);
        case QCD::TOP:
            return ( 0.0);
        case QCD::DOWN:
        case QCD::STRANGE:
            return ( -0.000160);
        case QCD::BOTTOM:
            return ( -0.000040);
        default:
            throw std::runtime_error("Error in StandardModel::Delta_EWQCD");
    }
}

double StandardModel::RVq(const QCD::quark q) const
{
    if (q == QCD::TOP) return 0.0;

    double mcMz, mbMz;
    mcMz = myEWSMcache->mf(getQuarks(CHARM), Mz, FULLNNLO);
    mbMz = myEWSMcache->mf(getQuarks(BOTTOM), Mz, FULLNNLO);
    //mcMz = 0.56381685; /* for debug */
    //mbMz = 2.8194352; /* for debug */

    double MtPole = mtpole;

    /* electric charge squared */
    double Qf2 = pow(quarks[q].getCharge(), 2.0);

    /* s = Mz^2 */
    double s = Mz * Mz;

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch (q) {
        case QCD::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case QCD::BOTTOM:
            mqMz2 = mbMz*mbMz;
            mqdash4 = mcMz2*mcMz2;
            break;
        default:
            mqMz2 = 0.0;
            mqdash4 = 0.0;
            break;
    }

    /* Logarithms */
    //double log_t = log(pow(quarks[TOP].getMass(),2.0)/s);
    double log_t = log(MtPole * MtPole / s); // the pole mass
    double log_c = log(mcMz2 / s);
    double log_b = log(mbMz2 / s);
    double log_q;
    switch (q) {
        case QCD::CHARM:
        case QCD::BOTTOM:
            log_q = log(mqMz2 / s);
            break;
        default:
            log_q = 0.0;
            break;
    }

    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = getMyEWSMcache()->getZeta2();
    double zeta3 = getMyEWSMcache()->getZeta3();
    //double zeta4 = getMyCache()->GetZeta4();
    double zeta5 = getMyEWSMcache()->getZeta5();

    /* massless non-singlet corrections */
    double C02 = 365.0 / 24.0 - 11.0 * zeta3 + (-11.0 / 12.0 + 2.0 / 3.0 * zeta3) * nf;
    double C03 = 87029.0 / 288.0 - 121.0 / 8.0 * zeta2 - 1103.0 / 4.0 * zeta3
            + 275.0 / 6.0 * zeta5
            + (-7847.0 / 216.0 + 11.0 / 6.0 * zeta2 + 262.0 / 9.0 * zeta3
            - 25.0 / 9.0 * zeta5) * nf
            + (151.0 / 162.0 - zeta2 / 18.0 - 19.0 / 27.0 * zeta3) * nf*nf;
    double C04 = -156.61 + 18.77 * nf - 0.7974 * nf * nf + 0.0215 * nf * nf*nf;
    //std::cout << "TEST: C02 = " << C02 << std::endl;// TEST (should be 1.40923)
    //std::cout << "TEST: C03 = " << C03 << std::endl;// TEST (should be -12.7671)
    //std::cout << "TEST: C04 = " << C04 << std::endl;// TEST (should be -80.0075)

    /* quadratic massive corrections */
    double C23 = -80.0 + 60.0 * zeta3 + (32.0 / 9.0 - 8.0 / 3.0 * zeta3) * nf;
    double C21V = 12.0;
    double C22V = 253.0 / 2.0 - 13.0 / 3.0 * nf;
    double C23V = 2522.0 - 855.0 / 2.0 * zeta2 + 310.0 / 3.0 * zeta3 - 5225.0 / 6.0 * zeta5
            + (-4942.0 / 27.0 + 34.0 * zeta2 - 394.0 / 27.0 * zeta3
            + 1045.0 / 27.0 * zeta5) * nf
            + (125.0 / 54.0 - 2.0 / 3.0 * zeta2) * nf*nf;

    /* quartic massive corrections */
    double C42 = 13.0 / 3.0 - 4.0 * zeta3;
    double C40V = -6.0;
    double C41V = -22.0;
    double C42V = -3029.0 / 12.0 + 162.0 * zeta2 + 112.0 * zeta3
            + (143.0 / 18.0 - 4.0 * zeta2 - 8.0 / 3.0 * zeta3) * nf;
    double C42VL = -11.0 / 2.0 + nf / 3.0;

    /* power suppressed top-mass correction */
    //double xt = s/pow(quarks[TOP].getMass(),2.0);
    double xt = s / MtPole / MtPole; // the pole mass
    double C2t = xt * (44.0 / 675.0 - 2.0 / 135.0 * (-log_t));

    /* rescaled strong coupling constant */
    double AlsMzPi = AlsMz / M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    /* electromagnetic coupling at Mz */
    double alpMz = alphaMz();

    /* radiator function to the vector current */
    double RVf;
    RVf = 1.0 + 3.0 / 4.0 * Qf2 * alpMz / M_PI + AlsMzPi - Qf2 / 4.0 * alpMz / M_PI * AlsMzPi
            + (C02 + C2t) * AlsMzPi2 + C03 * AlsMzPi3 + C04 * AlsMzPi4
            + (mcMz2 + mbMz2) / s * C23 * AlsMzPi3
            + mqMz2 / s * (C21V * AlsMzPi + C22V * AlsMzPi2 + C23V * AlsMzPi3)
            + mcMz2 * mcMz2 / s / s * (C42 - log_c) * AlsMzPi2
            + mbMz2 * mbMz2 / s / s * (C42 - log_b) * AlsMzPi2
            + mqMz2 * mqMz2 / s / s * (C40V + C41V * AlsMzPi + (C42V + C42VL * log_q) * AlsMzPi2)
            + 12.0 * mqdash4 / s / s * AlsMzPi2
            - mqMz2 * mqMz2 * mqMz2 / s / s / s
            * (8.0 + 16.0 / 27.0 * (155.0 + 6.0 * log_q) * AlsMzPi);
    return RVf;
}

double StandardModel::RAq(const QCD::quark q) const
{
    if (q == QCD::TOP) return 0.0;

    double mcMz, mbMz;
    mcMz = myEWSMcache->mf(getQuarks(CHARM), Mz, FULLNNLO);
    mbMz = myEWSMcache->mf(getQuarks(BOTTOM), Mz, FULLNNLO);
    //mcMz = 0.56381685; /* for debug */
    //mbMz = 2.8194352; /* for debug */

    double MtPole = mtpole;

    /* z-component of isospin */
    double I3q = quarks[q].getIsospin();
    /* electric charge squared */
    double Qf2 = pow(quarks[q].getCharge(), 2.0);

    /* s = Mz^2 */
    double s = Mz * Mz;

    /* products of the charm and bottom masses at Mz */
    double mcMz2 = mcMz*mcMz;
    double mbMz2 = mbMz*mbMz;
    double mqMz2, mqdash4;
    switch (q) {
        case QCD::CHARM:
            mqMz2 = mcMz*mcMz;
            mqdash4 = mbMz2*mbMz2;
            break;
        case QCD::BOTTOM:
            mqMz2 = mbMz*mbMz;
            mqdash4 = mcMz2*mcMz2;
            break;
        default:
            mqMz2 = 0.0;
            mqdash4 = 0.0;
            break;
    }

    /* Logarithms */
    //double log_t = log(pow(quarks[TOP].getMass(),2.0)/s);
    double log_t = log(MtPole * MtPole / s); // the pole mass
    double log_c = log(mcMz2 / s);
    double log_b = log(mbMz2 / s);
    double log_q;
    switch (q) {
        case QCD::CHARM:
        case QCD::BOTTOM:
            log_q = log(mqMz2 / s);
            break;
        default:
            log_q = 0.0;
            break;
    }

    /* the active number of flavour */
    double nf = 5.0;

    /* zeta functions */
    double zeta2 = getMyEWSMcache()->getZeta2();
    double zeta3 = getMyEWSMcache()->getZeta3();
    double zeta4 = getMyEWSMcache()->getZeta4();
    double zeta5 = getMyEWSMcache()->getZeta5();

    /* massless non-singlet corrections */
    double C02 = 365.0 / 24.0 - 11.0 * zeta3 + (-11.0 / 12.0 + 2.0 / 3.0 * zeta3) * nf;
    double C03 = 87029.0 / 288.0 - 121.0 / 8.0 * zeta2 - 1103.0 / 4.0 * zeta3
            + 275.0 / 6.0 * zeta5
            + (-7847.0 / 216.0 + 11.0 / 6.0 * zeta2 + 262.0 / 9.0 * zeta3
            - 25.0 / 9.0 * zeta5) * nf
            + (151.0 / 162.0 - zeta2 / 18.0 - 19.0 / 27.0 * zeta3) * nf*nf;
    double C04 = -156.61 + 18.77 * nf - 0.7974 * nf * nf + 0.0215 * nf * nf*nf;
    //std::cout << "TEST: C02 = " << C02 << std::endl;// TEST (should be 1.40923)
    //std::cout << "TEST: C03 = " << C03 << std::endl;// TEST (should be -12.7671)
    //std::cout << "TEST: C04 = " << C04 << std::endl;// TEST (should be -80.0075)

    /* quadratic massive corrections */
    double C23 = -80.0 + 60.0 * zeta3 + (32.0 / 9.0 - 8.0 / 3.0 * zeta3) * nf;
    double C20A = -6.0;
    double C21A = -22.0;
    double C22A = -8221.0 / 24.0 + 57.0 * zeta2 + 117.0 * zeta3
            + (151.0 / 12.0 - 2.0 * zeta2 - 4.0 * zeta3) * nf;
    double C23A = -4544045.0 / 864.0 + 1340.0 * zeta2 + 118915.0 / 36.0 * zeta3
            - 127.0 * zeta5
            + (71621.0 / 162.0 - 209.0 / 2.0 * zeta2 - 216.0 * zeta3
            + 5.0 * zeta4 + 55.0 * zeta5) * nf
            + (-13171.0 / 1944.0 + 16.0 / 9.0 * zeta2 + 26.0 / 9.0 * zeta3) * nf*nf;

    /* quartic massive corrections */
    double C42 = 13.0 / 3.0 - 4.0 * zeta3;
    double C40A = 6.0;
    double C41A = 10.0;
    double C42A = 3389.0 / 12.0 - 162.0 * zeta2 - 220.0 * zeta3
            + (-41.0 / 6.0 + 4.0 * zeta2 + 16.0 / 3.0 * zeta3) * nf;
    double C42AL = 77.0 / 2.0 - 7.0 / 3.0 * nf;

    /* power suppressed top-mass correction */
    //double xt = s/pow(quarks[TOP].getMass(),2.0);
    double xt = s / MtPole / MtPole; // the pole mass
    double C2t = xt * (44.0 / 675.0 - 2.0 / 135.0 * (-log_t));

    /* singlet axial-vector corrections */
    double I2 = -37.0 / 12.0 + (-log_t) + 7.0 / 81.0 * xt + 0.0132 * xt*xt;
    double I3 = -5075.0 / 216.0 + 23.0 / 6.0 * zeta2 + zeta3 + 67.0 / 18.0 * (-log_t)
            + 23.0 / 12.0 * log_t*log_t;
    double I4 = 49.0309 - 17.6637 * (-log_t) + 14.6597 * log_t * log_t
            + 3.6736 * (-log_t * log_t * log_t);

    /* rescaled strong coupling constant */
    double AlsMzPi = AlsMz / M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    /* electromagnetic coupling at Mz */
    double alpMz = alphaMz();

    /* radiator function to the axial-vector current */
    double RAf;
    RAf = 1.0 + 3.0 / 4.0 * Qf2 * alpMz / M_PI + AlsMzPi - Qf2 / 4.0 * alpMz / M_PI * AlsMzPi
            + (C02 + C2t - 2.0 * I3q * I2) * AlsMzPi2
            + (C03 - 2.0 * I3q * I3) * AlsMzPi3
            + (C04 - 2.0 * I3q * I4) * AlsMzPi4
            + (mcMz2 + mbMz2) / s * C23 * AlsMzPi3
            + mqMz2 / s * (C20A + C21A * AlsMzPi + C22A * AlsMzPi2
            + 6.0 * (3.0 + log_t) * AlsMzPi2 + C23A * AlsMzPi3)
            //- 10.0*mqMz2/pow(quarks[TOP].getMass(),2.0)
            - 10.0 * mqMz2 / MtPole / MtPole // the pole mass
            * (8.0 / 81.0 + log_t / 54.0) * AlsMzPi2
            + mcMz2 * mcMz2 / s / s * (C42 - log_c) * AlsMzPi2
            + mbMz2 * mbMz2 / s / s * (C42 - log_b) * AlsMzPi2
            + mqMz2 * mqMz2 / s / s * (C40A + C41A * AlsMzPi
            + (C42A + C42AL * log_q) * AlsMzPi2)
            - 12.0 * mqdash4 / s / s*AlsMzPi2;
    return RAf;
}

double StandardModel::RVh() const
{
    /* rescaled strong coupling constant */
    double AlsMzPi = AlsMz / M_PI;
    double AlsMzPi2 = AlsMzPi*AlsMzPi;
    double AlsMzPi3 = AlsMzPi2*AlsMzPi;
    double AlsMzPi4 = AlsMzPi3*AlsMzPi;

    gslpp::complex gV_sum(0.0, 0.0);
    gslpp::complex gV_q;
    for (int q = 0; q < 6; q++) {
        gV_q = gV_f(QCD::quarks[(QCD::quark)q]);
        if (q == (int) (QCD::TOP))
            gV_q = 0.0;
        gV_sum += gV_q;
    }

    // singlet vector corrections
    return ( gV_sum.abs2()*(-0.4132 * AlsMzPi3 - 4.9841 * AlsMzPi4));
}

/* BEGIN: REMOVE FROM THE PACKAGE */
////////////////////////////////////////////////////////////////////////////////////
//LEP2 Observables


double StandardModel::LEP2sigmaMu(const double s) const
{
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-8;
    double abserr = 1.e-20;
    
    if(s == 130.*130.){
    
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        } else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 136.*136.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        } else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 161.*161.){
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 172.*172.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 183.*183.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 189.*189.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 192.*192.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 196.*196.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 200.*200.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 202.*202.) {
         if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 205.*205.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 207.*207.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(MU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_mu207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2AFBmu!");
    }
        
   
    double sigma_mu = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    
    return sigma_mu;    

}


double StandardModel::LEP2sigmaTau(const double s) const
{
    
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;
    
    if(s == 130.*130.){
    
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        } else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 136.*136.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        } else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 161.*161.){
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 172.*172.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 183.*183.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 189.*189.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 192.*192.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 196.*196.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 200.*200.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 202.*202.) {
         if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 205.*205.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 207.*207.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_l(QCD::lepton(TAU), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_tau207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2sigmaTau!");
    }

    
    double sigma_tau = SMresult_cache;
    
    gsl_set_error_handler(old_handler);

    
    return sigma_tau;
}


double StandardModel::LEP2sigmaCharm(const double s) const
{
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-8;
    double abserr = 1.e-20;
    
    if(s == 133.*133.){
    
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        } else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 167.*167.){
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 183.*183.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 189.*189.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 192.*192.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 196.*196.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 200.*200.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 202.*202.) {
         if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 205.*205.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 207.*207.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(CHARM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2sigmaCharm!");
    }
        
   
    double sigma_mu = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    
    return sigma_mu;    

}


double StandardModel::LEP2sigmaBottom(const double s) const
{
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-8;
    double abserr = 1.e-20;
    
    if(s == 133.*133.){
    
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 167.*167.){
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 183.*183.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 189.*189.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 192.*192.) { 
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 196.*196.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 200.*200.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 202.*202.) {
         if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 205.*205.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 207.*207.) {
        if (!flagLEP2[ISR])
            SMresult_cache = sigma_NoISR_q(QCD::quark(BOTTOM), s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            SMresult_cache += sigma_box;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2sigmaBottom!");
    }
        
   
    double sigma_mu = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    
    return sigma_mu;    

}


double StandardModel::LEP2sigmaHadron(const double s) const
{
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-8;
    double abserr = 1.e-20;
    
    if(s == 130.*130.){
    
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
      
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 133.*133.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    }  else if (s == 136.*136.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 161.*161.){
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, 1.e-12, 1.e-6, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 200, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    }  else if (s == 167.*167.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, 1.e-15, 1.e-9, 200, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, 1.e-15, 1.e-9, 200, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, 1.e-15, 1.e-9, 200, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, 1.e-15, 1.e-9, 200, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 172.*172.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;      
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 183.*183.) { 
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 189.*189.) { 
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 192.*192.) { 
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm192, &(*this), _1));

            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange192, &(*this), _1));

            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 196.*196.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange196, &(*this), _1));

            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 200.*200.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 202.*202.) {
         if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom202, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 205.*205.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom205, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else if (s == 207.*207.) {
        if (!flagLEP2[ISR]){
            SMresult_cache = sigma_NoISR_q(QCD::quark(UP), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(DOWN), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(CHARM), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(STRANGE), s);
            SMresult_cache += sigma_NoISR_q(QCD::quark(BOTTOM), s);
        } else {
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_up207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_down207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_charm207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_strange207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_sigmaWithISR_bottom207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            SMresult_cache += average;
        }
        
        if (flagLEP2[WeakBox]) {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_up207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box = average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_down207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_strange207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom207, &(*this), _1));
 


            if (gsl_integration_qags(&f_GSL, -1., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            sigma_box += average;
            SMresult_cache += sigma_box;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2sigmaHadron!");
    }
        
   
    double sigma_had = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    
    return sigma_had;    

}


double StandardModel::LEP2AFBbottom(const double s) const
{

    bSigmaForAFB = true;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;
    
    if(s == 133.*133.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 167.*167.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 183.*183.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 189.*189.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 192.*192.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 196.*196.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 200.*200.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 202.*202.) {
         double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 205.*205.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 207.*207.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(BOTTOM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_bottom207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaBottom(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_bottom207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaBottom(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2AFBbottom!");
    }
        
   
    double AFBbottom = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    bSigmaForAFB = false;
    return AFBbottom;    

}


double StandardModel::LEP2AFBcharm(const double s) const
{

    bSigmaForAFB = true;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;
    
    if(s == 133.*133.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm133, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 167.*167.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm167, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 183.*183.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 189.*189.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 192.*192.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 196.*196.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 200.*200.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 202.*202.) {
         double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 205.*205.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 207.*207.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_q(QCD::quark(CHARM),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_charm205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaCharm(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_charm207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaCharm(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2AFBcharm!");
    }
        
   
    double AFBcharm = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    bSigmaForAFB = false;
    return AFBcharm;    

}


double StandardModel::LEP2AFBmu(const double s) const
{

    bSigmaForAFB = true;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;
    
    if(s == 130.*130.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 136.*136.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 161.*161.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 172.*172.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 183.*183.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 189.*189.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 192.*192.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 196.*196.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 200.*200.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 202.*202.) {
         double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 205.*205.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 207.*207.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(MU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_mu207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaMu(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_mu207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaMu(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::AFBmu!");
    }
        
   
    double AFBmu = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    bSigmaForAFB = false;
    return AFBmu;    

}


double StandardModel::LEP2AFBtau(const double s) const
{

    bSigmaForAFB = true;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;
    
    if(s == 130.*130.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau130, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 136.*136.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau136, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 161.*161.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau161, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 172.*172.){
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau172, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 183.*183.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau183, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 189.*189.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau189, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 192.*192.) { 
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau192, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 196.*196.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau196, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 200.*200.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau200, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 202.*202.) {
         double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau202, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 205.*205.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau205, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else if (s == 207.*207.) {
        double AFB_noBox, sigma = 0.0;
        if (!flagLEP2[ISR])
            AFB_noBox = AFB_NoISR_l(QCD::lepton(TAU),s);
        else {
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_AFBnumeratorWithISR_tau207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1.-0.85*0.85, abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double numerator = average; // interval
          
             
            sigma = LEP2sigmaTau(s);
            
            AFB_noBox = numerator/sigma;
        }
        SMresult_cache = AFB_noBox;
        
        if (flagLEP2[WeakBox]) {
            // numerator
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, 0., 1., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_F = average; // interval
            f_GSL = convertToGslFunction(bind(&StandardModel::getIntegrand_dsigmaBox_tau207, &(*this), _1));
            if (gsl_integration_qags(&f_GSL, -1., 0., abserr, relerr, 100, w_GSL1, &average, &error) != 0){ 
                SMresult_cache = std::numeric_limits<double>::quiet_NaN();
            }
            double sigma_box_B = average; // interval
            
            // denominator
            if (!flagLEP2[ISR]) {
                
                sigma = LEP2sigmaTau(s);
            }
            
            SMresult_cache += (sigma_box_F - sigma_box_B)/sigma;
        }
    } else {
        throw std::runtime_error("ERROR: wrong LEP2 energy in StandardModel::LEP2AFBtau!");
    }
        
   
    double AFBtau = SMresult_cache;
    
    gsl_set_error_handler(old_handler);
    bSigmaForAFB = false;
    return AFBtau;    

}


double StandardModel::LEP2Rbottom(const double s) const
{

    double sigma_b = LEP2sigmaBottom(s);
    double sigma_had = LEP2sigmaHadron(s);
    SMresult_cache =  sigma_b / sigma_had;
    double R_bottom = SMresult_cache;
    
    return R_bottom;    
}


double StandardModel::LEP2Rcharm(const double s) const
{

    double sigma_c = LEP2sigmaCharm(s);
    double sigma_had = LEP2sigmaHadron(s);
    SMresult_cache =  sigma_c / sigma_had;
    double R_charm = SMresult_cache;
    
    return R_charm;    
}


double StandardModel::sigma_NoISR_l(const QCD::lepton l_flavor, const double s) const
{
    double ml = getLeptons(l_flavor).getMass();
    double l_charge = getLeptons(l_flavor).getCharge();
    double sigma = myTwoFermionsLEP2->sigma_l(l_flavor, ml, s, Mw(), Gamma_Z(), flagLEP2[Weak]);

    if (!bSigmaForAFB && flagLEP2[QEDFSR])
        sigma *= myTwoFermionsLEP2->QED_FSR_forSigma(s, l_charge);

    return sigma;
}   
    
double StandardModel::sigma_NoISR_q(const QCD::quark q_flavor, const double s) const
{
    double mq = m_q(q_flavor, sqrt(s));
    double q_charge = getQuarks(q_flavor).getCharge();
    double sigma = myTwoFermionsLEP2->sigma_q(q_flavor, mq, s, Mw(), Gamma_Z(), flagLEP2[Weak]);
    
    if (!bSigmaForAFB && flagLEP2[QEDFSR])
        sigma *= myTwoFermionsLEP2->QED_FSR_forSigma(s, q_charge);
        
    if (!bSigmaForAFB && flagLEP2[QCDFSR])
        sigma *= myTwoFermionsLEP2->QCD_FSR_forSigma(s);

    return sigma;
}      
    
double StandardModel::AFB_NoISR_l(const QCD::lepton l_flavor, const double s) const
{
    double ml = getLeptons(l_flavor).getMass();
    double AFB = myTwoFermionsLEP2->AFB_l(l_flavor, ml, s, Mw(), Gamma_Z(), flagLEP2[Weak]);

    return AFB;
}
    
double StandardModel::AFB_NoISR_q(const QCD::quark q_flavor, const  double s) const
{
    double mq = m_q(q_flavor, sqrt(s));
    double AFB = myTwoFermionsLEP2->AFB_q(q_flavor, mq, s, Mw(), Gamma_Z(), flagLEP2[Weak]);
        
    if (flagLEP2[QCDFSR])
        AFB *= myTwoFermionsLEP2->QCD_FSR_forAFB(q_flavor, mq, s);
        
    return AFB;
}
    
double StandardModel::Integrand_sigmaWithISR_l(double x, const QCD::lepton l_flavor, const  double s) const
{
    double sprime = (1.0 - x)*s;
    double ml = getLeptons(l_flavor).getMass();
    double l_charge = getLeptons(l_flavor).getCharge();
    double sigma = myTwoFermionsLEP2->sigma_l(l_flavor, ml, sprime, Mw(), Gamma_Z(), 
                                             flagLEP2[Weak]);
    double H = myTwoFermionsLEP2->H_ISR(x, s);
    
    if (!bSigmaForAFB && flagLEP2[QEDFSR])
        sigma *= myTwoFermionsLEP2->QED_FSR_forSigma(sprime, l_charge);

    return ( H*sigma );
}   
    
double StandardModel::getIntegrand_sigmaWithISR_mu130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_mu207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}
    

double StandardModel::getIntegrand_sigmaWithISR_tau130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_tau136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_sigmaWithISR_tau207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}
    
double StandardModel::Integrand_sigmaWithISR_q(double x, const QCD::quark q_flavor, const  double s) const
{
    double sprime = (1.0 - x)*s;
    double mq = m_q(q_flavor, sqrt(s));
    double q_charge = getQuarks(q_flavor).getCharge();
    double sigma = myTwoFermionsLEP2->sigma_q(q_flavor, mq, sprime, Mw(), Gamma_Z(), 
                                             flagLEP2[Weak]);
    double H = myTwoFermionsLEP2->H_ISR(x, s);
    
    if (!bSigmaForAFB && flagLEP2[QEDFSR])
        sigma *= myTwoFermionsLEP2->QED_FSR_forSigma(sprime, q_charge);
        
    if (!bSigmaForAFB && flagLEP2[QCDFSR])
        sigma *= myTwoFermionsLEP2->QCD_FSR_forSigma(sprime);

    return ( H*sigma );
}    
    



//up


double StandardModel::getIntegrand_sigmaWithISR_up130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_up133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_sigmaWithISR_up207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}


//down

double StandardModel::getIntegrand_sigmaWithISR_down130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_down133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}


double StandardModel::getIntegrand_sigmaWithISR_down136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_sigmaWithISR_down207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}


//charm


double StandardModel::getIntegrand_sigmaWithISR_charm130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_charm133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_charm136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_charm207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}


//strange


double StandardModel::getIntegrand_sigmaWithISR_strange130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_strange133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_sigmaWithISR_strange207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}


//bottom


double StandardModel::getIntegrand_sigmaWithISR_bottom130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_bottom133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}    

double StandardModel::getIntegrand_sigmaWithISR_bottom136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_sigmaWithISR_bottom207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}





double StandardModel::Integrand_dsigmaBox_l(double cosTheta, const QCD::lepton l_flavor, const double s) const
{
    double ml = getLeptons(l_flavor).getMass();
    return ( myTwoFermionsLEP2->dsigma_l_box(l_flavor, ml, s, cosTheta, Mw(), Gamma_Z()) );
}       
    
double StandardModel::getIntegrand_dsigmaBox_mu130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_dsigmaBox_mu207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}





double StandardModel::getIntegrand_dsigmaBox_tau130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_dsigmaBox_tau207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}






double StandardModel::Integrand_dsigmaBox_q(double cosTheta, const QCD::quark q_flavor, const  double s) const
{
    double mq = m_q(q_flavor, sqrt(s)); 
    return ( myTwoFermionsLEP2->dsigma_q_box(q_flavor, mq, s, cosTheta, Mw(), Gamma_Z()) );
}       




//up


double StandardModel::getIntegrand_dsigmaBox_up130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}    

double StandardModel::getIntegrand_dsigmaBox_up133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

double StandardModel::getIntegrand_dsigmaBox_up207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}


//down

double StandardModel::getIntegrand_dsigmaBox_down130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}    

double StandardModel::getIntegrand_dsigmaBox_down133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

double StandardModel::getIntegrand_dsigmaBox_down207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}


//charm


double StandardModel::getIntegrand_dsigmaBox_charm130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}    

double StandardModel::getIntegrand_dsigmaBox_charm133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}    

double StandardModel::getIntegrand_dsigmaBox_charm136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_dsigmaBox_charm207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}


//strange


double StandardModel::getIntegrand_dsigmaBox_strange130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}    

double StandardModel::getIntegrand_dsigmaBox_strange133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}



double StandardModel::getIntegrand_dsigmaBox_strange172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

double StandardModel::getIntegrand_dsigmaBox_strange207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}


//bottom


double StandardModel::getIntegrand_dsigmaBox_bottom130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}    

double StandardModel::getIntegrand_dsigmaBox_bottom133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}    

double StandardModel::getIntegrand_dsigmaBox_bottom136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_dsigmaBox_bottom207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}












    
double StandardModel::Integrand_AFBnumeratorWithISR_l(double x, const QCD::lepton l_flavor, const  double s) const
{
    double sprime = (1.0 - x)*s;
    double Ncf = 1.0;
    double ml = getLeptons(l_flavor).getMass();
    double G3prime = myTwoFermionsLEP2->G_3prime_l(l_flavor, ml, sprime, Mw(), Gamma_Z(),flagLEP2[Weak]);
    double H = myTwoFermionsLEP2->H_ISR_FB(x, s);

    return ( M_PI*ale*ale*Ncf*H*G3prime/sprime );
}
    
    
    double StandardModel::getIntegrand_AFBnumeratorWithISR_mu130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_mu207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}
    

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}    

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_tau207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

    
    
double StandardModel::Integrand_AFBnumeratorWithISR_q(double x, const QCD::quark q_flavor, const  double s) const
{
    double sprime = (1.0 - x)*s;
    double Ncf = 3.0;
    double mq = m_q(q_flavor, sqrt(s)); 
    double G3prime = myTwoFermionsLEP2->G_3prime_q(q_flavor, mq, sprime, Mw(), Gamma_Z(),flagLEP2[Weak]);
    double H = myTwoFermionsLEP2->H_ISR_FB(x, s);
        
    if (flagLEP2[QCDFSR])
        G3prime *= myTwoFermionsLEP2->QCD_FSR_forAFB(q_flavor, mq, sprime);
        
    return ( M_PI*ale*ale*Ncf*H*G3prime/sprime );
}


double StandardModel::getIntegrand_AFBnumeratorWithISR_charm133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}    

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_charm207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}    
    
    


double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}    

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}
/* END: REMOVE FROM THE PACKAGE */