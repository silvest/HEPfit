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
    "delMw", "delSin2th_l", "delSin2th_q", "delSin2th_b", "delGammaZ", "delsigma0H", "delR0l", "delR0c", "delR0b", "delGammaWlv", "delGammaWqq",
    "mneutrino_1", "mneutrino_2", "mneutrino_3", "melectron", "mmu", "mtau", "muw"
};

const double StandardModel::GeVminus2_to_nb = 389379.338;
const double StandardModel::Mw_error = 0.00001; /* 0.01 MeV */

StandardModel::StandardModel()
: QCD(), Yu(3, 3, 0.), Yd(3, 3, 0.), Yn(3, 3, 0.),
SMM(*this), SMFlavour(*this), Ye(3, 3, 0.)
{
    setModelName("StandardModel");
    requireCKM = false;
    requireYe = false;
    requireYn = false;

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
    ModelParamMap.insert(std::make_pair("delGammaWlv", std::cref(delGammaWlv)));
    ModelParamMap.insert(std::make_pair("delGammaWqq", std::cref(delGammaWqq)));    
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

    SMSuccess = true;
    /* Set the CKM and PMNS matrices if not already set in the derived classes */
    if(requireCKM)
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
    else if (name.compare("delGammaWlv") == 0)
        delGammaWlv = value;
    else if (name.compare("delGammaWqq") == 0)
        delGammaWqq = value;    
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





void StandardModel::computeYukawas()
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
    // 21 parameters in StandardModel
    // GF, ale, dAle5Mz, mHl,
    // mneutrino_1, mneutrino_2, mneutrino_3, melectron, mmu, mtau,
    // delMw, delSin2th_l, delSin2th_q, delSin2th_b, delGammaZ, delsigma0H, delR0l, delR0c, delR0b, delGammaWlv, delGammaWqq,
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
        delMw, delSin2th_l, delSin2th_q, delSin2th_b, delGammaZ, delsigma0H, delR0l, delR0c, delR0b, delGammaWlv, delGammaWqq,
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

const double StandardModel::ale_OS(const double mu, orders order) const
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

const double StandardModel::Beta_s(int nm, unsigned int nf) const
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

const double StandardModel::Beta_e(int nm, unsigned int nf) const
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

const double StandardModel::AlsE(double mu, orders order, bool Nf_thr) const
{
    switch (order)
    {
        case FULLNNNLO:
            realorder = order;
            return (AlsByOrder(mu, LO, Nf_thr) + AlsByOrder(mu, NLO, Nf_thr) + AlsByOrder(mu, NNLO, Nf_thr) + AlsEByOrder(mu, NNNLO, Nf_thr));
        default:
            throw std::runtime_error("StandardModel::AlsE(): " + orderToString(order) + " is not implemented.");
    }
}

const double StandardModel::AlsEByOrder(double mu, orders order, bool Nf_thr) const
{
    int i, nfAls = (int) Nf(Mz), nfmu = Nf_thr ? (int) Nf(mu) : nfAls;
    double als, alstmp, mutmp;
    orders fullord;
    
    for (i = 0; i < CacheSize; ++i)
        if ((mu == als_cache[0][i]) && ((double) order == als_cache[1][i]) &&
                (AlsMz == als_cache[2][i]) && (Mz == als_cache[3][i]) &&
                (mut == als_cache[4][i]) && (mub == als_cache[5][i]) &&
                (muc == als_cache[6][i]) && (double) true == als_cache[7][i]
                && (double) Nf_thr == als_cache[8][i] && alphaMz() == als_cache[9][i])
            return als_cache[10][i];

    switch (order)
    {
        case NNNLO:
            if (nfAls == nfmu) 
                als = AlsEWithInit(mu, AlsMz, Mz, nfAls, order);
            fullord = FullOrder(order);
            if (nfAls > nfmu) {
                mutmp = BelowTh(Mz);
                alstmp = AlsEWithInit(mutmp, AlsMz, Mz, nfAls, realorder);
                alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(nfAls), alstmp, nfAls, fullord)); // WARNING: QED threshold corrections not implemented yet
                for (i = nfAls - 1; i > nfmu; i--) {
                    mutmp = BelowTh(mutmp - MEPS);
                    alstmp = AlsEWithInit(mutmp, alstmp, AboveTh(mutmp) - MEPS, i, realorder);
                    alstmp *= (1. - NfThresholdCorrections(mutmp, MassOfNf(i), alstmp, i, fullord)); // WARNING: QED threshold corrections not implemented yet
                }
                als = AlsEWithInit(mu, alstmp, AboveTh(mu) - MEPS, nfmu, order);
            }

            if (nfAls < nfmu) {
                mutmp = AboveTh(Mz) - MEPS;
                alstmp = AlsEWithInit(mutmp, AlsMz, Mz, nfAls, realorder);
                alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(nfAls + 1), alstmp, nfAls + 1, fullord)); // WARNING: QED threshold corrections not implemented yet
                for (i = nfAls + 1; i < nfmu; i++) {
                    mutmp = AboveTh(mutmp) - MEPS;
                    alstmp = AlsEWithInit(mutmp, alstmp, BelowTh(mutmp) + MEPS, i, realorder); 
                    alstmp *= (1. + NfThresholdCorrections(mutmp, MassOfNf(i + 1), alstmp, i + 1, fullord)); // WARNING: QED threshold corrections not implemented yet
                }
                als = AlsEWithInit(mu, alstmp, BelowTh(mu) + MEPS, nfmu, order);
            }

            CacheShift(als_cache, 11);
            als_cache[0][0] = mu;
            als_cache[1][0] = (double) order;       
            als_cache[2][0] = AlsMz;
            als_cache[3][0] = Mz;
            als_cache[4][0] = mut;
            als_cache[5][0] = mub;
            als_cache[6][0] = muc;
            als_cache[7][0] = (double) true;
            als_cache[8][0] = (double) Nf_thr;
            als_cache[9][0] = alphaMz();
            als_cache[10][0] = als;

             return als;
        default:
            throw std::runtime_error("StandardModel::AlsEByOrder(): " + orderToString(order) + " is not implemented.");
    }
}

const double StandardModel::AlsEWithInit(double mu, double alsi, double mu_i, const int nf_i, orders order) const
{
    double nf = (double) nf_i, alei = Ale(mu_i, FULLNLO); // CHANGE ME!
    double b00s = Beta_s(00, nf), b00e = Beta_e(00, nf);
    double v = 1. + b00s * alsi / 2. / M_PI * log(mu / mu_i);
    double ve = 1. - b00e * alei / 2. / M_PI * log(mu / mu_i);
    double logv = log(v), logve = log(ve);
    double rho = 1. / (1. + b00e * alei / b00s / alsi);
    double als = AlsWithInit(mu, alsi, mu_i, nf, order);
    double b01s = Beta_s(01,nf), b01s00e = b01s / b00e;

        switch (order)
        {
            case NNNLO:
                als += alsi * alsi * alei / 4. / 4. / M_PI / M_PI / v / v / ve * (Beta_s(02, nf) / b00e *
                        (ve - 1.) + Beta_s(11, nf) / b00s * rho * ve * (logve - logv) + b01s00e * Beta_e(10, nf) /
                        b00e * (logve - ve + 1.) + b01s * Beta_s(10, nf) / b00s / b00s * rho * logv +
                        b01s00e * Beta_e(01, nf) / b00s * (rho * ve * (logv - logve) - logv));
                break;
            case FULLNNNLO:
                return (AlsWithInit(mu, alsi, mu_i, nf_i, LO) + AlsWithInit(mu, alsi, mu_i, nf_i, NLO)+ AlsWithInit(mu, alsi, mu_i, nf_i, NNLO) + AlsEWithInit(mu, alsi, mu_i, nf_i, NNNLO));
            default:
                throw std::runtime_error("StandardModel::AlsEWithInit(): " + orderToString(order) + " is not implemented.");
        }

    return (als);
}

const double StandardModel::Ale(const double mu, orders order, bool Nf_thr) const
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

const double StandardModel::AleWithInit(double mu, double alei, double mu_i, orders order) const
{
    if (fabs(mu - mu_i) < MEPS) return(alei);

    double nf = Nf(mu), alsi = (mu_i == Mz ? AlsMz : Als(mu_i, FULLNNNLO, true, true));
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

const double StandardModel::DeltaAlphaLepton(const double s) const
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

const double StandardModel::DeltaAlphaL5q() const
{
    double Mz2 = Mz*Mz;
    return (DeltaAlphaLepton(Mz2) + dAl5hMz);
}

const double StandardModel::DeltaAlphaTop(const double s) const
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

const double StandardModel::DeltaAlpha() const
{
    if (FlagCacheInStandardModel)
        if (useDeltaAlpha_cache)
            return DeltaAlpha_cache;

    double Mz2 = Mz*Mz;
    DeltaAlpha_cache = DeltaAlphaL5q() + DeltaAlphaTop(Mz2);
    useDeltaAlpha_cache = true;
    return DeltaAlpha_cache;
}

const double StandardModel::alphaMz() const
{
    return (ale / (1.0 - DeltaAlpha()));
//    return(1./127.918); // FOR HEFFDF1 TEST: VALUE IN hep-ph/0512066
//    return(1./127.955); // FOR HEFFDF1 TEST: VALUE IN 2007.04191
}

const double StandardModel::Alstilde5(const double mu) const
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

const double StandardModel::v() const
{
    return ( 1. / sqrt(sqrt(2.) * GF));
}


///////////////////////////////////////////////////////////////////////////

const double StandardModel::Mw_tree() const
{
    if (FlagMWinput){
        return Mw_inp;
    } else
        return ( Mz / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0 * M_PI * ale / sqrt(2.0) / GF / Mz / Mz)));
}

const double StandardModel::s02() const
{
    double tmp = 1.0 - 4.0 * M_PI * alphaMz() / sqrt(2.0) / GF / Mz / Mz;
    if (tmp < 0.0)
        throw std::runtime_error("Error in s02()");

    return ( (1.0 - sqrt(tmp)) / 2.0);
}

const double StandardModel::c02() const
{
    return ( 1.0 - s02());
}

const double StandardModel::Mw() const
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

const double StandardModel::Dalpha5hMz() const
{
    if (FlagMWinput){
        return (myApproximateFormulae->dAlpha5hMw());
    } else
        return dAle5Mz;
}

const double StandardModel::cW2(double Mw_i) const
{
    return ( Mw_i * Mw_i / Mz / Mz);
}

const double StandardModel::cW2() const
{
    return ( cW2(Mw()));
//    return (1.0 - 0.2312); // FOR HEFFDF1 TEST
}

const double StandardModel::sW2(double Mw_i) const
{
    return ( 1.0 - cW2(Mw_i));
}

const double StandardModel::sW2() const
{
    return ( 1.0 - cW2());
}

const double StandardModel::sW2_MSbar_Approx() const 
{
    //double rho_t= 3. * getGF() * getMtpole() * getMtpole() / (8. * sqrt(2.) * M_PI * M_PI ); 
    return ( sW2()*1.0351 ); //PDG 22 electroweak review eq. (10.19)
}
 
const double StandardModel::sW2_ND() const 
{
    double d = 1. / 3. * (1. / sW2_MSbar_Approx() - 8. / 3.) * 
               ( (1 + getAlsMz()/M_PI)*log(getMtpole()/getMz()) - 15.*getAlsMz()/(8.*M_PI) ); 
    
    return sW2_MSbar_Approx()*(1. + Ale(getMz(),FULLNLO)*d/M_PI);

}

const double StandardModel::DeltaR() const
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

void StandardModel::ComputeDeltaRho(double Mw_i,
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

const double StandardModel::MwbarFromMw(const double Mw) const
{
    double AlsMw = Als(Mw, FULLNLO);
    double Gw_SM = 3.0 * GF * pow(Mw, 3.0) / 2.0 / sqrt(2.0) / M_PI
            * (1.0 + 2.0 * AlsMw / 3.0 / M_PI);

    return ( Mw - Gw_SM * Gw_SM / 2.0 / Mw);
}

const double StandardModel::MwFromMwbar(const double Mwbar) const
{
    double AlsMw = Als(Mwbar, FULLNNLO);
    double Gw_SM = 3.0 * GF * pow(Mwbar, 3.0) / 2.0 / sqrt(2.0) / M_PI
            * (1.0 + 2.0 * AlsMw / 3.0 / M_PI);

    return (Mwbar + Gw_SM * Gw_SM / 2.0 / Mwbar);
}

const double StandardModel::DeltaRbar() const
{
    double Mwbar_SM = MwbarFromMw(Mw());
    double sW2bar = 1.0 - Mwbar_SM * Mwbar_SM / Mzbar() / Mzbar();
    double tmp = sqrt(2.0) * GF * sW2bar * Mwbar_SM * Mwbar_SM / M_PI / ale;

    return (tmp - 1.0);
}


////////////////////////////////////////////////////////////////////////

const double StandardModel::rho_GammaW(const Particle fi, const Particle fj) const
{
    double rhoW = 0.0;
    if (flag_order[EW1])
        rhoW = myOneLoopEW->rho_GammaW(fi, fj, Mw());
    return rhoW;
}

const double StandardModel::GammaW(const Particle fi, const Particle fj) const
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
        return ( V.abs2() * G0 * rho_GammaW(fi, fj) * ( 1.0 + delGammaWlv ) );
    else {
        double AlsMw = AlsWithInit(Mw(), AlsMz, Mz, 5, FULLNLO);
        return ( 3.0 * V.abs2() * G0 * rho_GammaW(fi, fj) * (1.0 + AlsMw / M_PI) * ( 1.0 + delGammaWqq ) );
    }
}

const double StandardModel::GammaW() const
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


const double StandardModel::BrW(const Particle fi, const Particle fj) const
{
    double GammW = GammaW();
    double GammWij = GammaW(fi, fj);

    return GammWij/GammW;
}


const double StandardModel::RWlilj(const Particle li, const Particle lj) const
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

const double StandardModel::Ruc() const       //AG:added
{
    return 0.5 * ( R0_f(quarks[UP]) + R0_f(quarks[CHARM]) );
}

const double StandardModel::RWc() const
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

const double StandardModel::A_f(const Particle f) const
{
    double Re_kappa = kappaZ_f(f).real();
    double Re_gV_over_gA = 1.0 - 4.0 * fabs(f.getCharge()) * Re_kappa * sW2();
    return ( 2.0 * Re_gV_over_gA / (1.0 + pow(Re_gV_over_gA, 2.0)));
}

const double StandardModel::AFB(const Particle f) const
{
    return (3.0 / 4.0 * A_f(leptons[ELECTRON]) * A_f(f));
}

const double StandardModel::sin2thetaEff(const Particle f) const
{
    double Re_kappa = kappaZ_f(f).real();
    return ( Re_kappa * sW2());
}

const double StandardModel::GammaZ(const Particle f) const
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

const double StandardModel::Gamma_inv() const
{
    return ( GammaZ(leptons[NEUTRINO_1]) + GammaZ(leptons[NEUTRINO_2])
            + GammaZ(leptons[NEUTRINO_3]));
}

const double StandardModel::Gamma_had() const
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

const double StandardModel::Gamma_Z() const
{
    if (!IsFlagNoApproximateGammaZ()){
            
        /* SM contribution with the approximate formula */
        return myApproximateFormulae->X_full("GammaZ");

    } else {
        return ( GammaZ(leptons[ELECTRON]) + GammaZ(leptons[MU]) + GammaZ(leptons[TAU])
            + Gamma_inv() + Gamma_had());
    }
}


const double StandardModel::RZlilj(const Particle li, const Particle lj) const
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


const double StandardModel::sigma0_had() const
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

const double StandardModel::R0_f(const Particle f) const
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

const double StandardModel::R_inv() const
{
    return (Gamma_inv() / GammaZ(leptons[ELECTRON]));

}

const double StandardModel::N_nu() const
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

const gslpp::complex StandardModel::gV_f(const Particle f) const
{
    return ( gA_f(f)
            *(1.0 - 4.0 * fabs(f.getCharge())*(kappaZ_f(f)) * sW2()));
}

const gslpp::complex StandardModel::gA_f(const Particle f) const
{
    return ( sqrt(rhoZ_f(f)) * f.getIsospin());
}

const gslpp::complex StandardModel::rhoZ_f(const Particle f) const
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

const gslpp::complex StandardModel::kappaZ_f(const Particle f) const
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

const gslpp::complex StandardModel::deltaRhoZ_f(const Particle f) const
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

const gslpp::complex StandardModel::deltaKappaZ_f(const Particle f) const
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

const double StandardModel::epsilon1() const
{
    double rhoZe = rhoZ_f(leptons[ELECTRON]).real();
    double DeltaRhoPrime = 2.0 * (sqrt(rhoZe) - 1.0);

    return DeltaRhoPrime;
}

const double StandardModel::epsilon2() const
{
    double rhoZe = rhoZ_f(leptons[ELECTRON]).real();
    double sin2thetaEff = kappaZ_f(leptons[ELECTRON]).real() * sW2();
    double DeltaRhoPrime = 2.0 * (sqrt(rhoZe) - 1.0);
    double DeltaKappaPrime = sin2thetaEff / s02() - 1.0;
    double DeltaRW = 1.0 - M_PI * alphaMz() / (sqrt(2.0) * GF * Mz * Mz * sW2() * cW2());

    return ( c02() * DeltaRhoPrime + s02() * DeltaRW / (c02() - s02())
            - 2.0 * s02() * DeltaKappaPrime);
}

const double StandardModel::epsilon3() const
{
    double rhoZe = rhoZ_f(leptons[ELECTRON]).real();
    double sin2thetaEff = kappaZ_f(leptons[ELECTRON]).real() * sW2();
    double DeltaRhoPrime = 2.0 * (sqrt(rhoZe) - 1.0);
    double DeltaKappaPrime = sin2thetaEff / s02() - 1.0;

    return ( c02() * DeltaRhoPrime + (c02() - s02()) * DeltaKappaPrime);
}

const double StandardModel::epsilonb() const
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
        const double DeltaRbar_rem, bool bool_Zbb) const
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

////////////////////////////////////////////////////////////////////////     
// EW low-energy observables: Parity violation


//    The anomalous magnetic moment of the muon a_mu=(g_mu-2)/2

const double StandardModel::amuon() const
{
      
//      output
      double amu;
      
//      -----------------------------------------------------------------
//      qed contributions
      double amuqed,alfa0pi;
      
//      ew contributions
      double amuew,amuew1,amuew2b,amuew2f,amuew2,amuew3,cft,cf,corr1amuew2, corr2amuew2,corrwaamuew2,al,aq,b1; //,b2;

//      qcd contributions
      double amuhad,amuhhovp,amuhholbl,amuhho,amuhlo;
      
//      -----------------------------------------------------------------
//     numerical constants
      const double sn2=0.2604341;

//      -----------------------------------------------------------------
//     SM parameters

//     light quark masses. constituent masses
      const double umass=0.3;
      const double dmass=0.3;
      const double smass=0.5;
      
      const double mum=leptons[MU].getMass(),taum=leptons[TAU].getMass();
      const double cqm=quarks[CHARM].getMass(),bqm=quarks[BOTTOM].getMass();

//     all fermion masses (constituent masses for u,d,s. for the other from model)
      double fermmass[9]={leptons[ELECTRON].getMass(),mum,taum,
            dmass,umass,
            smass,cqm,
            bqm,mtpole};
     
//      w mass and on-shell weak angle
      double MwSM, s2;
      
//      running of alfa_qed and dummy variable
      double aqed;

//     for the 2-loop bosonic corrections
      double a2l[4]={0.,0.,0.,0.},b2l[4]={0.,0.,0.,0.},sw2l[4]={0.,0.,0.,0.};

//     for the 2-loop corrections from the renormalization of weak angle
      double c2lren[6]={0.,0.,0.,0.,0.,0.};
      
//      w mass
      MwSM=Mw();
      
      s2=1.0 - MwSM*MwSM/Mz/Mz;
      
//------------------------------------------------------------------
//      qed contribution to amu (arxiv: hep-ph/0606174)
      alfa0pi=ale/M_PI;
      
      amuqed=alfa0pi*(0.5+alfa0pi*(0.765857410+alfa0pi*(24.05050964+
     + alfa0pi*(130.8055+663.0*alfa0pi))));

//-----------------------------------------------------------------
//      one-loop ew correction(phys.rev.lett. 76,3267 (1996))

      amuew1=5.0*GF*mum*mum/(24.0*sqrt(2.0)*M_PI*M_PI)*(1.0+
     + 0.2*(1.0-4.0*s2)*(1.0-4.0*s2));

//-----------------------------------------------------------------
//      two-loop computation

//      these depend on aqed and since we are going to include also three-loop
//      effects we need to include in the two-loop results the running of aqed at
//      1-loop up to the scale mum
//-----------------------------------------------------------------
//      running of alpha em down to mu mass (1-loop)
            
      aqed = 1.0/ale + 2.0 * log(fermmass[0]/mum)/3.0/M_PI;
      
      aqed = 1.0/aqed;

//-----------------------------------------------------------------
//      two-loop ew bosonic correction(phys.rev.lett. 76,3267 (1996))

//      previous definitions
      a2l[0]=19.0/36.0-99.0*sn2/8.0-1.0*2.0*log(mHl/MwSM)/24.0;
      
      b2l[0]=155.0/192.0+3.0*M_PI*M_PI/8.0-9.0*sn2/8.0+3.0*2.0*pow(log(mHl/MwSM),2)/2.0-21.0*2.0*log(mHl/MwSM)/16.0;

      sw2l[0]=1.0/s2;

      a2l[1]=-859.0/18.0+11.0*M_PI/sqrt(3.0)+20.0*M_PI*M_PI/9.0+ 393.0*sn2/8.0-65.0*2.0*log(MwSM/mum)/9.0+ 31.0*2.0*log(mHl/MwSM)/72.0;

      b2l[1]=433.0/36.0+5.0*M_PI*M_PI/24.0-51.0*sn2/8.0+ 3.0*4.0*pow(log(mHl/MwSM),2)/8.0+9.0*2.0*log(mHl/MwSM)/4.0;

      sw2l[1]=1.0;

      a2l[2]=165169.0/1080.0-385.0*M_PI/(6.0*sqrt(3.0))-29.0*M_PI*M_PI/6.0+ 33.0*sn2/8.0+92.0*2.0*log(MwSM/mum)/9.0- 133.0*2.0*log(mHl/MwSM)/72.0;

      b2l[2]=-431.0/144.0+3.0*M_PI*M_PI/8.0+315.0*sn2/8.0+ 3.0*4.0*pow(log(mHl/MwSM),2)/2.0-11.0*2.0*log(mHl/MwSM)/8.0;

      sw2l[2]=s2;

      a2l[3]=-195965.0/864.0+265.0*M_PI/(3.0*sqrt(3.0))+163.0*M_PI*M_PI/18.0+ 223.0*sn2/12.0-184.0*2.0*log(MwSM/mum)/9.0- 5.0*2.0*log(mHl/MwSM)/8.0;

      b2l[3]=433.0/216.0+13.0*M_PI*M_PI/24.0+349.0*sn2/24.0+ 21.0*4.0*pow(log(mHl/MwSM),2)/8.0-49.0*2.0*log(mHl/MwSM)/12.0;

      sw2l[3]=s2*s2;

//      computation

      amuew2b=0.0;

      for (int i = 0; i < 4; ++i) {
            amuew2b=amuew2b+a2l[i]*sw2l[i]+(MwSM*MwSM/mHl/mHl)*b2l[i]*sw2l[i];
      }

//      the contribution with the running of aqed up to the mu scale
      amuew2b=mum*mum*aqed*GF*amuew2b/(8.0*sqrt(2.0)*M_PI*M_PI*M_PI);

//-----------------------------------------------------------------
//      two-loop ew fermionic correction(phys.rev.d 52,r2619(1995)

//      contribution from higgs boson diagram
      if (mHl < (mtpole-10.0)) {
            cft=-104.0/45.0-16.0*2.0*log(mtpole/mHl)/15.0;
      } else if (mHl > (mtpole+10)) {
            cft=-(mtpole*mtpole/mHl/mHl)*(24.0/5.0+8.0*M_PI*M_PI/15.0+
                                     + 8.0/5.0*pow(2.0*log(mHl/mtpole)-1.0,2));
      } else {
            cft=-(32.0/5.0)*(1.0-9.0*sn2/4.0);
      }

      cf=pow((umass*cqm*Mz),(4.0/3.0));

      cf=cf/(pow((dmass*smass*bqm),(1.0/3.0))*mum*mum*taum);

      cf=-18.0*log(cf)/5.0-3.0*mtpole*mtpole/(16.0*s2*MwSM*MwSM)- 3.0*2.0*log(mtpole/MwSM)/(10.0*s2)- 8.0*2.0*log(mtpole/Mz)/5.0-41.0/5.0-7.0/(10.0*s2)+ 8.0*M_PI*M_PI/15.0+cft;

//      the contribution with the running of aqed up to the mu scale
      amuew2f=5.0*GF*mum*mum*cf*aqed/(24.0*sqrt(2.0)*M_PI*M_PI*M_PI);

//-----------------------------------------------------------------
//      corrections from hadronic loops (phys.rev.d 67,073006(2003))
//      i also include the running here even though in the previous reference seems that it is not included
//      first family (eqs. (60) and (61))
      corr1amuew2=-aqed*GF*mum*mum/(8.0*M_PI*M_PI*M_PI*sqrt(2.0))*(8.41- log(pow(umass,8)/(pow(mum,6)*pow(dmass,2)))-17.0/2.0);
//      second family (eqs. (65) and (66))
      corr2amuew2=-aqed*GF*mum*mum/(8.0*M_PI*M_PI*M_PI*sqrt(2.0))*(17.1- log(pow(cqm,8)/(pow(mum,6)*pow(smass,2)))-47.0/6.0+8.0*M_PI*M_PI/9.0);

//-----------------------------------------------------------------
//      corrections from the renormalization of the weak mixing
//      terms prop. to (1-4s2) included in eq. (7) of phys.rev.d 67,073006(2003)
//      and neglected in the previous references
      
      corrwaamuew2=-43.0*31.0*(1.0-4.0*s2)*(1.0-4.0*s2)/(215.0*3.0)*log(Mz/mum);

      c2lren[0]=(72.0/135.0)*(-1.0+2.0*s2)*(1.0-4.0*s2); //leptons
      c2lren[1]=(72.0/135.0)*(-1.0+2.0*s2/3.0)*(1.0-4.0*s2); //d-quark
      c2lren[2]=-(144.0/135.0)*(1.0-4.0*s2/3.0)*(1.0-4.0*s2); //u-quark
      c2lren[3]=c2lren[1];//d-quark
      c2lren[4]=c2lren[2]; //u-quark
      c2lren[5]=c2lren[1]; //d-quark

      for (int i = 2; i < 8; ++i) {
            corrwaamuew2=corrwaamuew2+c2lren[i-2]*log(Mz/fermmass[i]);
      }

      corrwaamuew2=5*GF*mum*mum*aqed/(24.0*sqrt(2.0)*M_PI*M_PI*M_PI)*corrwaamuew2;

//      finally i also add the small correction to the eq.8
      corrwaamuew2=corrwaamuew2-0.2e-11;

//-----------------------------------------------------------------
//      total 2-loop ew contribution
      amuew2=amuew2b+amuew2f+corr1amuew2+corr2amuew2+corrwaamuew2;

//-----------------------------------------------------------------
//      three-loop ew correction(phys.rev.d 67,073006(2003)
      
      al=2789.0*log(Mz/mum)*log(Mz/mum)/90.0- 302.0*log(Mz/taum)*log(Mz/taum)/45.0+ 72.0*log(Mz/taum)*log(Mz/mum)/5.0;

      aq=-2662.0*log(Mz/bqm)*log(Mz/bqm)/1215.0+11216.0*log(Mz/cqm)*log(Mz/cqm)/1215.0+1964.0*log(Mz/umass)*log(Mz/umass)/405.0+24.0*log(Mz/bqm)*log(Mz/mum)/5.0-96.0*log(Mz/cqm)*log(Mz/mum)/5.0-48.0*log(Mz/umass)*log(Mz/mum)/5.0+32.0*log(Mz/bqm)*log(Mz/cqm)/405.0+32.0*log(Mz/bqm)*log(Mz/umass)/135.0;

      b1=-179.0/45.0*(log(Mz/bqm)*log(Mz/bqm)/3.0+log(Mz/taum)*log(Mz/taum)+4.0*log(Mz/cqm)*log(Mz/cqm)/3.0+2.0*log(Mz/umass)*log(Mz/umass)+2.0*log(Mz/mum)*log(Mz/mum))+2.0/5.0*(log(bqm/taum)*log(bqm/taum)+4.0/3.0*log(bqm/cqm)*log(bqm/cqm)+2.0*log(bqm/umass)*log(bqm/umass)+2.0*log(bqm/mum)*log(bqm/mum) )-8.0/5.0*(2.0*log(cqm/umass)*log(cqm/umass)+2.0*log(cqm/mum)*log(cqm/mum))+6.0/5.0*(4.0/3.0*log(taum/cqm)*log(taum/cqm)+2.0*log(taum/umass)*log(taum/umass)+2.0*log(taum/mum)*log(taum/mum))-8.0*log(umass/mum)*log(umass/mum)/5.0;

      // b2 is not used, as it can be absorved in the two loop part if alpha(m_mu) is used instead of alpha(Mz), as done above
      // b2=2.0/5.0*(2.0*log(Mz/mum)+2.0*log(Mz/umass)+4.0*log(Mz/cqm)/3.0+log(Mz/taum)+log(Mz/bqm)/3.0)*(215.0*log(Mz/mum)/9.0-4.0*log(Mz/umass)-8.0*log(Mz/cqm)+6.0*log(Mz/taum)+2.0*log(Mz/bqm));

//      the final correction(it is implied aqed at mum for the 2-loop
//      correction

      amuew3=amuew1*(ale*ale/M_PI/M_PI)*(al+aq+b1);

//-----------------------------------------------------------------
//      total ew correction

      amuew=amuew1+amuew2+amuew3;

//-----------------------------------------------------------------
//      hadronic contributions (arxiv: 0908.4300 & 1001.5401 [hep-ph])

//      leading order: vacuum polarization (arxiv: 0908.4300 [hep-ph])
      amuhlo=6955.e-11;

//      higher order: vacuum polarization
      amuhhovp=-97.9e-11;
      
//      higher order: light-by-light
      amuhholbl=105.e-11;

      amuhho=amuhhovp+amuhholbl;

//      total hadronic contribution

      amuhad=amuhlo+amuhho;

//-----------------------------------------------------------------
//      final value for the muon (g-2)/2

      amu=amuqed+amuew+amuhad;
      
//-----------------------------------------------------------------

      return amu;

}


//      The electron's weak charge

const double StandardModel::Qwemoller(const double q2, const double y) const
{
      //      Weak charge
      double Qwe;
      
      //      definitions
      double MwSM,f1,fy,f2,af2;
      const double mpion=134.9766e-3;
      
      //      -----------------------------------------------------------------

      double dalfos, dalfms, alfams;
      double rhoNC, kappa0, s2MSbar,c2MSbar;
      double xi;
      double leptk0,quarkk0;
      double elm=leptons[ELECTRON].getMass(), mum=leptons[MU].getMass(), taum=leptons[TAU].getMass();
      
      //      -----------------------------------------------------------------
      
      //      w mass
      MwSM=Mw();
      
      // xi factor
      xi=mHl*mHl/Mz/Mz;
      
      //      -----------------------------------------------------------------
      
      //      universal corrections
      //      ---------------------
            
      //      obtaining alfa(mz)_msbar from alfa(mz)_on-shell
      //      -----------------------------------------------
      
      //      on-shell value of delta alpha(mz)
      dalfos=1.0-ale/alphaMz();
      //      msbar value of delta alpha(mz) (formula from PDG, Erler & Langacker ew review)
      dalfms=dalfos+ale/M_PI*(100.0/27.0-1.0/6.0-7.0*2.0*log(Mz/MwSM)/4.0);
      //      msbar value of alfa(mz)
      alfams=ale/(1.0-dalfms);
            
      //      ms bar weinberg's angle from the effective leptonic angle
      //      (formula from PDG, Erler & Langacker ew review)
      //	---------------------------------------------------------
      s2MSbar=(myApproximateFormulae->sin2thetaEff_l_full())-0.00029;
      c2MSbar=1.0-s2MSbar;
      
      //      rho parameter (expansion in alfams)
      //      -------------
      
      rhoNC=1.0+alfams/(4.0*M_PI)*(3.0/(4.0*s2MSbar*s2MSbar)*log(c2MSbar)-7.0/(4.0*s2MSbar)+3.0*mtpole*mtpole/(4.0*s2MSbar*MwSM*MwSM) + 3.0*xi/(4.0*s2MSbar)*(log(c2MSbar/xi)/(c2MSbar-xi)+(1.0/c2MSbar)*log(xi)/(1.0-xi)));
           
      //      kappa at zero momentum (expansion in alfa)
      //      ----------------------
      
      //      lepton contribution to kappa0
      leptk0=((-0.5)*(-1)-2.0*s2MSbar)*2.0*(log(elm/Mz)+log(mum/Mz)+log(taum/Mz))/3.0;
      
      //      quark contribution to kappa0 (updated from hep-ph/0302149)
      quarkk0=-6.802;

      kappa0=1.0-ale/(2.0*M_PI*s2MSbar)*(leptk0+quarkk0-(7.0*c2MSbar/2.0+1.0/12.0)*log(c2MSbar)+(7.0/9.0-s2MSbar/3.0));
      
      //      -----------------------------------------------------------------
      
      //      f1(y,q2) (expansion in alfa)
      //      --------
      
      //      f(y)
      fy=-2.0*log(y*(1.0-y))/3.0+1.0/pow((1.0-y+y*y),2)*(-2.0*(1.0-y)*(3.0-3.0*y+4.0*y*y*y- 3.0*y*y*y*y)*log(1.0-y)-2.0*y*(1.0+3.0*y-6.0*y*y+8.0*y*y*y-3.0*y*y*y*y)*log(y)+ (1.0-y)*(2.0-2.0*y-7.0*y*y+10.0*y*y*y-8.0*y*y*y*y+3.0*y*y*y*y*y)*log(1.0-y)*log(1.0-y)- y*(2.0-3.0*y-5.0*y*y+8.0*y*y*y-7.0*y*y*y*y+3.0*y*y*y*y*y)*log(y)*log(y)+ (2.0-4.0*y+11.0*y*y*y-13.0*y*y*y*y+9.0*y*y*y*y*y-3.0*y*y*y*y*y*y)*(M_PI*M_PI-2.0*log(1.0-y)*log(y)));
      
      f1=-ale/(4.0*M_PI)*(1.0-4.0*kappa0*s2MSbar)*(22.0*log(y*Mz*Mz/q2)/3.0+85.0/9.0+fy);
      
      //      note that i have used 1-4*kappa*s2MSbar instead of 1-4*s2MSbar or an average as suggested in the
      //      reference
      
      
      //      f2(y,q2) (expansion in alfa)
      //      --------
      //      (y=1/2 approximattion using a pion loop calculation)
      
      //      af2
      af2=sqrt(1.0+4.0*mpion*mpion/q2);
      f2=ale/(4.0*M_PI)*(af2*af2*af2/3.0*log((af2+1.0)/(af2-1.0))-2.0/9.0-2.0*af2*af2/3.0);
      
      
      //      electron's weak charge
      //      ----------------------
      Qwe=-rhoNC*(1.0-4.0*kappa0*s2MSbar+alfams/(4.0*M_PI*s2MSbar)+f1+f2- 3.0*alfams*(1.0-4.0*kappa0*s2MSbar)*(1.0+(1.0-4.0*kappa0*s2MSbar)*(1.0-4.0*kappa0*s2MSbar))/(32.0*M_PI*s2MSbar*c2MSbar));
      
      //      again, i have used 1-4*kappa*s2MSbar even in the loop contributions
      
      return Qwe;
}



//     The parity violating asymmetry in Moller scattering

const double StandardModel::alrmoller(const double q2, const double y) const
{
      //      functions and inputs
      double alrmoller;
      
      // which alfa is this? => alpha(0). is this ale?
      
      //      parity violation asymmetry
      //      --------------------------
      alrmoller=-GF*q2*(1.0-y)/(sqrt(2.0)*M_PI*ale*(1.0+pow(y,4)+pow(1.0-y,4)))*Qwemoller(q2,y);
      
      return alrmoller;
}



//    The computation of the proton and neutron weak charge: Qwp,Qwn

const double StandardModel::Qwp() const
{
      //      Definitions
      double qwproton;

      double MwSM,alfapi,asMw,dkappa5h,s2MSbar0,deltae,deltaep,boxpww,boxpzz,boxpaz;
      //      I choose as lambda m_rho (pdg rho(770)) --> caz=3/2
      const double lambda=775.49e-3;
      const double caz=1.5;
      
      //      lepton masses
      double mlept[3]={leptons[ELECTRON].getMass(),leptons[MU].getMass(),leptons[TAU].getMass()};
      
      //      -----------------------------------------------------------------
      double dalfos, dalfms, alfams;
      double rhoNC, s2MSbar,c2MSbar;
      double xi;
      double elm=leptons[ELECTRON].getMass();
      //      -----------------------------------------------------------------
      
      //      W mass
      MwSM=Mw();
      
      // xi factor
      xi=mHl*mHl/Mz/Mz;
      
      //      alfa/pi
      alfapi=ale/M_PI;
      
      //      alfa_s(Mw)
      asMw = Als(MwSM, FULLNLO);
      
      //      -----------------------------------------------------------------
      
      //      Universal corrections
      //      ---------------------
            
      //      Obtaining alfa(mz)_msbar from alfa(mz)_on-shell
      //      -----------------------------------------------
      
      //      on-shell value of delta alpha(mz)
      dalfos=1.0-ale/alphaMz();
      //      MSbar value of delta alpha(mz) (formula from PDG, Erler & Langacker ew review)
      dalfms=dalfos+ale/M_PI*(100.0/27.0-1.0/6.0-7.0*2.0*log(Mz/MwSM)/4.0);
      //      MSbar value of alfa(mz)
      alfams=ale/(1.0-dalfms);
            
      //      MS bar weinberg's angle from the effective leptonic angle
      //      (formula from PDG, Erler & Langacker ew review)
      //	---------------------------------------------------------
      s2MSbar=(myApproximateFormulae->sin2thetaEff_l_full())-0.00029;
      c2MSbar=1.0-s2MSbar;
      
      //      rho parameter (expansion in alfams)
      //      -------------
      
      rhoNC=1.0+alfams/(4.0*M_PI)*(3.0/(4.0*s2MSbar*s2MSbar)*log(c2MSbar)-7.0/(4.0*s2MSbar)+3.0*mtpole*mtpole/(4.0*s2MSbar*MwSM*MwSM) + 3.0*xi/(4.0*s2MSbar)*(log(c2MSbar/xi)/(c2MSbar-xi)+(1.0/c2MSbar)*log(xi)/(1.0-xi)));
      
      //      -----------------------------------------------------------------
      
      //      sin2w_ms(0) eq.14
      //      -----------------
      
      //      hadronic contribution
      dkappa5h=7.9e-3;
      
      s2MSbar0=0.0;
      
      for (int i = 0; i < 3; ++i) {
            s2MSbar0=s2MSbar0+2.0*log(Mz/mlept[i]);
      }
      
      s2MSbar0=s2MSbar+dkappa5h+alfapi*((s2MSbar0*(1.0+0.75*alfapi)+135.0*alfapi/32.0)*(1.0-4.0*s2MSbar)/12.0- (7.0*c2MSbar/4.0+1.0/24.0)*2.0*log(Mz/MwSM)+s2MSbar/6.0-7.0/18.0);
      
      //      -----------------------------------------------------------------
      
      //      external leg corrections
      
      deltae=-0.5*alfapi;
      
      deltaep=-alfapi/3.0*(1.0-4.0*s2MSbar)*(2.0*log(Mz/elm)+1.0/6.0);
      
      //      -----------------------------------------------------------------
      
      //      boxes
      //      -----
      
      boxpww=alfams*(2.0+5.0*(1.0-asMw/M_PI))/(4.0*M_PI*s2MSbar);
       
      //      pure zz and az boxes from prd 17 3055 app.a
      
      boxpzz=alfams*(9.0/4.0-14.0*s2MSbar+38.0*s2MSbar*s2MSbar-40.0*s2MSbar*s2MSbar*s2MSbar)*(1.0-AlsMz/M_PI)/(4.0*M_PI*s2MSbar*c2MSbar);
      
      boxpaz=5.0*alfams*(1.0-4.0*s2MSbar)*(2.0*log(Mz/lambda)+caz)/(2.0*M_PI);
      
      //      i assumme the same caz as in the proton enters for the neutron
      //      -----------------------------------------------------------------
      
      //      weak charges
      //      ------------
      
      qwproton=(rhoNC+deltae)*(1.0-4.0*s2MSbar0+deltaep)+boxpww+boxpzz+boxpaz;
      
      return qwproton;
      
}


const double StandardModel::Qwn() const
{
      //      Definitions
      double qwneutron;

      double MwSM,alfapi,asMw,dkappa5h,s2MSbar0,deltae,deltaep,boxnww,boxnzz,boxnaz;
      //      I choose as lambda m_rho (pdg rho(770)) --> caz=3/2
      const double lambda=775.49e-3;
      const double caz=1.5;
      
      //      lepton masses
      double mlept[3]={leptons[ELECTRON].getMass(),leptons[MU].getMass(),leptons[TAU].getMass()};
      
      //      -----------------------------------------------------------------
      double dalfos, dalfms, alfams;
      double rhoNC, s2MSbar,c2MSbar;
      double xi;
      double elm=leptons[ELECTRON].getMass();
      //      -----------------------------------------------------------------
      
      //      W mass
      MwSM=Mw();
      
      // xi factor
      xi=mHl*mHl/Mz/Mz;
      
      //      alfa/pi
      alfapi=ale/M_PI;
      
      //      alfa_s(Mw)
      asMw = Als(MwSM, FULLNLO);
      
      //      -----------------------------------------------------------------
      
      //      Universal corrections
      //      ---------------------
            
      //      Obtaining alfa(mz)_msbar from alfa(mz)_on-shell
      //      -----------------------------------------------
      
      //      on-shell value of delta alpha(mz)
      dalfos=1.0-ale/alphaMz();
      //      MSbar value of delta alpha(mz) (formula from PDG, Erler & Langacker ew review)
      dalfms=dalfos+ale/M_PI*(100.0/27.0-1.0/6.0-7.0*2.0*log(Mz/MwSM)/4.0);
      //      MSbar value of alfa(mz)
      alfams=ale/(1.0-dalfms);
            
      //      MS bar weinberg's angle from the effective leptonic angle
      //      (formula from PDG, Erler & Langacker ew review)
      //	---------------------------------------------------------
      s2MSbar=(myApproximateFormulae->sin2thetaEff_l_full())-0.00029;
      c2MSbar=1.0-s2MSbar;
      
      //      rho parameter (expansion in alfams)
      //      -------------
      
      rhoNC=1.0+alfams/(4.0*M_PI)*(3.0/(4.0*s2MSbar*s2MSbar)*log(c2MSbar)-7.0/(4.0*s2MSbar)+3.0*mtpole*mtpole/(4.0*s2MSbar*MwSM*MwSM) + 3.0*xi/(4.0*s2MSbar)*(log(c2MSbar/xi)/(c2MSbar-xi)+(1.0/c2MSbar)*log(xi)/(1.0-xi)));
      
      //      -----------------------------------------------------------------
      
      //      sin2w_ms(0) eq.14
      //      -----------------
      
      //      hadronic contribution
      dkappa5h=7.9e-3;
      
      s2MSbar0=0.0;
      
      for (int i = 0; i < 3; ++i) {
            s2MSbar0=s2MSbar0+2.0*log(Mz/mlept[i]);
      }
      
      s2MSbar0=s2MSbar+dkappa5h+alfapi*((s2MSbar0*(1.0+0.75*alfapi)+135.0*alfapi/32.0)*(1.0-4.0*s2MSbar)/12.0- (7.0*c2MSbar/4.0+1.0/24.0)*2.0*log(Mz/MwSM)+s2MSbar/6.0-7.0/18.0);
      
      //      -----------------------------------------------------------------
      
      //      external leg corrections
      
      deltae=-0.5*alfapi;
      
      deltaep=-alfapi/3.0*(1.0-4.0*s2MSbar)*(2.0*log(Mz/elm)+1.0/6.0);
      
      //      -----------------------------------------------------------------
      
      //      boxes
      //      -----
      
      boxnww=alfams*(-2.0+4.0*(1.0-asMw/M_PI))/(4.0*M_PI*s2MSbar);
      
      //      pure zz and az boxes from prd 17 3055 app.a
      
      boxnzz=alfams*(9.0/4.0-13.0*s2MSbar+34.0*s2MSbar*s2MSbar-32.0*s2MSbar*s2MSbar*s2MSbar)*(1.0-AlsMz/M_PI)/(4.0*M_PI*s2MSbar*c2MSbar);
      
      //      i assumme the same caz as in the proton enters for the neutron
      boxnaz=alfams*(4.0-16.0*s2MSbar)*(2.0*log(Mz/lambda)+caz)/(2.0*M_PI);
      
      //      -----------------------------------------------------------------
      
      //      weak charges
      //      ------------
      
      qwneutron=-(rhoNC+deltae)*(1.0+deltaep)+boxnww+boxnzz+boxnaz;
      
      return qwneutron;
      
}


    ////////////////////////////////////////////////////////////////////////     
    // EW low-energy observables: neutrino-scattering


const double StandardModel::gLnuN2() const
{    
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEgLnuN2Approx());

    } else {
        throw std::runtime_error("ERROR: StandardModel::gLnuN2, prediction implemented only via semianalytical approximate formula. Check flags!");
    }      
}


const double StandardModel::gRnuN2() const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEgRnuN2Approx());

    } else {
        throw std::runtime_error("ERROR: StandardModel::gRnuN2, prediction implemented only via semianalytical approximate formula. Check flags!");
    }   
}

const double StandardModel::ThetaLnuN() const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEThetaLnuNApprox());

    } else {
        throw std::runtime_error("ERROR: StandardModel::ThetaLnuN, prediction implemented only via semianalytical approximate formula. Check flags!");
    }   
}


const double StandardModel::ThetaRnuN() const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEThetaRnuNApprox());

    } else {
        throw std::runtime_error("ERROR: StandardModel::ThetaRnuN, prediction implemented only via semianalytical approximate formula. Check flags!");
    }   
}

const double StandardModel::gVnue() const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEgVnueApprox());

    } else {
        throw std::runtime_error("ERROR: StandardModel::gVnue, prediction implemented only via semianalytical approximate formula. Check flags!");
    }   
}

const double StandardModel::gAnue() const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEgAnueApprox());

    } else {
        throw std::runtime_error("ERROR: StandardModel::gAnue, prediction implemented only via semianalytical approximate formula. Check flags!");
    }   
}



////////////////////////////////////////////////////////////////////////     
// Lepton decays

// Muon decay

const double StandardModel::Gamma_muon() const
{
    double Gamma;
    double me, mmu, x, Fx, H1x, H2x, H3x, zeta3;
    double alpha, rEW;
    double pi2;
      
    me = leptons[ELECTRON].getMass();
    mmu = leptons[MU].getMass();
    pi2 = M_PI*M_PI;
      
    x = me*me/mmu/mmu;
    Fx = 1. - 8. * x + 8. * x*x*x - x*x*x*x -12. * x*x * log(x);
      
    H1x = 25./8. - pi2/2. - (9. + 4. *pi2 + 12. * log(x) )*x + 16. * pi2 * pow(x,3./2.);
      
    zeta3 = 1.2020569031595942;
      
    H2x= 156815./5184. - 518. * pi2/81. - 895. *zeta3/36. + 67.*pi2*pi2/720. + 53. *pi2*log(2.)/6. - 0.042 - (5./4.) * pi2*sqrt(x);
      
    H3x = -15.3;
      
    // alpha(m_mu)
    alpha = 1./ale - log(x)/3./M_PI; // + 1./6./M_PI;
    alpha = 1./alpha;
      
    // Rad. corrections
    rEW = 1. + H1x * alpha/M_PI + H2x * alpha*alpha/pi2 + H3x * alpha * alpha *alpha/pi2/M_PI;

    // Gamma: PDG formula
    Gamma = GF*GF*pow(mmu,5)*Fx*rEW/192./pow(M_PI,3);
                      
    return Gamma;
}


// Tau decays

// Leptonic decays

const double StandardModel::Gamma_tau_l_nunu(const Particle l) const
{
    double Gamma;
    double ml, mtau, x, Fx, H1x, H2x, H3x, zeta3;
    double alpha, rEW;
    double pi2;
      
    ml = l.getMass();
    mtau = leptons[TAU].getMass();
    pi2 = M_PI*M_PI;
      
    x = ml*ml/mtau/mtau;
    Fx = 1. - 8. * x + 8. * x*x*x - x*x*x*x -12. * x*x * log(x);
      
    H1x = 25./8. - pi2/2. - (9. + 4. *pi2 + 12. * log(x) )*x + 16. * pi2 * pow(x,3./2.);
      
    zeta3 = 1.2020569031595942;
      
    H2x= 156815./5184. - 518. * pi2/81. - 895. *zeta3/36. + 67.*pi2*pi2/720. + 53. *pi2*log(2.)/6. - 0.042 - (5./4.) * pi2*sqrt(x);
      
    H3x = -15.3;
      
    // alpha(m_tau)
    alpha = 1./133.29; // Improve
      
    // Rad. corrections
    rEW = 1. + H1x * alpha/M_PI + H2x * alpha*alpha/pi2 + H3x * alpha * alpha *alpha/pi2/M_PI;

    // Gamma: PDG formula
    Gamma = GF*GF*pow(mtau,5)*Fx*rEW/192./pow(M_PI,3);
                      
    return Gamma;
}


// Lepton universality tests

const double StandardModel::TauLFU_gmuge() const
{
    double g2LFU;
    
    double me, mmu, mtau, xe, Fxe, xmu, Fxmu;
      
    me = leptons[ELECTRON].getMass();
    mmu = leptons[MU].getMass();
    mtau = leptons[TAU].getMass();
    
    xe = me*me/mtau/mtau;
    Fxe = 1. - 8. * xe + 8. * xe*xe*xe - xe*xe*xe*xe -12. * xe*xe * log(xe);     
      
    xmu = mmu*mmu/mtau/mtau;
    Fxmu = 1. - 8. * xmu + 8. * xmu*xmu*xmu - xmu*xmu*xmu*xmu -12. * xmu*xmu * log(xmu); 
    
    g2LFU = (Gamma_tau_l_nunu(leptons[MU])/Gamma_tau_l_nunu(leptons[ELECTRON]));
    
    g2LFU = g2LFU * (Fxe/Fxmu);
    
    return sqrt(g2LFU);
}

const double StandardModel::TauLFU_gtaugmu() const
{
    double g2LFU;
    
    double me, mmu, mtau, xtau, Fxtau, xmu, Fxmu;
      
    me = leptons[ELECTRON].getMass();
    mmu = leptons[MU].getMass();
    mtau = leptons[TAU].getMass();
    
    xtau = me*me/mtau/mtau;
    Fxtau = 1. - 8. * xtau + 8. * xtau*xtau*xtau - xtau*xtau*xtau*xtau -12. * xtau*xtau * log(xtau);     
      
    xmu = me*me/mmu/mmu;
    Fxmu = 1. - 8. * xmu + 8. * xmu*xmu*xmu - xmu*xmu*xmu*xmu -12. * xmu*xmu * log(xmu); 
      
    g2LFU = (Gamma_tau_l_nunu(leptons[ELECTRON])/Gamma_muon());
    
    g2LFU = g2LFU * (pow(mmu,5)*Fxmu/pow(mtau,5)/Fxtau);
    
    return sqrt(g2LFU);
}

const double StandardModel::TauLFU_gtauge() const
{
    double g2LFU;
    
    double me, mmu, mtau, xtau, Fxtau, xmu, Fxmu;
      
    me = leptons[ELECTRON].getMass();
    mmu = leptons[MU].getMass();
    mtau = leptons[TAU].getMass();
    
    xtau = mmu*mmu/mtau/mtau;
    Fxtau = 1. - 8. * xtau + 8. * xtau*xtau*xtau - xtau*xtau*xtau*xtau -12. * xtau*xtau * log(xtau);     
      
    xmu = me*me/mmu/mmu;
    Fxmu = 1. - 8. * xmu + 8. * xmu*xmu*xmu - xmu*xmu*xmu*xmu -12. * xmu*xmu * log(xmu); 
      
    g2LFU = (Gamma_tau_l_nunu(leptons[MU])/Gamma_muon());
    
    g2LFU = g2LFU * (pow(mmu,5)*Fxmu/pow(mtau,5)/Fxtau);
    
    return sqrt(g2LFU);
}


const double StandardModel::TauLFU_gtaugmuPi() const
{
    // 1st approx. 
    
    return 1.0;
}

const double StandardModel::TauLFU_gtaugmuK() const
{
    // 1st approx. 
    
    return 1.0;
}


////////////////////////////////////////////////////////////////////////     
// Higgs processes
////////////////////////////////////////////////////////////////////////

//  Integrals

gslpp::complex StandardModel::f_triangle(const double tau) const {
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = asin(1.0 / sqrt(tau));
        return (tmp * tmp);
    } else {
        tmp = log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i();
        return (-0.25 * tmp * tmp);
    }
}

gslpp::complex StandardModel::g_triangle(const double tau) const {
    gslpp::complex tmp;
    if (tau >= 1.0) {
        tmp = sqrt(tau - 1.0) * asin(1.0 / sqrt(tau));
        return tmp;
    } else {
        tmp = sqrt(1.0 - tau) * (log((1.0 + sqrt(1.0 - tau)) / (1.0 - sqrt(1.0 - tau))) - M_PI * gslpp::complex::i());
        return 0.5 * tmp;
    }
}

gslpp::complex StandardModel::I_triangle_1(const double tau, const double lambda) const {
    gslpp::complex tmp;

    tmp = (tau * lambda * (f_triangle(tau) - f_triangle(lambda)) + 2.0 * tau * (g_triangle(tau) - g_triangle(lambda))) / (tau - lambda);

    tmp = tau * lambda * (1.0 + tmp) / (2.0 * (tau - lambda));

    return tmp;
}

gslpp::complex StandardModel::I_triangle_2(const double tau, const double lambda) const {
    gslpp::complex tmp;

    tmp = -0.5 * tau * lambda * (f_triangle(tau) - f_triangle(lambda)) / (tau - lambda);

    return tmp;
}

gslpp::complex StandardModel::AH_f(const double tau) const {
    return (2.0 * tau * (1.0 + (1.0 - tau) * f_triangle(tau)));
}

gslpp::complex StandardModel::AH_W(const double tau) const {
    return -(2.0 + 3.0 * tau + 3.0 * tau * (2.0 - tau) * f_triangle(tau));
}

gslpp::complex StandardModel::AHZga_f(const double tau, const double lambda) const {
    return I_triangle_1(tau, lambda) - I_triangle_2(tau, lambda);
}

gslpp::complex StandardModel::AHZga_W(const double tau, const double lambda) const {
    gslpp::complex tmp;

    double tan2w = sW2() / cW2();

    tmp = 4.0 * (3.0 - tan2w) * I_triangle_2(tau, lambda);

    tmp = tmp + ((1.0 + 2.0 / tau) * tan2w - (5.0 + 2.0 / tau)) * I_triangle_1(tau, lambda);

    return sqrt(cW2()) * tmp;
}

////////////////////////////////////////////////////////////////////////

const double StandardModel::SigmaeeZH(const double sqrt_s, const double Pe, const double Pp) const
{   
    double xsLH, xsRH;
    double gL,gR,lam,fact;
    double s = sqrt_s*sqrt_s;
    
    // From https://arxiv.org/pdf/hep-ph/9605437
    
    gL = -0.5 + sW2();
    
    gR = sW2();
    
    lam = (1.0-(mHl+Mz)*(mHl+Mz)/s)*(1.0-(mHl-Mz)*(mHl-Mz)/s);
    
    fact = (pow(GF*Mz*Mz,2.0)/96.0/M_PI/s) * sqrt(lam)*( lam + 12.0*Mz*Mz/s )/( 1.0 - Mz*Mz/s )/( 1.0 - Mz*Mz/s );
            
    xsLH = 32.0 * gL * gL * fact;
    xsRH = 32.0 * gR * gR * fact;

    return 0.25*( (1.0 - Pe)*(1.0 + Pp)*xsLH + (1.0 + Pe)*(1.0 - Pp)*xsRH );
}

const double StandardModel::SigmaeeHvv(const double sqrt_s, const double Pe, const double Pp) const
{   
    double xsLH=1.0, xsRH=0.0;

    return 0.25*( (1.0 - Pe)*(1.0 + Pp)*xsLH + (1.0 + Pe)*(1.0 - Pp)*xsRH ); 
}

const double StandardModel::SigmaeeHee(const double sqrt_s, const double Pe, const double Pp) const
{   
    double xsLH=0.0, xsRH=0.0;

    return 0.25*( (1.0 - Pe)*(1.0 + Pp)*xsLH + (1.0 + Pe)*(1.0 - Pp)*xsRH ); 
}

////////////////////////////////////////////////////////////////////////
// Higgs decay widths
////////////////////////////////////////////////////////////////////////

const double StandardModel::GammaHtogg() const
{   
    double gamma;
    double tau_t = 4.0 * pow(quarks[TOP].getMass(),2)/mHl/mHl; 
    double tau_b = 4.0 * pow(quarks[BOTTOM].getMass(),2)/mHl/mHl; 
    double tau_c = 4.0 * pow(quarks[CHARM].getMass(),2)/mHl/mHl; 
    double tau_s = 4.0 * pow(quarks[STRANGE].getMass(),2)/mHl/mHl; 
    double asMH,LH,Lt,nl,h0,h1,h2, h3,G0;
    
    //      alfa_s(MH)
    asMH = Als(mHl, FULLNLO);
    
    // NLO corrections ( See https://arxiv.org/pdf/0708.0916 and its REf. [25])
    // I only keep up to h3 in expr. (4), and use pole mass in tau factors for the moment
    nl = 5;
    LH = 0.; //  log(mu^2/MH^2) evaluated at mu=MH
    Lt = 2.0*log(mHl/(quarks[TOP].getMass()));
    
    h0 = (95./4.) + (11./2.)*LH + nl*(-7./6. - LH/3.);
    h1 = 5803./540. + 77.*LH/30. -14.*Lt/15. + nl * (-29./60. - 7. * LH / 45.);
    h2 = 1029839./189000. + 16973.*LH/12600. - 1543.*Lt/1575. + nl * ( - 89533./378000 - 1543.*LH/18900. );
    h3 = 9075763./2976750. + 1243*LH/1575. - 452.*Lt/575. + nl * ( - 3763./28350. -226. * LH / 4725. );
    G0 = GF * pow(mHl,3.0)/(36.*M_PI*sqrt(2.));
    
    gamma = asMH*asMH * (4.0 * GF /sqrt(2.0)) * (mHl*mHl*mHl /64.0/pow(M_PI,3.0)) * 
            ( AH_f(tau_t) + AH_f(tau_b) + AH_f(tau_c) + AH_f(tau_s) ).abs2()/4.0;
    
    gamma = gamma + G0 * (asMH/M_PI) * (asMH/M_PI) * (asMH/M_PI) * (h0 + h1/tau_t + h2/tau_t/tau_t + h3/tau_t/tau_t/tau_t );
    
    return gamma;
}

const double StandardModel::GammaHtoZZstar() const
{
    double x=Mz/mHl;
    double fx;
    double g2 = 4.0 * sqrt(2.0) * GF * Mz * Mz;
    double gamma;
    
    fx = -fabs(1.0-x*x)*( 47.0*x*x/2.0 - 13.0/2.0 +1.0/x/x ) + 
            3.0*( 1.0 - 6.0*x*x + 4.0*x*x*x*x )*fabs(log(x)) + 
            3.0*( 1.0 - 8.0*x*x + 20.0*x*x*x*x )*acos(( 3.0*x*x - 1.0 )/2.0/x/x/x)/sqrt( 4.0*x*x- 1.0);
    
    gamma = g2*g2 * mHl * fx * ( 7.0 - 40.0*sW2()/3.0 + 160.0 *sW2()*sW2()/9.0 ) / 2048.0 / pow(M_PI,3.0);
    
    return gamma;
}
    
const double StandardModel::GammaHtoWWstar() const
{
    double x=Mw()/mHl;
    double fx;
    double g2 = 4.0 * sqrt(2.0) * GF * pow(Mw(),2);
    double gamma;
    
    fx = -fabs(1.0-x*x)*( 47.0*x*x/2.0 - 13.0/2.0 +1.0/x/x ) + 
            3.0*( 1.0 - 6.0*x*x + 4.0*x*x*x*x )*fabs(log(x)) + 
            3.0*( 1.0 - 8.0*x*x + 20.0*x*x*x*x )*acos(( 3.0*x*x - 1.0 )/2.0/x/x/x)/sqrt( 4.0*x*x- 1.0);
    
    gamma = 3.0 * g2*g2 * mHl * fx / 512.0 / pow(M_PI,3.0);
    
    return gamma;
}

const double StandardModel::GammaHtoZga() const
{       
    double gamma;
    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double m_s = quarks[STRANGE].getMass();
    double m_tau = leptons[TAU].getMass();
    double m_mu = leptons[MU].getMass();

    double M_w_2 = pow(Mw(),2.0);

    double Qt = quarks[TOP].getCharge();
    double Qb = quarks[BOTTOM].getCharge();
    double Qc = quarks[CHARM].getCharge();
    double Qs = quarks[STRANGE].getCharge();
    double Qtau = leptons[TAU].getCharge();
    double Qmu = leptons[MU].getCharge();

    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double tau_s = 4.0 * m_s * m_s / mHl / mHl;
    double tau_tau = 4.0 * m_tau * m_tau / mHl / mHl;
    double tau_mu = 4.0 * m_mu * m_mu / mHl / mHl;
    double tau_W = 4.0 * M_w_2 / mHl / mHl;

    double lambda_t = 4.0 * m_t * m_t / Mz / Mz;
    double lambda_b = 4.0 * m_b * m_b / Mz / Mz;
    double lambda_c = 4.0 * m_c * m_c / Mz / Mz;
    double lambda_s = 4.0 * m_s * m_s / Mz / Mz;
    double lambda_tau = 4.0 * m_tau * m_tau / Mz / Mz;
    double lambda_mu = 4.0 * m_mu * m_mu / Mz / Mz;
    double lambda_W = 4.0 * M_w_2 / Mz / Mz;

    double sc = sqrt(sW2()*cW2());   
    double vSMt = (2.0 * (quarks[TOP].getIsospin()) - 4.0 * Qt * sW2())/sc;
    double vSMb = (2.0 * (quarks[BOTTOM].getIsospin()) - 4.0 * Qb * sW2())/sc;
    double vSMc = (2.0 * (quarks[CHARM].getIsospin()) - 4.0 * Qc * sW2())/sc;
    double vSMs = (2.0 * (quarks[STRANGE].getIsospin()) - 4.0 * Qs * sW2())/sc;
    double vSMtau = (2.0 * (leptons[TAU].getIsospin()) - 4.0 * Qtau * sW2())/sc;
    double vSMmu = (2.0 * (leptons[MU].getIsospin()) - 4.0 * Qmu * sW2())/sc;
    
    gslpp::complex MSM;

    MSM = (ale/4.0/M_PI) * ((3.0 * vSMt * Qt * AHZga_f(tau_t, lambda_t) +
            3.0 * vSMb * Qb * AHZga_f(tau_b, lambda_b) +
            3.0 * vSMc * Qc * AHZga_f(tau_c, lambda_c) +
            3.0 * vSMs * Qs * AHZga_f(tau_s, lambda_s) +
            vSMtau * Qtau * AHZga_f(tau_tau, lambda_tau) +
            vSMmu * Qmu * AHZga_f(tau_mu, lambda_mu)) +
            AHZga_W(tau_W, lambda_W)/sqrt(sW2()));

    gamma = (4.0*sqrt(2)*GF) * (MSM.abs2()) * pow(mHl*(1.0-Mz*Mz/mHl/mHl),3.0)/32.0/M_PI;

    return gamma;
}

const double StandardModel::GammaHtogaga() const
{   
    double gamma;
    
    double m_t = mtpole;
    double m_b = quarks[BOTTOM].getMass();
    double m_c = quarks[CHARM].getMass();
    double m_s = quarks[STRANGE].getMass();
    double m_tau = leptons[TAU].getMass();
    double m_mu = leptons[MU].getMass();

    double M_w_2 = pow(Mw(),2.0);

    double Qt = quarks[TOP].getCharge();
    double Qb = quarks[BOTTOM].getCharge();
    double Qc = quarks[CHARM].getCharge();
    double Qs = quarks[STRANGE].getCharge();
    double Qtau = leptons[TAU].getCharge();
    double Qmu = leptons[MU].getCharge();

    double tau_t = 4.0 * m_t * m_t / mHl / mHl;
    double tau_b = 4.0 * m_b * m_b / mHl / mHl;
    double tau_c = 4.0 * m_c * m_c / mHl / mHl;
    double tau_s = 4.0 * m_s * m_s / mHl / mHl;
    double tau_tau = 4.0 * m_tau * m_tau / mHl / mHl;
    double tau_mu = 4.0 * m_mu * m_mu / mHl / mHl;
    double tau_W = 4.0 * M_w_2 / mHl / mHl;

    gslpp::complex MSM;

    MSM = ale * (3.0 * Qt * Qt * AH_f(tau_t) +
            3.0 * Qb * Qb * AH_f(tau_b) +
            3.0 * Qc * Qc * AH_f(tau_c) +
            3.0 * Qs * Qs * AH_f(tau_s) +
            Qtau * Qtau * AH_f(tau_tau) +
            Qmu * Qmu * AH_f(tau_mu) +
            AH_W(tau_W));

    gamma = (4.0*GF/sqrt(2)) * (MSM.abs2()) * pow(mHl,3.0)/512.0/pow(M_PI,3);

    return gamma;
}
    
const double StandardModel::GammaHtomumu() const
{   
    double mf=leptons[MU].getMass();
    double beta=1.0-4.0*mf*mf/mHl/mHl;
    double Nc=1.0;
    double gamma;
    
    gamma = Nc * (4.0*GF/sqrt(2.0)) * (mf*mf/16.0/M_PI) * mHl * beta*beta*beta;
    
    return gamma;
}

const double StandardModel::GammaHtotautau() const
{   
    double mf=leptons[TAU].getMass();
    double beta=1.0-4.0*mf*mf/mHl/mHl;
    double Nc=1.0;
    double gamma;
    
    gamma = Nc * (4.0*GF/sqrt(2.0)) * (mf*mf/16.0/M_PI) * mHl * beta*beta*beta;
    
    return gamma;
}

const double StandardModel::GammaHtocc() const
{   
    double mf0=quarks[CHARM].getMass(), mf;
    double beta;
    double Nc=3.0;
    double gamma;
    double asMH,DeltaQCD,Deltamt,NF;
    
    //      alfa_s(MH)
    asMH = Als(mHl, FULLNLO);
    
    mf = Mrun(mHl, mf0, mf0, CHARM, FULLNLO);
    
    beta=1.0-4.0*mf*mf/mHl/mHl;
    
    NF=5;
    
    DeltaQCD = 1 + (asMH/M_PI) * ( 17.0/3.0 + (asMH/M_PI) * ( (35.94 - 1.36*NF) + (164.14 - 25.77*NF + 0.26*NF*NF)*(asMH/M_PI) ) );
    
    Deltamt = (asMH/M_PI) * (asMH/M_PI) * ( 1.57 + (4.0/3.0)*log(mHl/mtpole) + (4.0/9.0) * log(mf/mHl) * log(mf/mHl) );
    
    gamma = Nc * (4.0*GF/sqrt(2.0)) * (mf*mf/16.0/M_PI) * mHl * beta*beta*beta * (DeltaQCD + Deltamt);
    
    return gamma;
}
    
const double StandardModel::GammaHtoss() const
{   
    double mf0=quarks[STRANGE].getMass(), mf;
    double beta;
    double Nc=3.0;
    double gamma;
    double asMH,DeltaQCD,Deltamt,NF;
    
    //      alfa_s(MH)
    asMH = Als(mHl, FULLNLO);
    
    mf = Mrun(mHl, 2.0, mf0, STRANGE, FULLNLO);
    
    beta=1.0-4.0*mf*mf/mHl/mHl;
    
    NF=5;
    
    DeltaQCD = 1 + (asMH/M_PI) * ( 17.0/3.0 + (asMH/M_PI) * ( (35.94 - 1.36*NF) + (164.14 - 25.77*NF + 0.26*NF*NF)*(asMH/M_PI) ) );
    
    Deltamt = (asMH/M_PI) * (asMH/M_PI) * ( 1.57 + (4.0/3.0)*log(mHl/mtpole) + (4.0/9.0) * log(mf/mHl) * log(mf/mHl) );
    
    gamma = Nc * (4.0*GF/sqrt(2.0)) * (mf*mf/16.0/M_PI) * mHl * beta*beta*beta * (DeltaQCD + Deltamt);
    
    return gamma;
}

const double StandardModel::GammaHtobb() const
{   
    double mf0=quarks[BOTTOM].getMass(), mf;
    double beta;
    double Nc=3.0;
    double gamma;
    double asMH,DeltaQCD,Deltamt,NF;
    
    //      alfa_s(MH)
    asMH = Als(mHl, FULLNLO);
    
    mf = Mrun(mHl, mf0, mf0, BOTTOM, FULLNLO);
    
    beta=1.0-4.0*mf*mf/mHl/mHl;
    
    NF=5;
    
    DeltaQCD = 1 + (asMH/M_PI) * ( 17.0/3.0 + (asMH/M_PI) * ( (35.94 - 1.36*NF) + (164.14 - 25.77*NF + 0.26*NF*NF)*(asMH/M_PI) ) );
    
    Deltamt = (asMH/M_PI) * (asMH/M_PI) * ( 1.57 + (4.0/3.0)*log(mHl/mtpole) + (4.0/9.0) * log(mf/mHl) * log(mf/mHl) );
    
    gamma = Nc * (4.0*GF/sqrt(2.0)) * (mf*mf/16.0/M_PI) * mHl * beta*beta*beta * (DeltaQCD + Deltamt);
    
    return gamma;
}

const double StandardModel::GammaHTot() const
{   
    double gamma;
    
    gamma = GammaHtobb() + GammaHtocc() + GammaHtoss() + 
            GammaHtotautau() + GammaHtomumu() + 
            GammaHtoZZstar() + GammaHtoWWstar() +
            GammaHtogg() + GammaHtogaga() + GammaHtoZga();
    
    return gamma;
}

////////////////////////////////////////////////////////////////////////
// Higgs branching ratios
////////////////////////////////////////////////////////////////////////

const double StandardModel::BrHtogg() const
{       
    return GammaHtogg()/GammaHTot();
}

const double StandardModel::BrHtoZZstar() const
{       
    return GammaHtoZZstar()/GammaHTot();
}

const double StandardModel::BrHtoWWstar() const
{       
    return GammaHtoWWstar()/GammaHTot();
}

const double StandardModel::BrHtoZga() const
{       
    return GammaHtoZga()/GammaHTot();
}

const double StandardModel::BrHtogaga() const
{       
    return GammaHtogaga()/GammaHTot();
}    

const double StandardModel::BrHtomumu() const
{       
    return GammaHtomumu()/GammaHTot();
}

const double StandardModel::BrHtotautau() const
{       
    return GammaHtotautau()/GammaHTot();
}

const double StandardModel::BrHtocc() const
{       
    return GammaHtocc()/GammaHTot();
}   

const double StandardModel::BrHtoss() const
{       
    return GammaHtoss()/GammaHTot();
}

const double StandardModel::BrHtobb() const
{       
    return GammaHtobb()/GammaHTot();
}

////////////////////////////////////////////////////////////////////////
//Generic e+e- -> ff Inclusive Observables

//  For f!=e
//  (f=e also included to define a t-subtracted observable, like in LEP)
    
//  Helicity amplitudes squared
const double StandardModel::MLR2eeff(const Particle f, const double s) const {
    
    // Definitions      
    double Qf, geLSM, gfRSM, is2c2, GZ, Mz2s;
    
    double MLR2SM;

    // -------------------------------------------

    geLSM = (leptons[ELECTRON].getIsospin()) - (leptons[ELECTRON].getCharge()) * s02();

    is2c2 = 1. / s02() / c02();

    GZ = Gamma_Z();

    Mz2s = Mz * Mz - s;

    if (f.is("ELECTRON")) {
        Qf = leptons[ELECTRON].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("MU")) {
        Qf = leptons[MU].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("TAU")) {
        Qf = leptons[TAU].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("UP")) {
        Qf = quarks[UP].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("CHARM")) {
        Qf = quarks[CHARM].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("DOWN")) {
        Qf = quarks[DOWN].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("STRANGE")) {
        Qf = quarks[STRANGE].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("BOTTOM")) {
        Qf = quarks[BOTTOM].getCharge();
        gfRSM = - Qf * s02();
    } else
        throw std::runtime_error("StandardModel::MLR2eeff: wrong argument");
    
    // LR, RL, LL and RR SM squared amplitudes
    MLR2SM = Qf * Qf
            + (is2c2 * is2c2 * (geLSM * geLSM * gfRSM * gfRSM) * s * s
            + 2.0 * Qf * is2c2 * (geLSM * gfRSM) * Mz2s * s) / (Mz2s * Mz2s + Mz * Mz * GZ * GZ);

    return MLR2SM;
}
const double StandardModel::MRL2eeff(const Particle f, const double s) const{
    
    // Definitions      
    double Qf, geRSM, gfLSM, is2c2, GZ, Mz2s;
    
    double MRL2SM;

    // -------------------------------------------

    geRSM = - (leptons[ELECTRON].getCharge()) * s02();

    is2c2 = 1. / s02() / c02();

    GZ = Gamma_Z();

    Mz2s = Mz * Mz - s;

    if (f.is("ELECTRON")) {
        Qf = leptons[ELECTRON].getCharge();
        gfLSM = (leptons[ELECTRON].getIsospin()) - Qf * s02();
    } else if (f.is("MU")) {
        Qf = leptons[MU].getCharge();
        gfLSM = (leptons[MU].getIsospin()) - Qf * s02();
    } else if (f.is("TAU")) {
        Qf = leptons[TAU].getCharge();
        gfLSM = (leptons[TAU].getIsospin()) - Qf * s02();
    } else if (f.is("UP")) {
        Qf = quarks[UP].getCharge();
        gfLSM = (quarks[UP].getIsospin()) - Qf * s02();
    } else if (f.is("CHARM")) {
        Qf = quarks[CHARM].getCharge();
        gfLSM = (quarks[CHARM].getIsospin()) - Qf * s02();
    } else if (f.is("DOWN")) {
        Qf = quarks[DOWN].getCharge();
        gfLSM = (quarks[DOWN].getIsospin()) - Qf * s02();
    } else if (f.is("STRANGE")) {
        Qf = quarks[STRANGE].getCharge();
        gfLSM = (quarks[STRANGE].getIsospin()) - Qf * s02();
    } else if (f.is("BOTTOM")) {
        Qf = quarks[BOTTOM].getCharge();
        gfLSM = (quarks[BOTTOM].getIsospin()) - Qf * s02();
    } else
        throw std::runtime_error("StandardModel::MRL2eeff: wrong argument");
    
    // RL SM squared amplitude    
    MRL2SM = Qf * Qf
            + (is2c2 * is2c2 * (geRSM * geRSM * gfLSM * gfLSM) * s * s
            + 2.0 * Qf * is2c2 * (geRSM * gfLSM) * Mz2s * s) / (Mz2s * Mz2s + Mz * Mz * GZ * GZ); 
 
    return MRL2SM;
}

const double StandardModel::MLL2eeff(const Particle f, const double s, const double t) const{
    
    // Definitions      
    double Qf, geLSM, gfLSM, is2c2, GZ, Mz2s;
    
    double MLL2SM;

    // -------------------------------------------

    geLSM = (leptons[ELECTRON].getIsospin()) - (leptons[ELECTRON].getCharge()) * s02();

    is2c2 = 1. / s02() / c02();

    GZ = Gamma_Z();

    Mz2s = Mz * Mz - s;

    if (f.is("ELECTRON")) {
        Qf = leptons[ELECTRON].getCharge();
        gfLSM = (leptons[ELECTRON].getIsospin()) - Qf * s02();
    } else if (f.is("MU")) {
        Qf = leptons[MU].getCharge();
        gfLSM = (leptons[MU].getIsospin()) - Qf * s02();
    } else if (f.is("TAU")) {
        Qf = leptons[TAU].getCharge();
        gfLSM = (leptons[TAU].getIsospin()) - Qf * s02();
    } else if (f.is("UP")) {
        Qf = quarks[UP].getCharge();
        gfLSM = (quarks[UP].getIsospin()) - Qf * s02();
    } else if (f.is("CHARM")) {
        Qf = quarks[CHARM].getCharge();
        gfLSM = (quarks[CHARM].getIsospin()) - Qf * s02();
    } else if (f.is("DOWN")) {
        Qf = quarks[DOWN].getCharge();
        gfLSM = (quarks[DOWN].getIsospin()) - Qf * s02();
    } else if (f.is("STRANGE")) {
        Qf = quarks[STRANGE].getCharge();
        gfLSM = (quarks[STRANGE].getIsospin()) - Qf * s02();
    } else if (f.is("BOTTOM")) {
        Qf = quarks[BOTTOM].getCharge();
        gfLSM = (quarks[BOTTOM].getIsospin()) - Qf * s02();
    } else
        throw std::runtime_error("StandardModel::MLL2eeff: wrong argument");
    
    // LL SM squared amplitude    
    MLL2SM = Qf * Qf
            + (is2c2 * is2c2 * (geLSM * geLSM * gfLSM * gfLSM) * s * s
            + 2.0 * Qf * is2c2 * (geLSM * gfLSM) * Mz2s * s) / (Mz2s * Mz2s + Mz * Mz * GZ * GZ); 
    
    return MLL2SM;
        
}       
const double StandardModel::MRR2eeff(const Particle f, const double s, const double t) const{
    
    // Definitions      
    double Qf, geRSM, gfRSM, is2c2, GZ, Mz2s;
    
    double MRR2SM;

    // -------------------------------------------

    geRSM = - (leptons[ELECTRON].getCharge()) * s02();

    is2c2 = 1. / s02() / c02();

    GZ = Gamma_Z();

    Mz2s = Mz * Mz - s;

    if (f.is("ELECTRON")) {
        Qf = leptons[ELECTRON].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("MU")) {
        Qf = leptons[MU].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("TAU")) {
        Qf = leptons[TAU].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("UP")) {
        Qf = quarks[UP].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("CHARM")) {
        Qf = quarks[CHARM].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("DOWN")) {
        Qf = quarks[DOWN].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("STRANGE")) {
        Qf = quarks[STRANGE].getCharge();
        gfRSM = - Qf * s02();
    } else if (f.is("BOTTOM")) {
        Qf = quarks[BOTTOM].getCharge();
        gfRSM = - Qf * s02();
    } else
        throw std::runtime_error("StandardModel::MRR2eeff: wrong argument");
    
    // RR SM squared amplitude    
    MRR2SM = Qf * Qf
            + (is2c2 * is2c2 * (geRSM * geRSM * gfRSM * gfRSM) * s * s
            + 2.0 * Qf * is2c2 * (geRSM * gfRSM) * Mz2s * s) / (Mz2s * Mz2s + Mz * Mz * GZ * GZ); 
    
    return MRR2SM;
}

//  Some simple functions for cos \theta integrals 

const double StandardModel::tovers2(const double cosmin, const double cosmax) const {
    return 0.25 * (cosmax * (1.0 - cosmax * (1.0 - cosmax / 3.0)) - cosmin * (1.0 - cosmin * (1.0 - cosmin / 3.0)));
}

const double StandardModel::uovers2(const double cosmin, const double cosmax) const {
    return 0.25 * (cosmax * (1.0 + cosmax * (1.0 + cosmax / 3.0)) - cosmin * (1.0 + cosmin * (1.0 + cosmin / 3.0)));
}

//  Expressions for f=e   

//  Integrals of the SM squared amplitudes x (t/s)^2, (s/t)^2, (u/s)^2 in [t0, t1]    
const double StandardModel::intMLR2eeeets2(const double s, const double t0, const double t1) const {
    
    double intM2;
    double sw2cw2;
    double gLeSM,gReSM;
    double GammaZSM;
    double Mz2, s2;
    double propZSM2,propZSMRe,MeeLR2SM;
    
    sw2cw2 = s02() * c02();
    gLeSM = (leptons[ELECTRON].getIsospin()) - (leptons[ELECTRON].getCharge()) * s02();
    gReSM = - (leptons[ELECTRON].getCharge()) * s02();
    GammaZSM = Gamma_Z();
    Mz2 = Mz * Mz;
    s2 = s * s;
    
    propZSM2 = s2/((s - Mz2)*(s - Mz2) + Mz2*GammaZSM*GammaZSM);
    propZSMRe = (s*(s - Mz2))/((s - Mz2)*(s - Mz2) + Mz2*GammaZSM*GammaZSM);
    
    MeeLR2SM = 1.0 + (gLeSM*gLeSM*gReSM*gReSM/(sw2cw2*sw2cw2))*propZSM2 + 2.0*(gLeSM*gReSM/sw2cw2)*propZSMRe;

    intM2 = MeeLR2SM*(t1*t1*t1 - t0*t0*t0)/(3.0*s*s);

    return intM2;
}

const double StandardModel::intMLRtilde2eeeest2(const double s, const double t0, const double t1) const {
    
    double intM2;
    double sw2cw2; 
    double gLeSM,gReSM;
    double Mz2;
    
    sw2cw2 = s02() * c02();
    gLeSM = (leptons[ELECTRON].getIsospin()) - (leptons[ELECTRON].getCharge()) * s02();
    gReSM = - (leptons[ELECTRON].getCharge()) * s02();
    Mz2 = Mz * Mz;
    
    intM2 = s*s*(((gLeSM*gLeSM*gReSM*gReSM)/sw2cw2/sw2cw2)*(1.0/(Mz2 - t1) - 1.0/(Mz2 - t0)) - 1.0/t1 + 1.0/t0 + 
            (2.0*gLeSM*gReSM*(-log(t1/t0) + log((-Mz2 + t1)/(-Mz2 + t0))))/(Mz2*sw2cw2));

    return intM2;
}

const double StandardModel::intMLL2eeeeus2(const double s, const double t0, const double t1) const {
    
    double intM2;
    double sw2cw2; 
    double gLeSM;
    double GammaZSM;
    double Mz2, Mz4, s2;

    sw2cw2 = s02() * c02();
    gLeSM = (leptons[ELECTRON].getIsospin()) - (leptons[ELECTRON].getCharge()) * s02();
    GammaZSM = Gamma_Z();
    Mz2 = Mz * Mz;
    Mz4 = Mz2 * Mz2;
    s2 = s * s;
    
    intM2 = (gLeSM*gLeSM*gLeSM*gLeSM*s2 + 2.0*gLeSM*gLeSM*s*(-Mz2 + s)*sw2cw2 + sw2cw2*sw2cw2*(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM)))/(3.0*s2*sw2cw2*sw2cw2*(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM)))*(pow(s + t1,3.0) - pow(s + t0,3.0)) +
            ((2.0*(1.0 + (gLeSM*gLeSM*s*(-Mz2 + s))/(sw2cw2*(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM)))) )/s)*(2.0*s *(t1 - t0) + (t1*t1 - t0*t0)/2.0 + s2*log(t1/t0)) +
            (2.0*gLeSM*gLeSM* (-sw2cw2 + (gLeSM*gLeSM*(Mz2 - s)*s)/(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM))))/(s*sw2cw2*sw2cw2)* (-(1.0/2.0)*t1*(2.0*Mz2 + 4.0*s + t1) + (1.0/2.0)*t0*(2.0*Mz2 + 4.0*s + t0) - (Mz2 + s)*(Mz2 + s)*log((-Mz2 + t1)/(-Mz2 + t0)) ) +
            (2.0*(gLeSM*gLeSM) )/(Mz2*sw2cw2)*(Mz2 *(t1 - t0) - s2*log(t1/t0) + (Mz2 + s)*(Mz2 + s)*log((-Mz2 + t1)/(-Mz2 + t0))) +
            (-(s2/t1) + s2/t0 + t1 - t0 + 2.0*s*log(t1/t0)) +
            (gLeSM*gLeSM*gLeSM*gLeSM /sw2cw2/sw2cw2)*((Mz2 + s)*(Mz2 + s)*(1.0/(Mz2 - t1) - 1.0/(Mz2 - t0)) + t1 - t0 + 2.0*(Mz2 + s)*log((-Mz2 + t1)/(-Mz2 + t0)));            

    return intM2;
}

const double StandardModel::intMRR2eeeeus2(const double s, const double t0, const double t1) const {
    
    double intM2;
    double sw2cw2; 
    double gReSM;
    double GammaZSM;
    double Mz2, Mz4, s2;

    sw2cw2 = s02() * c02();
    gReSM = - (leptons[ELECTRON].getCharge()) * s02();
    GammaZSM = Gamma_Z();
    Mz2 = Mz * Mz;
    Mz4 = Mz2 * Mz2;
    s2 = s * s;
    
    intM2 = (gReSM*gReSM*gReSM*gReSM*s2 + 2.0*gReSM*gReSM*s*(-Mz2 + s)*sw2cw2 + sw2cw2*sw2cw2*(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM)))/(3.0*s2*sw2cw2*sw2cw2*(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM)))*(pow(s + t1,3.0) - pow(s + t0,3.0)) +
            ((2.0*(1.0 + (gReSM*gReSM*s*(-Mz2 + s))/(sw2cw2*(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM)))) )/s)*(2.0*s *(t1 - t0) + (t1*t1 - t0*t0)/2.0 + s2*log(t1/t0)) +
            (2.0*gReSM*gReSM* (-sw2cw2 + (gReSM*gReSM*(Mz2 - s)*s)/(Mz4 + s2 + Mz2*(-2.0*s + GammaZSM*GammaZSM))))/(s*sw2cw2*sw2cw2)* (-(1.0/2.0)*t1*(2.0*Mz2 + 4.0*s + t1) + (1.0/2.0)*t0*(2.0*Mz2 + 4.0*s + t0) - (Mz2 + s)*(Mz2 + s)*log((-Mz2 + t1)/(-Mz2 + t0)) ) +
            (2.0*(gReSM*gReSM) )/(Mz2*sw2cw2)*(Mz2 *(t1 - t0) - s2*log(t1/t0) + (Mz2 + s)*(Mz2 + s)*log((-Mz2 + t1)/(-Mz2 + t0))) +
            (-(s2/t1) + s2/t0 + t1 - t0 + 2.0*s*log(t1/t0)) +
            (gReSM*gReSM*gReSM*gReSM /sw2cw2/sw2cw2)*((Mz2 + s)*(Mz2 + s)*(1.0/(Mz2 - t1) - 1.0/(Mz2 - t0)) + t1 - t0 + 2.0*(Mz2 + s)*log((-Mz2 + t1)/(-Mz2 + t0)));            

    return intM2;
}

//  Cross sections

const double StandardModel::eeffsigmaEbin(const double pol_e, const double pol_p, const double s, const double cosmin, const double cosmax) const {
    
    double sumM2, sigma;
    double topb = 0.3894e+9; 
    double t0, t1, lambdaK;
    
    double pLH, pRH; //Polarization factors, minus the 1/4 average
    
    pLH = (1.0 - pol_e) * (1.0 + pol_p);
    pRH = (1.0 + pol_e) * (1.0 - pol_p);
    
    // t values for cosmin and cosmax
    t0 = 0.5 * s * ( -1.0 + cosmin );
    t1 = 0.5 * s * ( -1.0 + cosmax );
    
    // Kähllén function of (s,0,0)
    lambdaK = s*s;
    
    // Sum of the integrals of the amplitudes squared x (t/s)^2, (s/t)^2, (u/s)^2 
    sumM2 = (pLH + pRH) * ( intMLR2eeeets2(s, t0, t1) + intMLRtilde2eeeest2(s, t0, t1) ) + 
            pLH * intMLL2eeeeus2(s, t0, t1) + pRH * intMRR2eeeeus2(s, t0, t1);   
    
    // Build the cross section
    sigma = M_PI * (alphaMz())*(alphaMz()) * sumM2 / s / sqrt(lambdaK);
    
    return topb * sigma;
    
}
    
const double StandardModel::eeffsigma(const Particle f, const double pol_e, const double pol_p, const double s, const double cosmin, const double cosmax) const {
    //  Only valid for f=/=e (MLL2, MRR2 do not depend on t for f=/=e. Simply enter t=1 as argument)
    //  For f=e this corresponds to t-subtracted definition from LEP
    double sumM2, sigma;
    double tdumm = 1.;
    double topb = 0.3894e+9;
    
    //double cosmin = -1.0;
    //double cosmax = 1.0;

    double Nf;
    
    double pLH, pRH; //Polarization factors, minus the 1/4 average
    
    pLH = (1.0 - pol_e) * (1.0 + pol_p);
    pRH = (1.0 + pol_e) * (1.0 - pol_p);

    if (f.is("LEPTON")) {
        Nf = 1.0;
    } else {
        Nf = 3.0;
    }

    sumM2 = (pLH * MLR2eeff(f, s) + pRH * MRL2eeff(f, s)) * tovers2(cosmin, cosmax)
            + (pLH * MLL2eeff(f, s, tdumm) + pRH * MRR2eeff(f, s, tdumm)) * uovers2(cosmin, cosmax);

    sigma = Nf * 0.5 * M_PI * (alphaMz())*(alphaMz()) * sumM2 / s;

    return topb * sigma;
}

/* BEGIN: REMOVE FROM THE PACKAGE */
////////////////////////////////////////////////////////////////////////////////////
//LEP2 Observables

const double StandardModel::LEP2sigmaE(const double s) const
{
    return 0.;
}

const double StandardModel::LEP2sigmaMu(const double s) const
{
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-8;
    double abserr = 1.e-20;
    
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2sigmaMuApprox(s));

    } else {
    
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
}


const double StandardModel::LEP2sigmaTau(const double s) const
{
    
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;

    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2sigmaTauApprox(s));

    } else {   

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
}


const double StandardModel::LEP2sigmaCharm(const double s) const
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


const double StandardModel::LEP2sigmaBottom(const double s) const
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


const double StandardModel::LEP2sigmaHadron(const double s) const
{
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-8;
    double abserr = 1.e-20;
    
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2sigmaHadronApprox(s));

    } else {  
    
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
}


const double StandardModel::LEP2AFBbottom(const double s) const
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


const double StandardModel::LEP2AFBcharm(const double s) const
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

const double StandardModel::LEP2AFBe(const double s) const
{
    return 0.;
}

const double StandardModel::LEP2AFBmu(const double s) const
{

    bSigmaForAFB = true;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;

    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2AFBmuApprox(s));

    } else {  
         
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
}


const double StandardModel::LEP2AFBtau(const double s) const
{

    bSigmaForAFB = true;
    gsl_error_handler_t * old_handler = gsl_set_error_handler_off();
    double relerr = 1.e-7;
    double abserr = 1.e-17;
    
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2AFBtauApprox(s));

    } else {  

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
}


const double StandardModel::LEP2Rbottom(const double s) const
{

    double sigma_b = LEP2sigmaBottom(s);
    double sigma_had = LEP2sigmaHadron(s);
    SMresult_cache =  sigma_b / sigma_had;
    double R_bottom = SMresult_cache;
    
    return R_bottom;    
}


const double StandardModel::LEP2Rcharm(const double s) const
{

    double sigma_c = LEP2sigmaCharm(s);
    double sigma_had = LEP2sigmaHadron(s);
    SMresult_cache =  sigma_c / sigma_had;
    double R_charm = SMresult_cache;
    
    return R_charm;    
}


const double StandardModel::sigma_NoISR_l(const QCD::lepton l_flavor, const double s) const
{
    double ml = getLeptons(l_flavor).getMass();
    double l_charge = getLeptons(l_flavor).getCharge();
    double sigma = myTwoFermionsLEP2->sigma_l(l_flavor, ml, s, Mw(), Gamma_Z(), flagLEP2[Weak]);

    if (!bSigmaForAFB && flagLEP2[QEDFSR])
        sigma *= myTwoFermionsLEP2->QED_FSR_forSigma(s, l_charge);

    return sigma;
}   
    
const double StandardModel::sigma_NoISR_q(const QCD::quark q_flavor, const double s) const
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
    
const double StandardModel::AFB_NoISR_l(const QCD::lepton l_flavor, const double s) const
{
    double ml = getLeptons(l_flavor).getMass();
    double AFB = myTwoFermionsLEP2->AFB_l(l_flavor, ml, s, Mw(), Gamma_Z(), flagLEP2[Weak]);

    return AFB;
}
    
const double StandardModel::AFB_NoISR_q(const QCD::quark q_flavor, const  double s) const
{
    double mq = m_q(q_flavor, sqrt(s));
    double AFB = myTwoFermionsLEP2->AFB_q(q_flavor, mq, s, Mw(), Gamma_Z(), flagLEP2[Weak]);
        
    if (flagLEP2[QCDFSR])
        AFB *= myTwoFermionsLEP2->QCD_FSR_forAFB(q_flavor, mq, s);
        
    return AFB;
}
    
const double StandardModel::Integrand_sigmaWithISR_l(double x, const QCD::lepton l_flavor, const  double s) const
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
    
const double StandardModel::getIntegrand_sigmaWithISR_mu130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_mu207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(MU), s));
}
    

const double StandardModel::getIntegrand_sigmaWithISR_tau130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_tau136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_tau207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_l(x, QCD::lepton(TAU), s));
}
    
const double StandardModel::Integrand_sigmaWithISR_q(double x, const QCD::quark q_flavor, const  double s) const
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


const double StandardModel::getIntegrand_sigmaWithISR_up130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_up133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_up207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(UP), s));
}


//down

const double StandardModel::getIntegrand_sigmaWithISR_down130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_down133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}


const double StandardModel::getIntegrand_sigmaWithISR_down136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_down207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(DOWN), s));
}


//charm


const double StandardModel::getIntegrand_sigmaWithISR_charm130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_charm133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_charm136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_charm207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(CHARM), s));
}


//strange


const double StandardModel::getIntegrand_sigmaWithISR_strange130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_strange133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_strange207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(STRANGE), s));
}


//bottom


const double StandardModel::getIntegrand_sigmaWithISR_bottom130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_bottom133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}    

const double StandardModel::getIntegrand_sigmaWithISR_bottom136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_sigmaWithISR_bottom207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_sigmaWithISR_q(x, QCD::quark(BOTTOM), s));
}





const double StandardModel::Integrand_dsigmaBox_l(double cosTheta, const QCD::lepton l_flavor, const double s) const
{
    double ml = getLeptons(l_flavor).getMass();
    return ( myTwoFermionsLEP2->dsigma_l_box(l_flavor, ml, s, cosTheta, Mw(), Gamma_Z()) );
}       
    
const double StandardModel::getIntegrand_dsigmaBox_mu130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_mu207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(MU), s));
}





const double StandardModel::getIntegrand_dsigmaBox_tau130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_dsigmaBox_tau207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_l(x, QCD::lepton(TAU), s));
}






const double StandardModel::Integrand_dsigmaBox_q(double cosTheta, const QCD::quark q_flavor, const  double s) const
{
    double mq = m_q(q_flavor, sqrt(s)); 
    return ( myTwoFermionsLEP2->dsigma_q_box(q_flavor, mq, s, cosTheta, Mw(), Gamma_Z()) );
}       




//up


const double StandardModel::getIntegrand_dsigmaBox_up130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_up133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}

const double StandardModel::getIntegrand_dsigmaBox_up207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(UP), s));
}


//down

const double StandardModel::getIntegrand_dsigmaBox_down130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_down133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}

const double StandardModel::getIntegrand_dsigmaBox_down207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(DOWN), s));
}


//charm


const double StandardModel::getIntegrand_dsigmaBox_charm130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_charm133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_charm136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_charm207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(CHARM), s));
}


//strange


const double StandardModel::getIntegrand_dsigmaBox_strange130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_strange133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}



const double StandardModel::getIntegrand_dsigmaBox_strange172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}

const double StandardModel::getIntegrand_dsigmaBox_strange207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(STRANGE), s));
}


//bottom


const double StandardModel::getIntegrand_dsigmaBox_bottom130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_bottom133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}    

const double StandardModel::getIntegrand_dsigmaBox_bottom136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_dsigmaBox_bottom207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_dsigmaBox_q(x, QCD::quark(BOTTOM), s));
}












    
const double StandardModel::Integrand_AFBnumeratorWithISR_l(double x, const QCD::lepton l_flavor, const  double s) const
{
    double sprime = (1.0 - x)*s;
    double Ncf = 1.0;
    double ml = getLeptons(l_flavor).getMass();
    double G3prime = myTwoFermionsLEP2->G_3prime_l(l_flavor, ml, sprime, Mw(), Gamma_Z(),flagLEP2[Weak]);
    double H = myTwoFermionsLEP2->H_ISR_FB(x, s);

    return ( M_PI*ale*ale*Ncf*H*G3prime/sprime );
}
    
    
    const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_mu207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(MU), s));
}
    

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau130(double x) const
{
    double s = 130. * 130.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}    

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau136(double x) const
{
    double s = 136. * 136.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau161(double x) const
{
    double s = 161. * 161.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_tau207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_l(x, QCD::lepton(TAU), s));
}

    
    
const double StandardModel::Integrand_AFBnumeratorWithISR_q(double x, const QCD::quark q_flavor, const  double s) const
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


const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}    

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_charm207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(CHARM), s));
}    
    
    


const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom133(double x) const
{
    double s = 133. * 133.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}    

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom167(double x) const
{
    double s = 167. * 167.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom172(double x) const
{
    double s = 172. * 172.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom183(double x) const
{
    double s = 183. * 183.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom189(double x) const
{
    double s = 189. * 189.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom192(double x) const
{
    double s = 192. * 192.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom196(double x) const
{
    double s = 196. * 196.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom200(double x) const
{
    double s = 200. * 200.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom202(double x) const
{
    double s = 202. * 202.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom205(double x) const
{
    double s = 205. * 205.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}

const double StandardModel::getIntegrand_AFBnumeratorWithISR_bottom207(double x) const
{
    double s = 207. * 207.;
    return (Integrand_AFBnumeratorWithISR_q(x, QCD::quark(BOTTOM), s));
}



//  LEP2 differential observables

const double StandardModel::LEP2dsigmadcosE(const double s, const double cos) const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2dsigmadcosEApprox(s, cos));

    } else {
        throw std::runtime_error("ERROR: StandardModel::LEP2dsigmadcosE only implemented via semi-analytical approx");
    }
}

const double StandardModel::LEP2dsigmadcosMu(double s, double cos) const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2dsigmadcosMuApprox(s, cos));

    } else {
        throw std::runtime_error("ERROR: StandardModel::LEP2dsigmadcosMu only implemented via semi-analytical approx");
    }
}

const double StandardModel::LEP2dsigmadcosTau(double s, double cos) const
{
    // Use same flag as other Z pole observables for the moment to decide whether to use approx formulae
    if (!IsFlagNoApproximateGammaZ()){
            
    /* SM contribution with the approximate formula */
        return (myApproximateFormulae->LEP2dsigmadcosTauApprox(s, cos));

    } else {
        throw std::runtime_error("ERROR: StandardModel::LEP2dsigmadcosTau only implemented via semi-analytical approx");
    }
}


//  LEP2 differential observables: Defined in bins. SM prediction given already by the above

const double StandardModel::LEP2dsigmadcosBinE(const double s, const double cos, const double cosmin, const double cosmax) const
{
    return LEP2dsigmadcosE(s, cos);
}

const double StandardModel::LEP2dsigmadcosBinMu(double s, double cos, const double cosmin, const double cosmax) const
{
    return LEP2dsigmadcosMu(s, cos);
}

const double StandardModel::LEP2dsigmadcosBinTau(double s, double cos, const double cosmin, const double cosmax) const
{
    return LEP2dsigmadcosTau(s, cos);
}


/* END: REMOVE FROM THE PACKAGE */
