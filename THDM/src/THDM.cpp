/* 
 * Copyright (C) 2012 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDM.h"
#include "THDMcache.h"

std::string THDM::THDMvars[NTHDMvars] = {"logtb","bma",/*"mHh1","mA1","mHp1",*/"mHh2","mA2","mHp2","m12_2",
                                               "BDtaunu_SM","BDtaunu_A","BDtaunu_B","BDstartaunu_SM","BDstartaunu_A","BDstartaunu_B",
                                               "bsgamma_theoryerror","Q_THDM","Rpeps","NLOuniscale"};

THDM::THDM() : StandardModel(), THDMM(*this) {
    SMM.setObj((StandardModelMatching&) THDMM.getObj());
    ModelParamMap.insert(std::make_pair("logtb", std::cref(logtb)));
    ModelParamMap.insert(std::make_pair("bma", std::cref(bma)));
    ModelParamMap.insert(std::make_pair("mHh1", std::cref(mHh1)));
    ModelParamMap.insert(std::make_pair("mA1", std::cref(mA1)));
    ModelParamMap.insert(std::make_pair("mHp1", std::cref(mHp1)));
    ModelParamMap.insert(std::make_pair("mHh2", std::cref(mHh2)));
    ModelParamMap.insert(std::make_pair("mA2", std::cref(mA2)));
    ModelParamMap.insert(std::make_pair("mHp2", std::cref(mHp2)));
    ModelParamMap.insert(std::make_pair("m12_2", std::cref(m12_2)));
    ModelParamMap.insert(std::make_pair("BDtaunu_SM", std::cref(BDtaunu_SM)));
    ModelParamMap.insert(std::make_pair("BDtaunu_A", std::cref(BDtaunu_A)));
    ModelParamMap.insert(std::make_pair("BDtaunu_B", std::cref(BDtaunu_B)));
    ModelParamMap.insert(std::make_pair("BDstartaunu_SM", std::cref(BDstartaunu_SM)));
    ModelParamMap.insert(std::make_pair("BDstartaunu_A", std::cref(BDstartaunu_A)));
    ModelParamMap.insert(std::make_pair("BDstartaunu_B", std::cref(BDstartaunu_B)));
    ModelParamMap.insert(std::make_pair("bsgamma_theoryerror", std::cref(bsgamma_theoryerror)));
    ModelParamMap.insert(std::make_pair("Q_THDM", std::cref(Q_THDM)));
    ModelParamMap.insert(std::make_pair("Rpeps", std::cref(Rpeps)));
    ModelParamMap.insert(std::make_pair("NLOuniscale", std::cref(NLOuniscale)));
    flag_use_sq_masses=true;
}

THDM::~THDM(){
    if (IsModelInitialized()) {
            if (myTHDMcache != NULL) delete(myTHDMcache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool THDM::InitializeModel()
{
    myTHDMcache = new THDMcache(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelTHDM();
    return(true);
}

bool THDM::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool THDM::PreUpdate()
{
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
}

bool THDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    mHl2=myTHDMcache->updateCache();

    return (true);
}

void THDM::setParameter(const std::string name, const double& value){

    if(name.compare("logtb") == 0) {
        logtb = value;
        tanb = pow(10.,logtb);
        if(tanb > 0.) {
            sinb = tanb / sqrt(1. + tanb*tanb);
            cosb = 1. / sqrt(1. + tanb*tanb);
        }
        else {
            throw std::runtime_error("error in THDM::SetParameter, tanb < 0!"); 
          }
        }
    else if(name.compare("bma") == 0) {
        bma = value;
        sin_ba = sin(bma);
    }
    else if(name.compare("mHh1") == 0 && !flag_use_sq_masses)
        mHh1 = value;
    else if(name.compare("mA1") == 0 && !flag_use_sq_masses)
        mA1 = value;
    else if(name.compare("mHp1") == 0 && !flag_use_sq_masses)
        mHp1 = value;
    else if(name.compare("mHh2") == 0 && flag_use_sq_masses)
        mHh2 = value;
    else if(name.compare("mA2") == 0 && flag_use_sq_masses)
        mA2 = value;
    else if(name.compare("mHp2") == 0 && flag_use_sq_masses)
        mHp2 = value;
    else if(name.compare("m12_2") == 0)
        m12_2 = value;
    else if(name.compare("BDtaunu_SM") == 0)
        BDtaunu_SM = value;
    else if(name.compare("BDtaunu_A") == 0)
        BDtaunu_A = value;
    else if(name.compare("BDtaunu_B") == 0)
        BDtaunu_B = value;
    else if(name.compare("BDstartaunu_SM") == 0)
        BDstartaunu_SM = value;
    else if(name.compare("BDstartaunu_A") == 0)
        BDstartaunu_A = value;
    else if(name.compare("BDstartaunu_B") == 0)
        BDstartaunu_B = value;
    else if(name.compare("bsgamma_theoryerror") == 0)
        bsgamma_theoryerror = value;
    else if(name.compare("Q_THDM") == 0)
        Q_THDM = value;
    else if(name.compare("Rpeps") == 0)
        Rpeps = value;
    else if(name.compare("NLOuniscale") == 0)
        NLOuniscale = value;
    else
        StandardModel::setParameter(name,value);
}

bool THDM::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NTHDMvars; i++) {
        if (DPars.find(THDMvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory THDM parameter " << THDMvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(THDMvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool THDM::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if(name.compare("modelTypeflag") == 0)
    {
        if (checkmodelType(value))
        {
            flag_model = value;
            res = true;
        }
        else
        {
            throw std::runtime_error("THDM::setFlagStr(): Invalid flag "
                + name + "=" + value);
        }
    }
    else if(name.compare("RGEorder") == 0) {
        if (checkRGEorder(value))
        {
            flag_RGEorder = value;
            res = true;
        }
        else
        {
            throw std::runtime_error("THDM::setFlagStr(): Invalid flag "
                + name + "=" + value);
        }
    }
    else
    {
        res = StandardModel::setFlagStr(name,value);
    }

    return(res);
}

bool THDM::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if(name.compare("use_sq_masses") == 0) {
        flag_use_sq_masses = value;
        res = true;
        if (!flag_use_sq_masses) {
            THDMvars[std::distance(THDMvars,std::find(THDMvars,THDMvars+NTHDMvars,"mHh2"))] = "mHh1";
            THDMvars[std::distance(THDMvars,std::find(THDMvars,THDMvars+NTHDMvars,"mA2"))] = "mA1";
            THDMvars[std::distance(THDMvars,std::find(THDMvars,THDMvars+NTHDMvars,"mHp2"))] = "mHp1";
        }
    }
    else if(name.compare("wavefunctionrenormalization") == 0) {
        flag_wfr = value;
        res = true;
    }
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}
