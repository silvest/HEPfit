/* 
 * Copyright (C) 2017 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "THDMW.h"
#include "THDMWcache.h"

const std::string THDMW::THDMWvars[NTHDMWvars] = {"THDMW_logtb","THDMW_bma",
                                               "THDMW_lambda1","THDMW_lambda2","THDMW_lambda3","THDMW_lambda4","THDMW_lambda5",
                                               "THDMW_mS2","THDMW_mu1","THDMW_mu2","THDMW_mu3","THDMW_mu4","THDMW_mu5","THDMW_mu6",
                                               "THDMW_nu1","THDMW_nu2","THDMW_nu3","THDMW_nu4","THDMW_nu5",
                                               "THDMW_omega1","THDMW_omega2","THDMW_omega3","THDMW_omega4","THDMW_omega5",
                                               "THDMW_kappa1","THDMW_kappa2","THDMW_kappa3",
                                               "THDMW_etaU","THDMW_etaD",
                                               "THDMW_rho_b", "THDMW_S_b",
                                               "Q_THDMW","RpepsTHDMW","NLOuniscaleTHDMW"};

THDMW::THDMW() : NPbase(), THDMWM(*this) {

    SMM.setObj((StandardModelMatching&) THDMWM.getObj());
    ModelParamMap.insert(std::make_pair("THDMW_logtb", std::cref(THDMW_logtb)));
    ModelParamMap.insert(std::make_pair("THDMW_bma", std::cref(THDMW_bma)));
    ModelParamMap.insert(std::make_pair("THDMW_lambda1", std::cref(THDMW_lambda1)));
    ModelParamMap.insert(std::make_pair("THDMW_lambda2", std::cref(THDMW_lambda2)));
    ModelParamMap.insert(std::make_pair("THDMW_lambda3", std::cref(THDMW_lambda3)));
    ModelParamMap.insert(std::make_pair("THDMW_lambda4", std::cref(THDMW_lambda4)));
    ModelParamMap.insert(std::make_pair("THDMW_lambda5", std::cref(THDMW_lambda5)));
    ModelParamMap.insert(std::make_pair("THDMW_mS2", std::cref(THDMW_mS2)));
    ModelParamMap.insert(std::make_pair("THDMW_mu1", std::cref(THDMW_mu1)));
    ModelParamMap.insert(std::make_pair("THDMW_mu2", std::cref(THDMW_mu2)));
    ModelParamMap.insert(std::make_pair("THDMW_mu3", std::cref(THDMW_mu3)));
    ModelParamMap.insert(std::make_pair("THDMW_mu4", std::cref(THDMW_mu4)));
    ModelParamMap.insert(std::make_pair("THDMW_mu5", std::cref(THDMW_mu5)));
    ModelParamMap.insert(std::make_pair("THDMW_mu6", std::cref(THDMW_mu6)));
    ModelParamMap.insert(std::make_pair("THDMW_nu1", std::cref(THDMW_nu1)));
    ModelParamMap.insert(std::make_pair("THDMW_nu2", std::cref(THDMW_nu2)));
    ModelParamMap.insert(std::make_pair("THDMW_nu3", std::cref(THDMW_nu3)));
    ModelParamMap.insert(std::make_pair("THDMW_nu4", std::cref(THDMW_nu4)));
    ModelParamMap.insert(std::make_pair("THDMW_nu5", std::cref(THDMW_nu5)));
    ModelParamMap.insert(std::make_pair("THDMW_omega1", std::cref(THDMW_omega1)));
    ModelParamMap.insert(std::make_pair("THDMW_omega2", std::cref(THDMW_omega2)));
    ModelParamMap.insert(std::make_pair("THDMW_omega3", std::cref(THDMW_omega3)));
    ModelParamMap.insert(std::make_pair("THDMW_omega4", std::cref(THDMW_omega4)));
    ModelParamMap.insert(std::make_pair("THDMW_omega5", std::cref(THDMW_omega5)));
    ModelParamMap.insert(std::make_pair("THDMW_kappa1", std::cref(THDMW_kappa1)));
    ModelParamMap.insert(std::make_pair("THDMW_kappa2", std::cref(THDMW_kappa2)));
    ModelParamMap.insert(std::make_pair("THDMW_kappa3", std::cref(THDMW_kappa3)));
    ModelParamMap.insert(std::make_pair("THDMW_etaU", std::cref(THDMW_etaU)));
    ModelParamMap.insert(std::make_pair("THDMW_etaD", std::cref(THDMW_etaD)));
    ModelParamMap.insert(std::make_pair("THDMW_rho_b", std::cref(THDMW_rho_b)));
    ModelParamMap.insert(std::make_pair("THDMW_S_b", std::cref(THDMW_S_b)));
    ModelParamMap.insert(std::make_pair("THDMW_kappa3", std::cref(THDMW_kappa3)));
    ModelParamMap.insert(std::make_pair("Q_THDMW", std::cref(Q_THDMW)));
    ModelParamMap.insert(std::make_pair("RpepsTHDMW", std::cref(RpepsTHDMW)));
    ModelParamMap.insert(std::make_pair("NLOuniscaleTHDMW", std::cref(NLOuniscaleTHDMW)));
}

THDMW::~THDMW(){
    if (IsModelInitialized()) {
            if (myTHDMWcache != NULL) delete(myTHDMWcache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool THDMW::InitializeModel()
{
    myTHDMWcache = new THDMWcache(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelTHDMW();
    return(true);
}

bool THDMW::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool THDMW::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool THDMW::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool THDMW::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    myTHDMWcache->updateCache();

    return (true);
}

void THDMW::setParameter(const std::string name, const double& value){

    if(name.compare("THDMW_logtb") == 0) {
        THDMW_logtb = value;
        THDMW_tanb = pow(10.,THDMW_logtb);
        if(THDMW_tanb > 0.) {
            THDMW_sinb = THDMW_tanb / sqrt(1. + THDMW_tanb*THDMW_tanb);
            THDMW_cosb = 1. / sqrt(1. + THDMW_tanb*THDMW_tanb);
        }
        else {
            throw std::runtime_error("error in THDMW::SetParameter, THDMW_tanb < 0!"); 
          }
        }
    else if(name.compare("THDMW_bma") == 0) {
        THDMW_bma = value;
        THDMW_sin_ba = sin(THDMW_bma);
    }
    else if(name.compare("THDMW_lambda1") == 0)
        THDMW_lambda1 = value;
    else if(name.compare("THDMW_lambda2") == 0)
        THDMW_lambda2 = value;
    else if(name.compare("THDMW_lambda3") == 0)
        THDMW_lambda3 = value;
    else if(name.compare("THDMW_lambda4") == 0)
        THDMW_lambda4 = value;
    else if(name.compare("THDMW_lambda5") == 0 && flag_modelTHDMW!="custodial1")
        THDMW_lambda5 = value;
    else if(name.compare("THDMW_mS2") == 0)
        THDMW_mS2 = value;
    else if(name.compare("THDMW_mu1") == 0)
        THDMW_mu1 = value;
    else if(name.compare("THDMW_mu2") == 0)
        THDMW_mu2 = value;
    else if(name.compare("THDMW_mu3") == 0)
        THDMW_mu3 = value;
    else if(name.compare("THDMW_mu4") == 0)
        THDMW_mu4 = value;
    else if(name.compare("THDMW_mu5") == 0)
        THDMW_mu5 = value;
    else if(name.compare("THDMW_mu6") == 0)
        THDMW_mu6 = value;
    else if(name.compare("THDMW_nu1") == 0)
        THDMW_nu1 = value;
    else if(name.compare("THDMW_nu2") == 0)
        THDMW_nu2 = value;
    else if(name.compare("THDMW_nu3") == 0 && flag_modelTHDMW!="custodial1")
        THDMW_nu3 = value;
    else if(name.compare("THDMW_nu4") == 0)
        THDMW_nu4 = value;
    else if(name.compare("THDMW_nu5") == 0 && flag_modelTHDMW!="custodial1")
        THDMW_nu5 = value;
    else if(name.compare("THDMW_omega1") == 0)
        THDMW_omega1 = value;
    else if(name.compare("THDMW_omega2") == 0)
        THDMW_omega2 = value;
    else if(name.compare("THDMW_omega3") == 0 && flag_modelTHDMW!="custodial1")
        THDMW_omega3 = value;
    else if(name.compare("THDMW_omega4") == 0)
        THDMW_omega4 = value;
    else if(name.compare("THDMW_omega5") == 0 && flag_modelTHDMW!="custodial1")
        THDMW_omega5 = value;
    else if(name.compare("THDMW_kappa1") == 0)
        THDMW_kappa1 = value;
    else if(name.compare("THDMW_kappa2") == 0)
        THDMW_kappa2 = value;
    else if(name.compare("THDMW_kappa3") == 0 && flag_modelTHDMW!="custodial1")
        THDMW_kappa3 = value;
    else if(name.compare("THDMW_etaU") == 0)
        THDMW_etaU = value;
    else if(name.compare("THDMW_etaD") == 0)
        THDMW_etaD = value;
    else if(name.compare("THDMW_rho_b") == 0)
        THDMW_rho_b = value;
    else if(name.compare("THDMW_S_b") == 0)
        THDMW_S_b = value;
    else if(name.compare("Q_THDMW") == 0)
        Q_THDMW = value;
    else if(name.compare("RpepsTHDMW") == 0)
        RpepsTHDMW = value;
    else if(name.compare("NLOuniscaleTHDMW") == 0)
        NLOuniscaleTHDMW = value;
    else
        StandardModel::setParameter(name,value);
}

bool THDMW::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NTHDMWvars; i++) {
        if (DPars.find(THDMWvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory THDMW parameter " << THDMWvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(THDMWvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool THDMW::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if(name.compare("RGEorder") == 0) {
        if (checkRGEorder(value))
        {
            flag_RGEorder = value;
            res = true;
        }
        else
        {
            throw std::runtime_error("THDMW::setFlagStr(): Invalid flag "
                + name + "=" + value);
        }
    }
    else if(name.compare("modelTypeTHDMW") == 0) {
        if (checkmodelTypeTHDMW(value))
        {
            flag_modelTHDMW = value;
            res = true;
            if (flag_modelTHDMW == "custodial1") {
                //
            }
            else if (flag_modelTHDMW == "ManoharWise") {
                //Pure Manohar-Wise model
            }
            else if (flag_modelTHDMW == "custodialMW") {
                //Custodial Manohar-Wise model
            }
            else if (flag_modelTHDMW == "custodial2") {
                //
            }
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

const double THDMW::muggH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_gg;
}

const double THDMW::muVBF(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::mueeWBF(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::muWH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::muZH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::mueeZH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::muVH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::muVBFpVH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV;
}

const double THDMW::muttH(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_QuQu;
}

const double THDMW::computeGammaTotalRatio() const
{
    return getMyTHDMWCache()->sumModBRs;
}

const double THDMW::GammaTotal() const
{
    return getMyTHDMWCache()->Gamma_h;
}

const double THDMW::BrHggRatio() const
{
    return getMyTHDMWCache()->rh_gg / computeGammaTotalRatio();
}

const double THDMW::BrHWWRatio() const
{
    return getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::BrHZZRatio() const
{
    return getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::BrHZgaRatio() const
{
    return getMyTHDMWCache()->rh_Zga / computeGammaTotalRatio();
}

const double THDMW::BrHgagaRatio() const
{
    //std::cout<<"rh_gaga="<< getMyTHDMWCache()->rh_gaga /computeGammaTotalRatio() << std::endl;
    return getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

const double THDMW::BrHmumuRatio() const
{
    return getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMW::BrHtautauRatio() const
{
    return getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMW::BrHccRatio() const
{
    return getMyTHDMWCache()->rh_QuQu / computeGammaTotalRatio();
}

const double THDMW::BrHbbRatio() const
{
    return getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

const double THDMW::muggHgaga(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_gg * getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

const double THDMW::muVBFHgaga(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

const double THDMW::muVHgaga(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

const double THDMW::muttHgaga(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_QuQu * getMyTHDMWCache()->rh_gaga / computeGammaTotalRatio();
}

const double THDMW::muggHZZ(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_gg * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muVBFHZZ(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muVHZZ(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muttHZZ(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_QuQu * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muggHWW(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_gg * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muVBFHWW(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muVHWW(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muttHWW(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_QuQu * getMyTHDMWCache()->rh_VV / computeGammaTotalRatio();
}

const double THDMW::muggHtautau(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_gg * getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMW::muVBFHtautau(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMW::muVHtautau(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMW::muttHtautau(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_QuQu * getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
}

const double THDMW::muggHbb(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_gg * getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

const double THDMW::muVBFHbb(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

const double THDMW::muVHbb(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_VV * getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

const double THDMW::muttHbb(const double sqrt_s) const
{
    return getMyTHDMWCache()->rh_QuQu * getMyTHDMWCache()->rh_QdQd / computeGammaTotalRatio();
}

const double THDMW::muppHmumu(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * getMyTHDMWCache()->rh_gg + 0.122 * getMyTHDMWCache()->rh_VV + 0.006 * getMyTHDMWCache()->rh_QuQu) * getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
    }
    if(sqrt_s==13)
    {
        return (0.871 * getMyTHDMWCache()->rh_gg + 0.119 * getMyTHDMWCache()->rh_VV + 0.010 * getMyTHDMWCache()->rh_QuQu) * getMyTHDMWCache()->rh_ll / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHmumu is only defined for 8 or 13 TeV.");
    }
}

const double THDMW::muppHZga(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * getMyTHDMWCache()->rh_gg + 0.122 * getMyTHDMWCache()->rh_VV + 0.006 * getMyTHDMWCache()->rh_QuQu) * getMyTHDMWCache()->rh_Zga / computeGammaTotalRatio();
    }
    if(sqrt_s==13)
    {
        return (0.871 * getMyTHDMWCache()->rh_gg + 0.119 * getMyTHDMWCache()->rh_VV + 0.010 * getMyTHDMWCache()->rh_QuQu) * getMyTHDMWCache()->rh_Zga / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHZga is only defined for 8 or 13 TeV.");
    }
}

const double THDMW::Mw() const{
    double MZ = StandardModel::Mz;
    return ( MZ / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0 * M_PI * StandardModel::ale / (sqrt(2.0) * StandardModel::GF * MZ* MZ))));
}
