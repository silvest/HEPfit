/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "LeftRightSymmetricModel.h"
#include "LRSMquantities.h"

const std::string LeftRightSymmetricModel::LeftRightSymmetricModelvars[NLeftRightSymmetricModelvars] =
                        {"mH1p_2", "mH2p_2", "mdeltappR_2", "xi_LRSM", "mWR", 
                         "lambda1_LRSM", "lambda2_LRSM", "lambda3_LRSM", "lambda4_LRSM", 
                         "rho1_LRSM", "rho4_LRSM", "alpha1_LRSM", "alpha2_LRSM",
                         "alpha_LRSM", "delta2_LRSM",
                         "thetaR12_LRSM","thetaR13_LRSM","thetaR23_LRSM",
                         "phiR1_LRSM","phiR2_LRSM","phiR3_LRSM","phiR4_LRSM","phiR5_LRSM","phiR6_LRSM",
                         "Q_LRSM"};

LeftRightSymmetricModel::LeftRightSymmetricModel() : StandardModel(),
U(5, 5, 0.), LRSMM(*this)
{
    SMM.setObj((StandardModelMatching&) LRSMM.getObj());
    ModelParamMap.insert(std::make_pair("mH1p_2", std::cref(mH1p_2)));
    ModelParamMap.insert(std::make_pair("mH2p_2", std::cref(mH2p_2)));
    ModelParamMap.insert(std::make_pair("mdeltappR_2", std::cref(mdeltappR_2)));
    ModelParamMap.insert(std::make_pair("xi_LRSM", std::cref(xi_LRSM)));
    ModelParamMap.insert(std::make_pair("mWR", std::cref(mWR)));
    ModelParamMap.insert(std::make_pair("lambda1_LRSM", std::cref(lambda1_LRSM)));
    ModelParamMap.insert(std::make_pair("lambda2_LRSM", std::cref(lambda2_LRSM)));
    ModelParamMap.insert(std::make_pair("lambda3_LRSM", std::cref(lambda3_LRSM)));
    ModelParamMap.insert(std::make_pair("lambda4_LRSM", std::cref(lambda4_LRSM)));
    ModelParamMap.insert(std::make_pair("rho1_LRSM", std::cref(rho1_LRSM)));
    ModelParamMap.insert(std::make_pair("rho4_LRSM", std::cref(rho4_LRSM)));
    ModelParamMap.insert(std::make_pair("alpha1_LRSM", std::cref(alpha1_LRSM)));
    ModelParamMap.insert(std::make_pair("alpha2_LRSM", std::cref(alpha2_LRSM)));
    ModelParamMap.insert(std::make_pair("alpha_LRSM", std::cref(alpha_LRSM)));
    ModelParamMap.insert(std::make_pair("delta2_LRSM", std::cref(delta2_LRSM)));
    ModelParamMap.insert(std::make_pair("thetaR12_LRSM", std::cref(thetaR12_LRSM)));
    ModelParamMap.insert(std::make_pair("thetaR13_LRSM", std::cref(thetaR13_LRSM)));
    ModelParamMap.insert(std::make_pair("thetaR23_LRSM", std::cref(thetaR23_LRSM)));
    ModelParamMap.insert(std::make_pair("phiR1_LRSM", std::cref(phiR1_LRSM)));
    ModelParamMap.insert(std::make_pair("phiR2_LRSM", std::cref(phiR2_LRSM)));
    ModelParamMap.insert(std::make_pair("phiR3_LRSM", std::cref(phiR3_LRSM)));
    ModelParamMap.insert(std::make_pair("phiR4_LRSM", std::cref(phiR4_LRSM)));
    ModelParamMap.insert(std::make_pair("phiR5_LRSM", std::cref(phiR5_LRSM)));
    ModelParamMap.insert(std::make_pair("phiR6_LRSM", std::cref(phiR6_LRSM)));
    ModelParamMap.insert(std::make_pair("Q_LRSM", std::cref(Q_LRSM)));
}

LeftRightSymmetricModel::~LeftRightSymmetricModel(){
        if (myLRSMquantities != NULL) delete(myLRSMquantities);
//    if (IsModelInitialized()) {
//        }
}

/*----------------------------------------------------------------------------*/
// Initialization

bool LeftRightSymmetricModel::InitializeModel()
{
    myLRSMquantities = new LRSMquantities(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelLeftRightSymmetricModel();
    return(true);
}

bool LeftRightSymmetricModel::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool LeftRightSymmetricModel::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool LeftRightSymmetricModel::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);
    UpdateError = false;
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    if (UpdateError) return (false);
    if(!PostUpdate()) return (false);

    return (true);
}

bool LeftRightSymmetricModel::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);
    if(!myLRSMquantities->CalcNeutralMasses(U,mH0sq)) return (false);
    if(!myLRSMquantities->CalcNeutralMasses_app(mH0sq_app)) return (false);
    
    LRSMM.getObj().updateLeftRightSymmetricModelParameters();

    return (true);
}

void LeftRightSymmetricModel::setParameter(const std::string name, const double& value)
{

    if (name.compare("mH1p_2") == 0) {
        if (mH1p_2 >= 0.) mH1p_2 = value;
        else throw std::runtime_error("error in LeftRightSymmetricModel::SetParameter, mH1p_2 < 0!");
    } else if (name.compare("mH2p_2") == 0) {
        if (mH2p_2 >= 0.) mH2p_2 = value;
        else throw std::runtime_error("error in LeftRightSymmetricModel::SetParameter, mH2p_2 < 0!");
    } else if (name.compare("mdeltappR_2") == 0) {
        if (mdeltappR_2 >= 0.) mdeltappR_2 = value;
        else throw std::runtime_error("error in LeftRightSymmetricModel::SetParameter, mdeltappR_2 < 0!");
    } else if (name.compare("xi_LRSM") == 0) {
        if (xi_LRSM >= 0. && xi_LRSM < 1.) xi_LRSM = value;
        else throw std::runtime_error("error in LeftRightSymmetricModel::SetParameter, xi_LRSM < 0!");
    } else if (name.compare("mWR") == 0) {
        if (mWR >= 0.) mWR = value;
        else throw std::runtime_error("error in LeftRightSymmetricModel::SetParameter, mWR < 0!");
    } else if (name.compare("lambda1_LRSM") == 0)
        lambda1_LRSM = value;
    else if (name.compare("lambda2_LRSM") == 0)
        lambda2_LRSM = value;
    else if (name.compare("lambda3_LRSM") == 0)
        lambda3_LRSM = value;
    else if (name.compare("lambda4_LRSM") == 0)
        lambda4_LRSM = value;
    else if (name.compare("rho1_LRSM") == 0)
        rho1_LRSM = value;
    else if (name.compare("rho4_LRSM") == 0)
        rho4_LRSM = value;
    else if (name.compare("alpha1_LRSM") == 0)
        alpha1_LRSM = value;
    else if (name.compare("alpha2_LRSM") == 0)
        alpha2_LRSM = value;
    else if (name.compare("alpha_LRSM") == 0)
        alpha_LRSM = value;
    else if (name.compare("delta2_LRSM") == 0)
        delta2_LRSM = value;
    else if (name.compare("thetaR12_LRSM") == 0)
        thetaR12_LRSM = value;
    else if (name.compare("thetaR13_LRSM") == 0)
        thetaR13_LRSM = value;
    else if (name.compare("thetaR23_LRSM") == 0)
        thetaR23_LRSM = value;
    else if (name.compare("phiR1_LRSM") == 0)
        phiR1_LRSM = value;
    else if (name.compare("phiR2_LRSM") == 0)
        phiR2_LRSM = value;
    else if (name.compare("phiR3_LRSM") == 0)
        phiR3_LRSM = value;
    else if (name.compare("phiR4_LRSM") == 0)
        phiR4_LRSM = value;
    else if (name.compare("phiR5_LRSM") == 0)
        phiR5_LRSM = value;
    else if (name.compare("phiR6_LRSM") == 0)
        phiR6_LRSM = value;
    else if (name.compare("Q_LRSM") == 0) {
        if (Q_LRSM >= 0.) Q_LRSM = value;
        else throw std::runtime_error("error in LeftRightSymmetricModel::SetParameter, Q_LRSM < 0!");
    } else
        StandardModel::setParameter(name, value);
}

bool LeftRightSymmetricModel::CheckParameters(const std::map<std::string, double>& DPars) 
{
    for (int i = 0; i < NLeftRightSymmetricModelvars; i++) {
        if (DPars.find(LeftRightSymmetricModelvars[i]) == DPars.end()) {
            std::cout << "missing mandatory GeneralTHDM parameter " << LeftRightSymmetricModelvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(LeftRightSymmetricModelvars[i]);
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

/*----------------------------------------------------------------------------*/
/* Flags */

bool LeftRightSymmetricModel::setFlagStr(const std::string name, const std::string value) 
{
    bool res = false;
//    if (name.compare("modelTypeflag") == 0) {
//        if (checkmodelType(value)) {
//            flag_model = value;
//            res = true;
//        } else {
//            throw std::runtime_error("LeftRightSymmetricModel::setFlagStr(): Invalid flag "
//                    + name + "=" + value);
//        }
//    } else {
        res = StandardModel::setFlagStr(name, value);
//    }

    return (res);
}

bool LeftRightSymmetricModel::setFlag(const std::string name, const bool value) 
{
    bool res = false;
    if (name.compare("CPviolation") == 0) {
        flag_CPV = value;
        res = true;
    } else
        res = StandardModel::setFlag(name, value);

    return (res);
}
