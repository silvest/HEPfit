/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include <StandardModelMatching.h>
#include "THDM.h"

const std::string THDM::THDMvars[NTHDMvars] = {"logtb","bma","mHh2","mHh2mmA2","mHh2mmHp2","m12_2","lambda6","lambda7",
                                               "BDtaunu_SM","BDtaunu_A","BDtaunu_B","BDstartaunu_SM","BDstartaunu_A","BDstartaunu_B",
                                               "BHatBsTHDM","etaBsTHDM", "bsgamma_theoryerror","Q_THDM"};

THDM::THDM() : StandardModel() {
//    mycache = new THDMcache();
    
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bma", boost::cref(bma)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh2", boost::cref(mHh2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh2mmA2", boost::cref(mHh2mmA2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHh2mmHp2", boost::cref(mHh2mmHp2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("m12_2", boost::cref(m12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambda6", boost::cref(lambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("lambda7", boost::cref(lambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDtaunu_SM", boost::cref(BDtaunu_SM)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDtaunu_A", boost::cref(BDtaunu_A)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDtaunu_B", boost::cref(BDtaunu_B)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDstartaunu_SM", boost::cref(BDstartaunu_SM)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDstartaunu_A", boost::cref(BDstartaunu_A)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BDstartaunu_B", boost::cref(BDstartaunu_B)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("BHatBsTHDM", boost::cref(BHatBsTHDM)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("etaBsTHDM", boost::cref(etaBsTHDM)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("bsgamma_theoryerror", boost::cref(bsgamma_theoryerror)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Q_THDM", boost::cref(Q_THDM)));
}

THDM::~THDM(){
    if (IsModelInitialized()) {
            if (myTHDMMatching != NULL) delete(myTHDMMatching);
//            if (mycache != NULL) delete(mycache);
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
}

bool THDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);
    
    //In THDM U couple with v2, D with v1 and L with v1
//    if (requireYu || requireCKM) {
//        Yu = matrix<complex>::Id(3);
//        for (int i = 0; i < 3; i++)
//            Yu.assign(i, i, this->quarks[UP + 2 * i].getMass() / (v() * sinb) * sqrt(2.));
//        Yu = VCKM.transpose()*Yu;
//    }
//    if (requireYd) {
//        for (int i = 0; i < 3; i++)
//            Yd.assign(i, i, this->QCD::quarks[DOWN + 2 * i].getMass() / (v() * cosb) * sqrt(2.));
//    }
//    if (requireYe) {
//        for (int i = 0; i < 3; i++)
//            Ye.assign(i, i, this->leptons[ELECTRON + 2 * i].getMass() / (v() * cosb) * sqrt(2.));
//    }
//    if (requireYn) {
//        Yn = matrix<complex>::Id(3);
//        for (int i = 0; i < 3; i++)
//            Yn.assign(i, i, this->leptons[NEUTRINO_1 + 2 * i].getMass() / (v() * cosb) * sqrt(2.));
//        Yn = Yn * UPMNS.hconjugate();
//    }

    /* Necessary for updating StandardModel parameters in StandardModelMatching,
     * and THDM and THDM-derived parameters in THDMMatching */
    myTHDMMatching->StandardModelMatching::updateSMParameters();
    myTHDMMatching->updateTHDMParameters();

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
    else if(name.compare("mHh2") == 0)
        mHh2 = value;
    else if(name.compare("mHh2mmA2") == 0)
        mHh2mmA2 = value;
    else if(name.compare("mHh2mmHp2") == 0)
        mHh2mmHp2 = value;
    else if(name.compare("m12_2") == 0)
        m12_2 = value;
    else if(name.compare("lambda6") == 0)
        lambda6 = value;
    else if(name.compare("lambda7") == 0)
        lambda7 = value;
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
    else if(name.compare("BHatBsTHDM") == 0)
        BHatBsTHDM = value;
    else if(name.compare("etaBsTHDM") == 0)
        etaBsTHDM = value;
    else if(name.compare("bsgamma_theoryerror") == 0)
        bsgamma_theoryerror = value;
    else if(name.compare("Q_THDM") == 0)
        Q_THDM = value;
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

bool THDM::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if(name.compare("modelTypeflag") == 0)
    {
        flag_model = value;
        res = true;
    }
    else
    {
        res = StandardModel::setFlagStr(name,value);
    }

    return(res);
}
//
//double THDM::computeCosa() const
//{
//    return cos(atan(pow(10.,logtb))-bma);
//}
//
//double THDM::computeSina() const
//{
//    return sin(atan(pow(10.,logtb))-bma);
//}
//
//double THDM::v1() const
//{
//    return v() * cosb;
//}
//
//double THDM::v2() const
//{
//    return v() * sinb;
//}
