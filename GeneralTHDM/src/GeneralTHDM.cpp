/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

const std::string GeneralTHDM::GeneralTHDMvars[NGeneralTHDMvars] = {"logtb", "mHp2", "alpha1", "alpha2", "alpha3", "Rem12_2", "Imlambda5", "Relambda6", "Imlambda6", "Relambda7", "Imlambda7", 
"Ytu_11r", "Ytu_11i", "Ytu_12r", "Ytu_12i", "Ytu_13r", "Ytu_13i", 
"Ytu_21r", "Ytu_21i", "Ytu_22r", "Ytu_22i", "Ytu_23r", "Ytu_23i", 
"Ytu_31r", "Ytu_31i", "Ytu_32r", "Ytu_32i", "Ytu_33r", "Ytu_33i", 
"Ytd_11r", "Ytd_11i", "Ytd_12r", "Ytd_12i", "Ytd_13r", "Ytd_13i", 
"Ytd_21r", "Ytd_21i", "Ytd_22r", "Ytd_22i", "Ytd_23r", "Ytd_23i", 
"Ytd_31r", "Ytd_31i", "Ytd_32r", "Ytd_32i", "Ytd_33r", "Ytd_33i", 
"Ytl_11r", "Ytl_11i", "Ytl_12r", "Ytl_12i", "Ytl_13r", "Ytl_13i", 
"Ytl_21r", "Ytl_21i", "Ytl_22r", "Ytl_22i", "Ytl_23r", "Ytl_23i", 
"Ytl_31r", "Ytl_31i", "Ytl_32r", "Ytl_32i", "Ytl_33r", "Ytl_33i", 
"Q_GTHDM"};

GeneralTHDM::GeneralTHDM() : StandardModel()/*, GTHDMM(*this)*/ {

//    SMM.setObj((StandardModelMatching&) THDMM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHp2", boost::cref(mHp2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1", boost::cref(alpha1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2", boost::cref(alpha2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha3", boost::cref(alpha3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Rem12_2", boost::cref(Rem12_2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Imlambda5", boost::cref(Imlambda5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda6", boost::cref(Relambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Imlambda6", boost::cref(Imlambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda7", boost::cref(Relambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Imlambda7", boost::cref(Imlambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_11r", boost::cref(Ytu_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_11i", boost::cref(Ytu_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_12r", boost::cref(Ytu_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_12i", boost::cref(Ytu_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_13r", boost::cref(Ytu_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_13i", boost::cref(Ytu_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_21r", boost::cref(Ytu_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_21i", boost::cref(Ytu_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_22r", boost::cref(Ytu_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_22i", boost::cref(Ytu_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_23r", boost::cref(Ytu_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_23i", boost::cref(Ytu_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_31r", boost::cref(Ytu_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_31i", boost::cref(Ytu_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_32r", boost::cref(Ytu_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_32i", boost::cref(Ytu_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_33r", boost::cref(Ytu_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytu_33i", boost::cref(Ytu_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_11r", boost::cref(Ytd_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_11i", boost::cref(Ytd_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_12r", boost::cref(Ytd_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_12i", boost::cref(Ytd_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_13r", boost::cref(Ytd_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_13i", boost::cref(Ytd_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_21r", boost::cref(Ytd_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_21i", boost::cref(Ytd_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_22r", boost::cref(Ytd_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_22i", boost::cref(Ytd_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_23r", boost::cref(Ytd_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_23i", boost::cref(Ytd_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_31r", boost::cref(Ytd_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_31i", boost::cref(Ytd_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_32r", boost::cref(Ytd_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_32i", boost::cref(Ytd_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_33r", boost::cref(Ytd_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytd_33i", boost::cref(Ytd_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_11r", boost::cref(Ytl_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_11i", boost::cref(Ytl_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_12r", boost::cref(Ytl_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_12i", boost::cref(Ytl_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_13r", boost::cref(Ytl_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_13i", boost::cref(Ytl_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_21r", boost::cref(Ytl_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_21i", boost::cref(Ytl_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_22r", boost::cref(Ytl_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_22i", boost::cref(Ytl_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_23r", boost::cref(Ytl_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_23i", boost::cref(Ytl_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_31r", boost::cref(Ytl_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_31i", boost::cref(Ytl_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_32r", boost::cref(Ytl_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_32i", boost::cref(Ytl_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_33r", boost::cref(Ytl_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Ytl_33i", boost::cref(Ytl_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Q_GTHDM", boost::cref(Q_GTHDM)));
}

GeneralTHDM::~GeneralTHDM(){
    if (IsModelInitialized()) {
            if (myGTHDMcache != NULL) delete(myGTHDMcache);
        }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool GeneralTHDM::InitializeModel()
{
    myGTHDMcache = new GeneralTHDMcache(*this);
    setModelInitialized(StandardModel::InitializeModel());
    setModelGeneralTHDM();
    return(true);
}

bool GeneralTHDM::Init(const std::map<std::string, double>& DPars) {
    return(StandardModel::Init(DPars));
}

bool GeneralTHDM::PreUpdate()
{
    if(!StandardModel::PreUpdate()) return (false);

    return (true);
}

bool GeneralTHDM::Update(const std::map<std::string, double>& DPars) {

    if(!PreUpdate()) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if(!PostUpdate()) return (false);

    return (true);
}

bool GeneralTHDM::PostUpdate()
{
    if(!StandardModel::PostUpdate()) return (false);

    myGTHDMcache->updateCache();

    return (true);
}

void GeneralTHDM::setParameter(const std::string name, const double& value){

    if(name.compare("logtb") == 0) {
        logtb = value;
        tanb = pow(10.,logtb);
        if(tanb > 0.) {
            sinb = tanb / sqrt(1. + tanb*tanb);
            cosb = 1. / sqrt(1. + tanb*tanb);
        }
        else {
            throw std::runtime_error("error in GeneralTHDM::SetParameter, tanb < 0!");
          }
        }
    else if(name.compare("mHp2") == 0)
        mHp2 = value;
    else if(name.compare("alpha1") == 0) {
        alpha1 = value;
        cosalpha1 = cos(alpha1);
        sinalpha1 = sin(alpha1);
    }
    else if(name.compare("alpha2") == 0) {
        alpha2 = value;
        cosalpha2 = cos(alpha2);
        sinalpha2 = sin(alpha2);
    }
    else if(name.compare("alpha3") == 0) {
        alpha3 = value;
        cosalpha3 = cos(alpha3);
        sinalpha3 = sin(alpha3);
    }
    else if(name.compare("Rem12_2") == 0)
        Rem12_2 = value;
    else if(name.compare("Imlambda5") == 0)
        Imlambda5 = value;
    else if(name.compare("Relambda6") == 0)
        Relambda6 = value;
    else if(name.compare("Imlambda6") == 0)
        Imlambda6 = value;
    else if(name.compare("Relambda7") == 0)
        Relambda7 = value;
    else if(name.compare("Imlambda7") == 0)
        Imlambda7 = value;
    else if(name.compare("Ytu_11r") == 0)
        Ytu_11r = value;
    else if(name.compare("Ytu_11i") == 0)
        Ytu_11i = value;
    else if(name.compare("Ytu_12r") == 0)
        Ytu_12r = value;
    else if(name.compare("Ytu_12i") == 0)
        Ytu_12i = value;
    else if(name.compare("Ytu_13r") == 0)
        Ytu_13r = value;
    else if(name.compare("Ytu_13i") == 0)
        Ytu_13i = value;
    else if(name.compare("Ytu_21r") == 0)
        Ytu_21r = value;
    else if(name.compare("Ytu_21i") == 0)
        Ytu_21i = value;
    else if(name.compare("Ytu_22r") == 0)
        Ytu_22r = value;
    else if(name.compare("Ytu_22i") == 0)
        Ytu_22i = value;
    else if(name.compare("Ytu_23r") == 0)
        Ytu_23r = value;
    else if(name.compare("Ytu_23i") == 0)
        Ytu_23i = value;
    else if(name.compare("Ytu_31r") == 0)
        Ytu_31r = value;
    else if(name.compare("Ytu_31i") == 0)
        Ytu_31i = value;
    else if(name.compare("Ytu_32r") == 0)
        Ytu_32r = value;
    else if(name.compare("Ytu_32i") == 0)
        Ytu_32i = value;
    else if(name.compare("Ytu_33r") == 0)
        Ytu_33r = value;
    else if(name.compare("Ytu_33i") == 0)
        Ytu_33i = value;
    else if(name.compare("Ytd_11r") == 0)
        Ytd_11r = value;
    else if(name.compare("Ytd_11i") == 0)
        Ytd_11i = value;
    else if(name.compare("Ytd_12r") == 0)
        Ytd_12r = value;
    else if(name.compare("Ytd_12i") == 0)
        Ytd_12i = value;
    else if(name.compare("Ytd_13r") == 0)
        Ytd_13r = value;
    else if(name.compare("Ytd_13i") == 0)
        Ytd_13i = value;
    else if(name.compare("Ytd_21r") == 0)
        Ytd_21r = value;
    else if(name.compare("Ytd_21i") == 0)
        Ytd_21i = value;
    else if(name.compare("Ytd_22r") == 0)
        Ytd_22r = value;
    else if(name.compare("Ytd_22i") == 0)
        Ytd_22i = value;
    else if(name.compare("Ytd_23r") == 0)
        Ytd_23r = value;
    else if(name.compare("Ytd_23i") == 0)
        Ytd_23i = value;
    else if(name.compare("Ytd_31r") == 0)
        Ytd_31r = value;
    else if(name.compare("Ytd_31i") == 0)
        Ytd_31i = value;
    else if(name.compare("Ytd_32r") == 0)
        Ytd_32r = value;
    else if(name.compare("Ytd_32i") == 0)
        Ytd_32i = value;
    else if(name.compare("Ytd_33r") == 0)
        Ytd_33r = value;
    else if(name.compare("Ytd_33i") == 0)
        Ytd_33i = value;
    else if(name.compare("Ytl_11r") == 0)
        Ytl_11r = value;
    else if(name.compare("Ytl_11i") == 0)
        Ytl_11i = value;
    else if(name.compare("Ytl_12r") == 0)
        Ytl_12r = value;
    else if(name.compare("Ytl_12i") == 0)
        Ytl_12i = value;
    else if(name.compare("Ytl_13r") == 0)
        Ytl_13r = value;
    else if(name.compare("Ytl_13i") == 0)
        Ytl_13i = value;
    else if(name.compare("Ytl_21r") == 0)
        Ytl_21r = value;
    else if(name.compare("Ytl_21i") == 0)
        Ytl_21i = value;
    else if(name.compare("Ytl_22r") == 0)
        Ytl_22r = value;
    else if(name.compare("Ytl_22i") == 0)
        Ytl_22i = value;
    else if(name.compare("Ytl_23r") == 0)
        Ytl_23r = value;
    else if(name.compare("Ytl_23i") == 0)
        Ytl_23i = value;
    else if(name.compare("Ytl_31r") == 0)
        Ytl_31r = value;
    else if(name.compare("Ytl_31i") == 0)
        Ytl_31i = value;
    else if(name.compare("Ytl_32r") == 0)
        Ytl_32r = value;
    else if(name.compare("Ytl_32i") == 0)
        Ytl_32i = value;
    else if(name.compare("Ytl_33r") == 0)
        Ytl_33r = value;
    else if(name.compare("Ytl_33i") == 0)
        Ytl_33i = value;
    else if(name.compare("Q_GTHDM") == 0)
        Q_GTHDM = value;
    else
        StandardModel::setParameter(name,value);
}

bool GeneralTHDM::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NGeneralTHDMvars; i++) {
        if (DPars.find(GeneralTHDMvars[i]) == DPars.end()) {
            std::cout << "missing mandatory GeneralTHDM parameter " << GeneralTHDMvars[i] << std::endl;
            return false;
        }
    }
    return(StandardModel::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags
//
//bool THDM::setFlagStr(const std::string name, const std::string value)
//{
//    bool res = false;
//    if(name.compare("modelTypeflag") == 0)
//    {
//        if (checkmodelType(value))
//        {
//            flag_model = value;
//            res = true;
//        }
//        else
//        {
//            throw std::runtime_error("THDM::setFlagStr(): Invalid flag "
//                + name + "=" + value);
//        }
//    }
//    else
//    {
//        res = StandardModel::setFlagStr(name,value);
//    }
//
//    return(res);
//}
//
//bool THDM::setFlag(const std::string name, const bool value)
//{
//    bool res = false;
//    if(name.compare("wavefunctionrenormalization") == 0) {
//        flag_wfr = value;
//        res = true;
//    }
//    else
//        res = StandardModel::setFlag(name,value);
//
//    return(res);
//}
