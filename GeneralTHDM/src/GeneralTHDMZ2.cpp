/* 
 * Copyright (C) 2025 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDMZ2.h"

const std::string GeneralTHDMZ2::GeneralTHDMZ2vars[NGeneralTHDMZ2vars] = {"tanb", "bma", "m12_2"};

GeneralTHDMZ2::GeneralTHDMZ2() : GeneralTHDM()
{
    ModelParamMap.insert(std::make_pair("tanb",  std::cref(tanb)));
    ModelParamMap.insert(std::make_pair("bma",   std::cref(bma)));
    ModelParamMap.insert(std::make_pair("m12_2", std::cref(m12_2)));
}

GeneralTHDMZ2::~GeneralTHDMZ2()
{
    if (IsModelInitialized()) {
    }
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool GeneralTHDMZ2::InitializeModel()
{
    setModelInitialized(GeneralTHDM::InitializeModel());
    return (true);
}

bool GeneralTHDMZ2::PreUpdate()
{
    if (!GeneralTHDM::PreUpdate()) return (false);
    
    if(!getCPconservationflag())
        throw std::runtime_error("ERROR: GeneralTHDMZ2 class requires CP conservation, check corresponding flag");

    return (true);
}

bool GeneralTHDMZ2::Update(const std::map<std::string, double>& DPars)
{
    if (!PreUpdate()) return (false);

    if (!GeneralTHDM::Update(DPars)) return (false);

    UpdateError = false;

    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);

    if (UpdateError) return (false);

    if (!PostUpdate()) return (false);

    return (true);
}

bool GeneralTHDMZ2::PostUpdate()
{
    if (!GeneralTHDM::PostUpdate()) return (false);

    return (true);
}

void GeneralTHDMZ2::setParameter(const std::string name, const double& value)
{
    if (name.compare("tanb") == 0) {
        tanb = value;
        if(tanb > 0.) {
            beta  = atan(tanb);
            cosb  = cos(beta);
            cos2b = cos(2.*beta);
            cos4b = cos(4.*beta);
            cos6b = cos(6.*beta);
            sinb  = sin(beta);
            sin2b = sin(2.*beta);
            sin6b = sin(6.*beta);
        }
        else {
            throw std::runtime_error("ERROR: in GeneralTHDMZ2::setParameter, tanb < 0!");
        }
    }
    else if (name.compare("bma") == 0) {
        bma = value;
        cos2bma = cos(2.*bma);
        sin2bma = sin(2.*bma);
    }
    else if (name.compare("m12_2") == 0)
        m12_2 = value;
    else if (name.compare("alpha1") == 0)
        GeneralTHDM::setParameter("alpha1", beta - bma);
    else if (name.compare("lambda2") == 0) {
        mh_2 = getmH1sq();
        mH_2 = getmH2sq();
        mA_2 = getmH3sq();
        mHp_2 = getmHp2();
        vev = GeneralTHDM::v();
        GeneralTHDM::setParameter("lambda2", ((4.*(mh_2 - mH_2)*(3. + 5.*cos4b)*cos2bma/cosb/sinb +
        (-32.*m12_2*(1. + cos4b) + (mh_2 + mH_2)*(7.*sin2b + 3.*sin6b))/cosb/cosb/sinb/sinb +
        (mH_2 - mh_2)*(20. - cos2b + 12.*cos4b + cos6b)*2.*sin2bma)/64./cosb/sinb/vev/vev));
    }
    else if (name.compare("lambda3") == 0)
        GeneralTHDM::setParameter("lambda3", ((-32.*m12_2 + 8.*(mh_2 + mH_2 + mA_2)*sin2b +
        (mh_2 - mH_2)*(-4. + 3.*cos2b + 4.*cos4b + cos6b)*sin2bma)/16./cosb/sinb/vev/vev));
    else if (name.compare("Relambda7") == 0)
        GeneralTHDM::setParameter("Relambda7", 0.);
    else if (name.compare("Nu_11r") == 0) {
        if (flag_model == "type1")
            GeneralTHDM::setParameter("Nu_11r", 1./tanb);
        else if (flag_model == "type2")
            GeneralTHDM::setParameter("Nu_11r", 1./tanb);
        else if (flag_model == "typeX")
            GeneralTHDM::setParameter("Nu_11r", 1./tanb);
        else if (flag_model == "typeY")
            GeneralTHDM::setParameter("Nu_11r", 1./tanb);
    }
    else if (name.compare("Nd_11r") == 0) {
        if (flag_model == "type1")
            GeneralTHDM::setParameter("Nd_11r", 1./tanb);
        else if (flag_model == "type2")
            GeneralTHDM::setParameter("Nd_11r", -tanb);
        else if (flag_model == "typeX")
            GeneralTHDM::setParameter("Nd_11r", 1./tanb);
        else if (flag_model == "typeY")
            GeneralTHDM::setParameter("Nd_11r", -tanb);
    }
    else if (name.compare("Nl_11r") == 0) {
        if (flag_model == "type1")
            GeneralTHDM::setParameter("Nl_11r", 1./tanb);
        else if (flag_model == "type2")
            GeneralTHDM::setParameter("Nl_11r", -tanb);
        else if (flag_model == "typeX")
            GeneralTHDM::setParameter("Nl_11r", -tanb);
        else if (flag_model == "typeY")
            GeneralTHDM::setParameter("Nl_11r", 1./tanb);
    }
    else
        GeneralTHDM::setParameter(name, value);
}

bool GeneralTHDMZ2::CheckParameters(const std::map<std::string, double>& DPars)
{
    for (int i = 0; i < NGeneralTHDMZ2vars; i++) {
        if (DPars.find(GeneralTHDMZ2vars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory GeneralTHDMZ2 parameter " << GeneralTHDMZ2vars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(GeneralTHDMZ2vars[i]);
        }
    }
    return(GeneralTHDM::CheckParameters(DPars));
}


///////////////////////////////////////////////////////////////////////////
// Flags

bool GeneralTHDMZ2::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;

    if(name.compare("modelType") == 0) {
        if (CheckModelType(value))
        {
            flag_model = value;
            res = true;
        }
        else
        {
            throw std::runtime_error("GeneralTHDMZ2::setFlagStr(): invalid model-type flag"
                + name + "=" + value);
        }
    }
    else
        res = GeneralTHDM::setFlagStr(name,value);

    return res;
}

bool GeneralTHDMZ2::CheckModelType(const std::string modeltype) const
{
    if (modeltype.compare("type1") == 0 ||
        modeltype.compare("type2") == 0 ||
        modeltype.compare("typeX") == 0 ||
        modeltype.compare("typeY") == 0 )
        return true;
    else
        return false;
}
