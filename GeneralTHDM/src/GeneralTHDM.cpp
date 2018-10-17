/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"

const std::string GeneralTHDM::GeneralTHDMvars[NGeneralTHDMvars] = {"logtb", "mHp2", "mH2sq", "mH3sq", "alpha1", "alpha2", "alpha3", "Relambda5", "Imlambda5", "Relambda6", "Relambda7",
"Nu_11r", "Nu_11i", "Nu_12r", "Nu_12i", "Nu_13r", "Nu_13i", 
"Nu_21r", "Nu_21i", "Nu_22r", "Nu_22i", "Nu_23r", "Nu_23i", 
"Nu_31r", "Nu_31i", "Nu_32r", "Nu_32i", "Nu_33r", "Nu_33i", 
"Nd_11r", "Nd_11i", "Nd_12r", "Nd_12i", "Nd_13r", "Nd_13i", 
"Nd_21r", "Nd_21i", "Nd_22r", "Nd_22i", "Nd_23r", "Nd_23i", 
"Nd_31r", "Nd_31i", "Nd_32r", "Nd_32i", "Nd_33r", "Nd_33i", 
"Nl_11r", "Nl_11i", "Nl_12r", "Nl_12i", "Nl_13r", "Nl_13i", 
"Nl_21r", "Nl_21i", "Nl_22r", "Nl_22i", "Nl_23r", "Nl_23i", 
"Nl_31r", "Nl_31i", "Nl_32r", "Nl_32i", "Nl_33r", "Nl_33i", 
"Q_GTHDM"};

GeneralTHDM::GeneralTHDM() : StandardModel(), GTHDMM(*this) {

    SMM.setObj((StandardModelMatching&) GTHDMM.getObj());
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("logtb", boost::cref(logtb)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mHp2", boost::cref(mHp2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mH2sq", boost::cref(mH2sq)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("mH3sq", boost::cref(mH3sq)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha1", boost::cref(alpha1)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha2", boost::cref(alpha2)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("alpha3", boost::cref(alpha3)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda5", boost::cref(Relambda5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Imlambda5", boost::cref(Imlambda5)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda6", boost::cref(Relambda6)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Relambda7", boost::cref(Relambda7)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_11r", boost::cref(Nu_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_11i", boost::cref(Nu_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_12r", boost::cref(Nu_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_12i", boost::cref(Nu_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_13r", boost::cref(Nu_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_13i", boost::cref(Nu_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_21r", boost::cref(Nu_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_21i", boost::cref(Nu_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_22r", boost::cref(Nu_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_22i", boost::cref(Nu_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_23r", boost::cref(Nu_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_23i", boost::cref(Nu_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_31r", boost::cref(Nu_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_31i", boost::cref(Nu_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_32r", boost::cref(Nu_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_32i", boost::cref(Nu_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_33r", boost::cref(Nu_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nu_33i", boost::cref(Nu_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_11r", boost::cref(Nd_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_11i", boost::cref(Nd_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_12r", boost::cref(Nd_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_12i", boost::cref(Nd_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_13r", boost::cref(Nd_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_13i", boost::cref(Nd_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_21r", boost::cref(Nd_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_21i", boost::cref(Nd_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_22r", boost::cref(Nd_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_22i", boost::cref(Nd_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_23r", boost::cref(Nd_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_23i", boost::cref(Nd_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_31r", boost::cref(Nd_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_31i", boost::cref(Nd_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_32r", boost::cref(Nd_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_32i", boost::cref(Nd_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_33r", boost::cref(Nd_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nd_33i", boost::cref(Nd_33i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_11r", boost::cref(Nl_11r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_11i", boost::cref(Nl_11i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_12r", boost::cref(Nl_12r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_12i", boost::cref(Nl_12i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_13r", boost::cref(Nl_13r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_13i", boost::cref(Nl_13i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_21r", boost::cref(Nl_21r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_21i", boost::cref(Nl_21i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_22r", boost::cref(Nl_22r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_22i", boost::cref(Nl_22i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_23r", boost::cref(Nl_23r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_23i", boost::cref(Nl_23i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_31r", boost::cref(Nl_31r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_31i", boost::cref(Nl_31i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_32r", boost::cref(Nl_32r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_32i", boost::cref(Nl_32i)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_33r", boost::cref(Nl_33r)));
    ModelParamMap.insert(std::pair<std::string, boost::reference_wrapper<const double> >("Nl_33i", boost::cref(Nl_33i)));
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
    else if(name.compare("mH2sq") == 0)
        mH2sq = value;
    else if(name.compare("mH3sq") == 0)
        mH3sq = value;
    else if(name.compare("alpha1") == 0) {
        alpha1 = value;
        cosalpha1 = cos(alpha1);
        sinalpha1 = sin(alpha1);
      // std::cout<<"alpha1 before ordering = "<<alpha1<<std::endl;

    }
    else if(name.compare("alpha2") == 0) {
        alpha2 = value;
        cosalpha2 = cos(alpha2);
        sinalpha2 = sin(alpha2);
              //std::cout<<"alpha2 before ordering = "<<alpha2<<std::endl;
    }
    else if(name.compare("alpha3") == 0) {
        alpha3 = value;
        cosalpha3 = cos(alpha3);
        sinalpha3 = sin(alpha3);
              //std::cout<<"alpha3 before ordering = "<<alpha3<<std::endl;
    }
    else if(name.compare("Relambda5") == 0)
        Relambda5 = value;
    else if(name.compare("Imlambda5") == 0)
        Imlambda5 = value;
    else if(name.compare("Relambda6") == 0)
        Relambda6 = value;
    else if(name.compare("Relambda7") == 0)
        Relambda7 = value;
    else if(name.compare("Nu_11r") == 0)
        Nu_11r = value;
    else if(name.compare("Nu_11i") == 0)
        Nu_11i = value;
    else if(name.compare("Nu_12r") == 0)
        Nu_12r = value;
    else if(name.compare("Nu_12i") == 0)
        Nu_12i = value;
    else if(name.compare("Nu_13r") == 0)
        Nu_13r = value;
    else if(name.compare("Nu_13i") == 0)
        Nu_13i = value;
    else if(name.compare("Nu_21r") == 0)
        Nu_21r = value;
    else if(name.compare("Nu_21i") == 0)
        Nu_21i = value;
    else if(name.compare("Nu_22r") == 0)
        Nu_22r = value;
    else if(name.compare("Nu_22i") == 0)
        Nu_22i = value;
    else if(name.compare("Nu_23r") == 0)
        Nu_23r = value;
    else if(name.compare("Nu_23i") == 0)
        Nu_23i = value;
    else if(name.compare("Nu_31r") == 0)
        Nu_31r = value;
    else if(name.compare("Nu_31i") == 0)
        Nu_31i = value;
    else if(name.compare("Nu_32r") == 0)
        Nu_32r = value;
    else if(name.compare("Nu_32i") == 0)
        Nu_32i = value;
    else if(name.compare("Nu_33r") == 0)
        Nu_33r = value;
    else if(name.compare("Nu_33i") == 0)
        Nu_33i = value;
    else if(name.compare("Nd_11r") == 0)
        Nd_11r = value;
    else if(name.compare("Nd_11i") == 0)
        Nd_11i = value;
    else if(name.compare("Nd_12r") == 0)
        Nd_12r = value;
    else if(name.compare("Nd_12i") == 0)
        Nd_12i = value;
    else if(name.compare("Nd_13r") == 0)
        Nd_13r = value;
    else if(name.compare("Nd_13i") == 0)
        Nd_13i = value;
    else if(name.compare("Nd_21r") == 0)
        Nd_21r = value;
    else if(name.compare("Nd_21i") == 0)
        Nd_21i = value;
    else if(name.compare("Nd_22r") == 0)
        Nd_22r = value;
    else if(name.compare("Nd_22i") == 0)
        Nd_22i = value;
    else if(name.compare("Nd_23r") == 0)
        Nd_23r = value;
    else if(name.compare("Nd_23i") == 0)
        Nd_23i = value;
    else if(name.compare("Nd_31r") == 0)
        Nd_31r = value;
    else if(name.compare("Nd_31i") == 0)
        Nd_31i = value;
    else if(name.compare("Nd_32r") == 0)
        Nd_32r = value;
    else if(name.compare("Nd_32i") == 0)
        Nd_32i = value;
    else if(name.compare("Nd_33r") == 0)
        Nd_33r = value;
    else if(name.compare("Nd_33i") == 0)
        Nd_33i = value;
    else if(name.compare("Nl_11r") == 0)
        Nl_11r = value;
    else if(name.compare("Nl_11i") == 0)
        Nl_11i = value;
    else if(name.compare("Nl_12r") == 0)
        Nl_12r = value;
    else if(name.compare("Nl_12i") == 0)
        Nl_12i = value;
    else if(name.compare("Nl_13r") == 0)
        Nl_13r = value;
    else if(name.compare("Nl_13i") == 0)
        Nl_13i = value;
    else if(name.compare("Nl_21r") == 0)
        Nl_21r = value;
    else if(name.compare("Nl_21i") == 0)
        Nl_21i = value;
    else if(name.compare("Nl_22r") == 0)
        Nl_22r = value;
    else if(name.compare("Nl_22i") == 0)
        Nl_22i = value;
    else if(name.compare("Nl_23r") == 0)
        Nl_23r = value;
    else if(name.compare("Nl_23i") == 0)
        Nl_23i = value;
    else if(name.compare("Nl_31r") == 0)
        Nl_31r = value;
    else if(name.compare("Nl_31i") == 0)
        Nl_31i = value;
    else if(name.compare("Nl_32r") == 0)
        Nl_32r = value;
    else if(name.compare("Nl_32i") == 0)
        Nl_32i = value;
    else if(name.compare("Nl_33r") == 0)
        Nl_33r = value;
    else if(name.compare("Nl_33i") == 0)
        Nl_33i = value;
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

//bool GeneralTHDM::setFlagStr(const std::string name, const std::string value)
//{
//    bool res = false;
//    if(name.compare("ATHDMflag") == 0)
//    {
//        if (checkmodelType(value))
//        {
//            flag_model = value;
//            res = true;
//        }
//        else
//        {
//            throw std::runtime_error("GeneralTHDM::setFlagStr(): Invalid flag "
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

bool GeneralTHDM::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if(name.compare("ATHDMflag") == 0) {
    std::cout<<"ATHDMflag = "<< value<<std::endl;
        flag_ATHDM = value;
        res = true;
    }
    else if(name.compare("CPconservation") == 0) {
    std::cout<<"CPconservation = "<< value<<std::endl;
        flag_CPconservation = value;
        res = true;
    }
    else
        res = StandardModel::setFlag(name,value);

    return(res);
}

double GeneralTHDM::muggH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_gg;
}

double GeneralTHDM::muVBF(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::mueeWBF(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::muWH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::muZH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::mueeZH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::muVH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::muVBFpVH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV;
}

double GeneralTHDM::muttH(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t);
}

double GeneralTHDM::computeGammaTotalRatio() const
{
    return getMyGTHDMCache()->sumModBRs;
}

double GeneralTHDM::GammaTotal() const
{
    return getMyGTHDMCache()->Gamma_h;
}

double GeneralTHDM::BrHggRatio() const
{
    return getMyGTHDMCache()->rh_gg / computeGammaTotalRatio();
}

double GeneralTHDM::BrHWWRatio() const
{
    return getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::BrHZZRatio() const
{
    return getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::BrHZgaRatio() const
{
    return getMyGTHDMCache()->rh_Zga / computeGammaTotalRatio();
}

double GeneralTHDM::BrHgagaRatio() const
{
    return getMyGTHDMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeneralTHDM::BrHmumuRatio() const
{
    return getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_mu*getMyGTHDMCache()->beta_h_mu) / computeGammaTotalRatio();
}

double GeneralTHDM::BrHtautauRatio() const
{
    return getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_tau*getMyGTHDMCache()->beta_h_tau) / computeGammaTotalRatio();
}

double GeneralTHDM::BrHccRatio() const
{
    return getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_c*getMyGTHDMCache()->beta_h_c) / computeGammaTotalRatio();
}

double GeneralTHDM::BrHbbRatio() const
{
    return getMyGTHDMCache()->rh_QdQdE + getMyGTHDMCache()->rh_QdQdO/(getMyGTHDMCache()->beta_h_b*getMyGTHDMCache()->beta_h_b) / computeGammaTotalRatio();
}

double GeneralTHDM::muggHgaga(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_gg * getMyGTHDMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeneralTHDM::muVBFHgaga(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * getMyGTHDMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeneralTHDM::muVHgaga(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * getMyGTHDMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeneralTHDM::muttHgaga(const double sqrt_s) const
{
    return (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t)) * getMyGTHDMCache()->rh_gaga / computeGammaTotalRatio();
}

double GeneralTHDM::muggHZZ(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_gg * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muVBFHZZ(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muVHZZ(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muttHZZ(const double sqrt_s) const
{
    return (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t)) * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muggHWW(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_gg * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muVBFHWW(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muVHWW(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muttHWW(const double sqrt_s) const
{
    return (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t)) * getMyGTHDMCache()->rh_VV / computeGammaTotalRatio();
}

double GeneralTHDM::muggHtautau(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_gg * (getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_tau*getMyGTHDMCache()->beta_h_tau)) / computeGammaTotalRatio();
}

double GeneralTHDM::muVBFHtautau(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * (getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_tau*getMyGTHDMCache()->beta_h_tau)) / computeGammaTotalRatio();
}

double GeneralTHDM::muVHtautau(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * (getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_tau*getMyGTHDMCache()->beta_h_tau)) / computeGammaTotalRatio();
}

double GeneralTHDM::muttHtautau(const double sqrt_s) const
{
    return (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t))
            * (getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_tau*getMyGTHDMCache()->beta_h_tau)) / computeGammaTotalRatio();
}

double GeneralTHDM::muggHbb(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_gg * (getMyGTHDMCache()->rh_QdQdE + getMyGTHDMCache()->rh_QdQdO/(getMyGTHDMCache()->beta_h_b*getMyGTHDMCache()->beta_h_b)) / computeGammaTotalRatio();
}

double GeneralTHDM::muVBFHbb(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * (getMyGTHDMCache()->rh_QdQdE + getMyGTHDMCache()->rh_QdQdO/(getMyGTHDMCache()->beta_h_b*getMyGTHDMCache()->beta_h_b)) / computeGammaTotalRatio();
}

double GeneralTHDM::muVHbb(const double sqrt_s) const
{
    return getMyGTHDMCache()->rh_VV * (getMyGTHDMCache()->rh_QdQdE + getMyGTHDMCache()->rh_QdQdO/(getMyGTHDMCache()->beta_h_b*getMyGTHDMCache()->beta_h_b)) / computeGammaTotalRatio();
}

double GeneralTHDM::muttHbb(const double sqrt_s) const
{
    return (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t))
            * (getMyGTHDMCache()->rh_QdQdE + getMyGTHDMCache()->rh_QdQdO/(getMyGTHDMCache()->beta_h_b*getMyGTHDMCache()->beta_h_b)) / computeGammaTotalRatio();
}

double GeneralTHDM::muppHmumu(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * getMyGTHDMCache()->rh_gg + 0.122 * getMyGTHDMCache()->rh_VV + 0.006 * (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t)))
                * (getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_mu*getMyGTHDMCache()->beta_h_mu)) / computeGammaTotalRatio();
    }
    if(sqrt_s==13)
    {
        return (0.871 * getMyGTHDMCache()->rh_gg + 0.119 * getMyGTHDMCache()->rh_VV + 0.010 * (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t)))
                * (getMyGTHDMCache()->rh_QlQlE + getMyGTHDMCache()->rh_QlQlO/(getMyGTHDMCache()->beta_h_mu*getMyGTHDMCache()->beta_h_mu)) / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHmumu is only defined for 8 or 13 TeV.");
    }
}

double GeneralTHDM::muppHZga(const double sqrt_s) const
{
    if(sqrt_s==8)
    {
        return (0.872 * getMyGTHDMCache()->rh_gg + 0.122 * getMyGTHDMCache()->rh_VV + 0.006 * (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t))) * getMyGTHDMCache()->rh_Zga / computeGammaTotalRatio();
    }
    if(sqrt_s==13)
    {
        return (0.871 * getMyGTHDMCache()->rh_gg + 0.119 * getMyGTHDMCache()->rh_VV + 0.010 * (getMyGTHDMCache()->rh_QuQuE + getMyGTHDMCache()->rh_QuQuO/(getMyGTHDMCache()->beta_h_t*getMyGTHDMCache()->beta_h_t))) * getMyGTHDMCache()->rh_Zga / computeGammaTotalRatio();
    }
    else
    {
        throw std::runtime_error("The observable muppHZga is only defined for 8 or 13 TeV.");
    }
}

double GeneralTHDM::Mw() const{
    double MZ = StandardModel::Mz;
    return ( MZ / sqrt(2.0) * sqrt(1.0 + sqrt(1.0 - 4.0 * M_PI * StandardModel::ale / (sqrt(2.0) * StandardModel::GF * MZ* MZ))));
}
