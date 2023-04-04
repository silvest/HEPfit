/* 
 * Copyright (C) 2016 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "GeneralTHDM.h"
#include "GeneralTHDMcache.h"
#include "GeneralTHDMSTU.h"


std::string GeneralTHDM::GeneralTHDMvars[NGeneralTHDMvars] = {"mHp2", "mH2sq", "mH3sq", "alpha1", "alpha2", "alpha3", 
"lambda2", "lambda3", "Relambda7",
"Nu_11r", "Nu_11i", "Nu_12r", "Nu_12i", "Nu_13r", "Nu_13i", 
"Nu_21r", "Nu_21i", "Nu_22r", "Nu_22i", "Nu_23r", "Nu_23i", 
"Nu_31r", "Nu_31i", "Nu_32r", "Nu_32i", "Nu_33r", "Nu_33i", 
"Nd_11r", "Nd_11i", "Nd_12r", "Nd_12i", "Nd_13r", "Nd_13i", 
"Nd_21r", "Nd_21i", "Nd_22r", "Nd_22i", "Nd_23r", "Nd_23i", 
"Nd_31r", "Nd_31i", "Nd_32r", "Nd_32i", "Nd_33r", "Nd_33i", 
"Nl_11r", "Nl_11i", "Nl_12r", "Nl_12i", "Nl_13r", "Nl_13i", 
"Nl_21r", "Nl_21i", "Nl_22r", "Nl_22i", "Nl_23r", "Nl_23i", 
"Nl_31r", "Nl_31i", "Nl_32r", "Nl_32i", "Nl_33r", "Nl_33i", 
"Q_GTHDM","RpepsGTHDM","NLOuniscaleGTHDM"};

GeneralTHDM::GeneralTHDM() : NPbase(), GTHDMM(*this) {

    SMM.setObj((StandardModelMatching&) GTHDMM.getObj());
    //ModelParamMap.insert(std::make_pair("logtb", std::cref(logtb)));
    
    ModelParamMap.insert(std::make_pair("mH1", std::cref(mH1)));
    ModelParamMap.insert(std::make_pair("mHp", std::cref(mHp1)));
    ModelParamMap.insert(std::make_pair("mH2", std::cref(mH21)));
    ModelParamMap.insert(std::make_pair("mH3", std::cref(mH31)));

    ModelParamMap.insert(std::make_pair("mH1sq", std::cref(mH1sq)));
    ModelParamMap.insert(std::make_pair("mHp2", std::cref(mHp2)));
    ModelParamMap.insert(std::make_pair("mH2sq", std::cref(mH2sq)));
    ModelParamMap.insert(std::make_pair("mH3sq", std::cref(mH3sq)));
    ModelParamMap.insert(std::make_pair("alpha1", std::cref(alpha1)));
    ModelParamMap.insert(std::make_pair("alpha2", std::cref(alpha2)));
    ModelParamMap.insert(std::make_pair("alpha3", std::cref(alpha3)));
    ModelParamMap.insert(std::make_pair("lambda1", std::cref(lambda1)));
    ModelParamMap.insert(std::make_pair("lambda2", std::cref(lambda2)));
    ModelParamMap.insert(std::make_pair("lambda3", std::cref(lambda3)));
    ModelParamMap.insert(std::make_pair("lambda4", std::cref(lambda4)));
    ModelParamMap.insert(std::make_pair("Relambda5", std::cref(Relambda5)));
    ModelParamMap.insert(std::make_pair("Imlambda5", std::cref(Imlambda5)));
    ModelParamMap.insert(std::make_pair("Relambda6", std::cref(Relambda6)));
    ModelParamMap.insert(std::make_pair("Relambda7", std::cref(Relambda7)));
    
    ModelParamMap.insert(std::make_pair("yu1R_GTHDM", std::cref(yu1R_GTHDM)));
    ModelParamMap.insert(std::make_pair("yd1R_GTHDM", std::cref(yd1R_GTHDM)));
    ModelParamMap.insert(std::make_pair("yl1R_GTHDM", std::cref(yl1R_GTHDM)));


 
    ModelParamMap.insert(std::make_pair("Nu_11r", std::cref(Nu_11r)));
    ModelParamMap.insert(std::make_pair("Nu_11i", std::cref(Nu_11i)));
    ModelParamMap.insert(std::make_pair("Nu_12r", std::cref(Nu_12r)));
    ModelParamMap.insert(std::make_pair("Nu_12i", std::cref(Nu_12i)));
    ModelParamMap.insert(std::make_pair("Nu_13r", std::cref(Nu_13r)));
    ModelParamMap.insert(std::make_pair("Nu_13i", std::cref(Nu_13i)));
    ModelParamMap.insert(std::make_pair("Nu_21r", std::cref(Nu_21r)));
    ModelParamMap.insert(std::make_pair("Nu_21i", std::cref(Nu_21i)));
    ModelParamMap.insert(std::make_pair("Nu_22r", std::cref(Nu_22r)));
    ModelParamMap.insert(std::make_pair("Nu_22i", std::cref(Nu_22i)));
    ModelParamMap.insert(std::make_pair("Nu_23r", std::cref(Nu_23r)));
    ModelParamMap.insert(std::make_pair("Nu_23i", std::cref(Nu_23i)));
    ModelParamMap.insert(std::make_pair("Nu_31r", std::cref(Nu_31r)));
    ModelParamMap.insert(std::make_pair("Nu_31i", std::cref(Nu_31i)));
    ModelParamMap.insert(std::make_pair("Nu_32r", std::cref(Nu_32r)));
    ModelParamMap.insert(std::make_pair("Nu_32i", std::cref(Nu_32i)));
    ModelParamMap.insert(std::make_pair("Nu_33r", std::cref(Nu_33r)));
    ModelParamMap.insert(std::make_pair("Nu_33i", std::cref(Nu_33i)));
    ModelParamMap.insert(std::make_pair("Nd_11r", std::cref(Nd_11r)));
    ModelParamMap.insert(std::make_pair("Nd_11i", std::cref(Nd_11i)));
    ModelParamMap.insert(std::make_pair("Nd_12r", std::cref(Nd_12r)));
    ModelParamMap.insert(std::make_pair("Nd_12i", std::cref(Nd_12i)));
    ModelParamMap.insert(std::make_pair("Nd_13r", std::cref(Nd_13r)));
    ModelParamMap.insert(std::make_pair("Nd_13i", std::cref(Nd_13i)));
    ModelParamMap.insert(std::make_pair("Nd_21r", std::cref(Nd_21r)));
    ModelParamMap.insert(std::make_pair("Nd_21i", std::cref(Nd_21i)));
    ModelParamMap.insert(std::make_pair("Nd_22r", std::cref(Nd_22r)));
    ModelParamMap.insert(std::make_pair("Nd_22i", std::cref(Nd_22i)));
    ModelParamMap.insert(std::make_pair("Nd_23r", std::cref(Nd_23r)));
    ModelParamMap.insert(std::make_pair("Nd_23i", std::cref(Nd_23i)));
    ModelParamMap.insert(std::make_pair("Nd_31r", std::cref(Nd_31r)));
    ModelParamMap.insert(std::make_pair("Nd_31i", std::cref(Nd_31i)));
    ModelParamMap.insert(std::make_pair("Nd_32r", std::cref(Nd_32r)));
    ModelParamMap.insert(std::make_pair("Nd_32i", std::cref(Nd_32i)));
    ModelParamMap.insert(std::make_pair("Nd_33r", std::cref(Nd_33r)));
    ModelParamMap.insert(std::make_pair("Nd_33i", std::cref(Nd_33i)));
    ModelParamMap.insert(std::make_pair("Nl_11r", std::cref(Nl_11r)));
    ModelParamMap.insert(std::make_pair("Nl_11i", std::cref(Nl_11i)));
    ModelParamMap.insert(std::make_pair("Nl_12r", std::cref(Nl_12r)));
    ModelParamMap.insert(std::make_pair("Nl_12i", std::cref(Nl_12i)));
    ModelParamMap.insert(std::make_pair("Nl_13r", std::cref(Nl_13r)));
    ModelParamMap.insert(std::make_pair("Nl_13i", std::cref(Nl_13i)));
    ModelParamMap.insert(std::make_pair("Nl_21r", std::cref(Nl_21r)));
    ModelParamMap.insert(std::make_pair("Nl_21i", std::cref(Nl_21i)));
    ModelParamMap.insert(std::make_pair("Nl_22r", std::cref(Nl_22r)));
    ModelParamMap.insert(std::make_pair("Nl_22i", std::cref(Nl_22i)));
    ModelParamMap.insert(std::make_pair("Nl_23r", std::cref(Nl_23r)));
    ModelParamMap.insert(std::make_pair("Nl_23i", std::cref(Nl_23i)));
    ModelParamMap.insert(std::make_pair("Nl_31r", std::cref(Nl_31r)));
    ModelParamMap.insert(std::make_pair("Nl_31i", std::cref(Nl_31i)));
    ModelParamMap.insert(std::make_pair("Nl_32r", std::cref(Nl_32r)));
    ModelParamMap.insert(std::make_pair("Nl_32i", std::cref(Nl_32i)));
    ModelParamMap.insert(std::make_pair("Nl_33r", std::cref(Nl_33r)));
    ModelParamMap.insert(std::make_pair("Nl_33i", std::cref(Nl_33i)));
    ModelParamMap.insert(std::make_pair("Q_GTHDM", std::cref(Q_GTHDM)));
    ModelParamMap.insert(std::make_pair("RpepsGTHDM", std::cref(RpepsGTHDM)));
    ModelParamMap.insert(std::make_pair("NLOuniscaleGTHDM", std::cref(NLOuniscaleGTHDM)));
    flag_use_sq_masses=true;
    flag_sigma=true;
    flag_SM_Higgs= true;


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
    setModelInitialized(NPbase::InitializeModel());
    setModelGeneralTHDM();
    return(IsModelInitialized());
}


bool GeneralTHDM::Update(const std::map<std::string, double>& DPars)
{
    if(!NPbase::Update(DPars)) return (false);

    
    if(flag_CPconservation){
        mu2=mHp2-lambda3/(2*v()*v());
        
        mH1=mHl;
                
        mH1sq=mH1*mH1;
        
        lambda1 = (mH1sq+mH2sq*tanalpha1*tanalpha1)/(v()*v()*(1+tanalpha1*tanalpha1));
    
        lambda4 = (mH1sq+mH3sq-2*mHp2+(mH2sq-mH1sq)/(1+tanalpha1*tanalpha1))/(v()*v());
    
        Relambda5 = ((mH2sq+mH1sq*tanalpha1*tanalpha1)/(1+tanalpha1*tanalpha1)-mH3sq)/(v()*v());
    
        Imlambda5 = 0.;
    
        Relambda6 = ((mH1sq-mH2sq)*tanalpha1)/(v()*v()*(1+tanalpha1*tanalpha1));
    
        Imlambda6 = 0.;
        
        Imlambda7 = 0.;
    }
    else{
        throw std::runtime_error("\033[1;31m The CP-Violating GeneralTHDM is still not implemented, please use the CP-conserving model \033[0m ");
    }
    

    return (true);
}



bool GeneralTHDM::PostUpdate()
{
    if(!NPbase::PostUpdate()) return (false);

    myGTHDMcache->updateCache();

    return (true);
}

void GeneralTHDM::setParameter(const std::string name, const double& value){
    
//    std::cout<<"\033[1;33m name = \033[0m "<< name << std::endl;
//    std::cout<<"\033[1;33m value = \033[0m "<< value << std::endl;

    
    if(name.compare("mH21") == 0){
        if(!flag_use_sq_masses){
            mH21 = value;
            mH2sq=mH21*mH21;
        }
        else{
            throw std::runtime_error(" If flag_use_sq_masses is true you should not include"
                    " the linear mass as a parameter in the configuration file. Please comment that line");

            //mH21=sqrt(mH2sq);
        }
    }
    else if(name.compare("mH31") == 0){
        if(!flag_use_sq_masses){
            mH31 = value;
            mH3sq = mH31*mH31;
        }
        else{
            throw std::runtime_error(" If flag_use_sq_masses is true you should not include"
                    " the linear mass as a parameter in the configuration file. Please comment that line");
            //mH31=sqrt(mH3sq);
        }
    }
    else if(name.compare("mHp1") == 0){
        if(!flag_use_sq_masses){
            mHp1 = value;
            mHp2 = mHp1*mHp1;
        }
        else{
            throw std::runtime_error(" If flag_use_sq_masses is true you should not include"
                    " the linear mass as a parameter in the configuration file. Please comment that line");
            //mHp1=sqrt(mHp2);
        }
    }

    else if(name.compare("mH2sq") == 0){
        if(flag_use_sq_masses){
            mH2sq = value;
            mH21=sqrt(mH2sq);
        }
        else{
            throw std::runtime_error(" If flag_use_sq_masses is false you should not include"
                    " the quadratic mass as a parameter in the configuration file. Please comment that line");
            //mH2sq=mH21*mH21;
        }
    }
    else if(name.compare("mH3sq") == 0){
        if(flag_use_sq_masses){
            mH3sq = value;
            mH31=sqrt(mH3sq);
        }
        else{
            throw std::runtime_error(" If flag_use_sq_masses is false you should not include"
                    " the quadratic mass as a parameter in the configuration file. Please comment that line");
            //mH3sq=mH31*mH31;
        }
    }
    else if(name.compare("mHp2") == 0){
        if(flag_use_sq_masses){
            mHp2 = value;
            mHp1=sqrt(mHp2);
        }
        else{
            throw std::runtime_error(" If flag_use_sq_masses is false you should not include"
                    " the quadratic mass as a parameter in the configuration file. Please comment that line");
            //mHp2=mHp1*mHp1;
        }
    }

    else if(name.compare("alpha1") == 0) {
        alpha1 = value;
        cosalpha1 = cos(alpha1);
        sinalpha1 = sin(alpha1);
        if(cosalpha1 != 0){
            tanalpha1 = sinalpha1/cosalpha1;
        }
        else{
            std::cout<<"\033[1;33m Warning!!! The alpha1=pi/2 and it's tangent is infinity, reduce the range of alpha1 \033[0m "<<std::endl;
        }

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
    else if(name.compare("lambda2") == 0)
        lambda2 = value;
    else if(name.compare("lambda3") == 0)
        lambda3 = value;
    else if(name.compare("Relambda7") == 0) 
        Relambda7 = value;
    else if(name.compare("Nu_11r") == 0 && flag_sigma)
        Nu_11r = value;
   else if(name.compare("yu1R_GTHDM") == 0 && !flag_sigma){
       if(alpha1==0 && flag_SM_Higgs )
           yu1R_GTHDM=1;
       else
         yu1R_GTHDM = value;
   }
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
    else if(name.compare("Nd_11r") == 0 && flag_sigma)
        Nd_11r = value;
    else if(name.compare("yd1R_GTHDM") == 0 && !flag_sigma)
        if(alpha1==0)
           yd1R_GTHDM=1;
       else
         yd1R_GTHDM = value;
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
    else if(name.compare("Nl_11r") == 0 && flag_sigma)
        Nl_11r = value;
    else if(name.compare("yl1R_GTHDM") == 0 && !flag_sigma)
        if(alpha1==0)
           yl1R_GTHDM=1;
       else
         yl1R_GTHDM = value;
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
    else if(name.compare("RpepsGTHDM") == 0)
        RpepsGTHDM = value;
    else if(name.compare("NLOuniscaleGTHDM") == 0)
        NLOuniscaleGTHDM = value;
    else
        NPbase::setParameter(name,value);
    
}

bool GeneralTHDM::CheckParameters(const std::map<std::string, double>& DPars) {
    for (int i = 0; i < NGeneralTHDMvars; i++) {
        if (DPars.find(GeneralTHDMvars[i]) == DPars.end()) {
            std::cout << "ERROR: missing mandatory GeneralTHDM parameter " << GeneralTHDMvars[i] << std::endl;
            raiseMissingModelParameterCount();
            addMissingModelParameter(GeneralTHDMvars[i]);
        }
    }
    return(NPbase::CheckParameters(DPars));
}

///////////////////////////////////////////////////////////////////////////
// Flags

bool GeneralTHDM::setFlagStr(const std::string name, const std::string value)
{
    bool res = false;
    if(name.compare("RGEorder") == 0)
    {
        if (checkRGEorder(value))
        {
            flag_RGEorder = value;
            res = true;
        }
        else
        {
            throw std::runtime_error("GeneralTHDM::setFlagStr(): Invalid flag "
                + name + "=" + value);
        }
    }
    else
    {
        res = NPbase::setFlagStr(name,value);
    }

    return(res);
}

bool GeneralTHDM::setFlag(const std::string name, const bool value)
{
    bool res = false;
    if(name.compare("use_sq_masses") == 0) {
        flag_use_sq_masses = value;
        res = true;
        if (!flag_use_sq_masses) {
           GeneralTHDMvars[std::distance(GeneralTHDMvars,std::find(GeneralTHDMvars,GeneralTHDMvars+NGeneralTHDMvars,"mH2sq"))] = "mH21";
           GeneralTHDMvars[std::distance(GeneralTHDMvars,std::find(GeneralTHDMvars,GeneralTHDMvars+NGeneralTHDMvars,"mH3sq"))] = "mH31";
           GeneralTHDMvars[std::distance(GeneralTHDMvars,std::find(GeneralTHDMvars,GeneralTHDMvars+NGeneralTHDMvars,"mHp2"))] = "mHp1";
           
           ModelParamMap.insert(std::make_pair("mH21", std::cref(mH21)));
           ModelParamMap.insert(std::make_pair("mH31", std::cref(mH31)));
           ModelParamMap.insert(std::make_pair("mHp1", std::cref(mHp1)));
           
        }
    }
    else if(name.compare("use_sigma") == 0) {
    std::cout<<"use_sigma = "<< value << std::endl;
        flag_sigma = value;
        res = true;
        if (!flag_sigma && flag_CPconservation) {
           GeneralTHDMvars[std::distance(GeneralTHDMvars,std::find(GeneralTHDMvars,GeneralTHDMvars+NGeneralTHDMvars,"Nu_11r"))] = "yu1R_GTHDM";
           GeneralTHDMvars[std::distance(GeneralTHDMvars,std::find(GeneralTHDMvars,GeneralTHDMvars+NGeneralTHDMvars,"Nd_11r"))] = "yd1R_GTHDM";
           GeneralTHDMvars[std::distance(GeneralTHDMvars,std::find(GeneralTHDMvars,GeneralTHDMvars+NGeneralTHDMvars,"Nl_11r"))] = "yl1R_GTHDM";
           
           ModelParamMap.insert(std::make_pair("yu1R_GTHDM", std::cref(yu1R_GTHDM)));
           ModelParamMap.insert(std::make_pair("yd1R_GTHDM", std::cref(yd1R_GTHDM)));
           ModelParamMap.insert(std::make_pair("yl1R_GTHDM", std::cref(yl1R_GTHDM)));
        }
        if(!flag_sigma && !flag_CPconservation){
           throw std::runtime_error("sigma flag different from true only in the CP-conserving");
        }
    }
    else if(name.compare("ATHDMflag") == 0) {
    std::cout<<"ATHDMflag = "<< value<<std::endl;
        flag_ATHDM = value;
        res = true;
    }
    else if(name.compare("CPconservation") == 0) {
    std::cout<<"CPconservation = "<< value<<std::endl;
        flag_CPconservation = value;
        res = true;  
    }
    else if(name.compare("SMHiggs") == 0) {
    std::cout<<"SMHiggs = "<< value<<std::endl;
        flag_SM_Higgs = value;
        res = true;  
    }
    else
        res = NPbase::setFlag(name,value);

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




double GeneralTHDM::F(const double m02, const double m12) const {
    double F;

    if(m02 == 0. && m12 != 0.) {
        F=0.5 * m12;
    } else if(m02 != 0. && m12 == 0.){
        F=0.5 * m02;
    } else if((m02 == 0. && m12 == 0.) || (fabs(m02-m12) < LEPS)){
        F=0.;
    } else if (m02 != 0 && m12 != 0){
        F=0.5 * (m02 + m12) - (m02 * m12) / (m02 - m12) * log(m02 / m12);
    } else
        throw std::runtime_error("Error in GeneralTHDM::F()");
    return (F);
}



double GeneralTHDM::GTHDMDeltaS() const
{

    
    //MAYBE WE CAN USE DIRECTLY THE PARAMETERS mH1sq, mH2sq, mH3sq, mHp2. IT'S A NONESENSE TO DEFINE THIS AGAIN...
    double mH1_2 = mH1sq;
    double mH2_2 = mH2sq;
    double mH3_2 = mH3sq;
    double mHp2= getmHp2();
    
    double mHref_2 = getmH1sq();  
    double  R11 = 0.0;
    double  R21 = 0.0;
    double  R31 = cosalpha2*sinalpha2*cosalpha3 + sinalpha1*sinalpha3;

    //THIS NEEDS TO BE CHECKED, MAYBE BETTER TO DEFINE THIS MATRIX ELEMENTS AS PARAMETERS
    if(getSMHiggs()){
           R11 = cosalpha1*cosalpha2;;
           R21 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
    }
    else{
           R21 = cosalpha1*cosalpha2;
           R11 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
    }
    
    double  R11_2 = R11*R11;
    double  R21_2 = R21*R21;
    double  R31_2 = R31*R31;

    //std::cout<<"\033[1;33m R11 =  \033[0m "<<  R11 <<std::endl;
    //std::cout<<"\033[1;33m R21 =  \033[0m "<<  R21 <<std::endl;
    //std::cout<<"\033[1;33m R31 =  \033[0m "<<  R31 <<std::endl;
    
    
    double MZ= getMz();
    double MZ2 = MZ*MZ;
    
    gslpp::complex B00prime_MZ2_MZ2_mH_2_mH_3;
    gslpp::complex B00prime_MZ2_MZ2_mHp2_mHp2;
    gslpp::complex B00prime_MZ2_MZ2_mH_1_mH_3;
    gslpp::complex B00prime_MZ2_MZ2_mH_1_mH_2;

    
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mH2_2;
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mH1_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mH2_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mH1_2;
    
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mH3_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mH3_2;
    
    gslpp::complex B00prime_MZ2_MZ2_MZ2_mHref_2;
    gslpp::complex B0prime_MZ2_MZ2_MZ2_mHref_2;



    B00prime_MZ2_MZ2_mH_2_mH_3 = myGTHDMcache->B00_MZ2_MZ2_mHh2_mA2(MZ2,mH2_2,mH3_2) - myGTHDMcache->B00_MZ2_0_mHh2_mA2(MZ2,mH2_2,mH3_2);
    B00prime_MZ2_MZ2_mHp2_mHp2 = myGTHDMcache->B00_MZ2_MZ2_mHp2_mHp2(MZ2,mHp2) - myGTHDMcache->B00_MZ2_0_mHp2_mHp2(MZ2,mHp2);
    B00prime_MZ2_MZ2_mH_1_mH_3 = myGTHDMcache->B00_MZ2_MZ2_mHl2_mA2(MZ2,mH1_2,mH3_2) - myGTHDMcache->B00_MZ2_0_mHl2_mA2(MZ2,mH1_2,mH3_2);
    B00prime_MZ2_MZ2_mH_1_mH_2 = myGTHDMcache->B00_MZ2_MZ2_mHl2_mA2(MZ2,mH1_2,mH2_2) - myGTHDMcache->B00_MZ2_0_mHl2_mA2(MZ2,mH1_2,mH2_2);

    
    B00prime_MZ2_MZ2_MZ2_mH2_2 = myGTHDMcache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mH2_2) - myGTHDMcache->B00_MZ2_0_MZ2_mHh2(MZ2,mH2_2);
    B00prime_MZ2_MZ2_MZ2_mH1_2 = myGTHDMcache->B00_MZ2_MZ2_MZ2_mHl2(MZ2,mH1_2) - myGTHDMcache->B00_MZ2_0_MZ2_mHl2(MZ2,mH1_2);
    B0prime_MZ2_MZ2_MZ2_mH2_2 = myGTHDMcache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mH2_2) - myGTHDMcache->B0_MZ2_0_MZ2_mHh2(MZ2,mH2_2);
    B0prime_MZ2_MZ2_MZ2_mH1_2 = myGTHDMcache->B0_MZ2_MZ2_MZ2_mHl2(MZ2,mH1_2) - myGTHDMcache->B0_MZ2_0_MZ2_mHl2(MZ2,mH1_2);
    B00prime_MZ2_MZ2_MZ2_mH3_2 = myGTHDMcache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mH3_2) - myGTHDMcache->B00_MZ2_0_MZ2_mHh2(MZ2,mH3_2);
    B0prime_MZ2_MZ2_MZ2_mH3_2 = myGTHDMcache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mH3_2) - myGTHDMcache->B0_MZ2_0_MZ2_mHh2(MZ2,mH3_2);

    B00prime_MZ2_MZ2_MZ2_mHref_2 = myGTHDMcache->B00_MZ2_MZ2_MZ2_mHh2(MZ2,mHref_2) - myGTHDMcache->B00_MZ2_0_MZ2_mHh2(MZ2,mHref_2);
    B0prime_MZ2_MZ2_MZ2_mHref_2 = myGTHDMcache->B0_MZ2_MZ2_MZ2_mHh2(MZ2,mHref_2) - myGTHDMcache->B0_MZ2_0_MZ2_mHh2(MZ2,mHref_2);



    
    return 1/MZ2/M_PI*( R11_2*(B00prime_MZ2_MZ2_MZ2_mH1_2.real() - MZ2*B0prime_MZ2_MZ2_MZ2_mH1_2.real())+
            R21_2*(B00prime_MZ2_MZ2_MZ2_mH2_2.real() - MZ2*B0prime_MZ2_MZ2_MZ2_mH2_2.real())+
            R31_2*(B00prime_MZ2_MZ2_MZ2_mH3_2.real() - MZ2*B0prime_MZ2_MZ2_MZ2_mH3_2.real())+
             R11_2*B00prime_MZ2_MZ2_mH_2_mH_3.real() + R21_2*B00prime_MZ2_MZ2_mH_1_mH_3.real()
             +R31_2*B00prime_MZ2_MZ2_mH_1_mH_2.real() - B00prime_MZ2_MZ2_mHp2_mHp2.real()
              - B00prime_MZ2_MZ2_MZ2_mHref_2.real() + MZ2*B0prime_MZ2_MZ2_MZ2_mHref_2.real());
     
}



double GeneralTHDM::GTHDMDeltaT() const
{

    gslpp::complex I = gslpp::complex::i();
    double mH1_2 = mH1sq;
    double mH2_2 = mH2sq;
    double mH3_2 = mH3sq;
    double mHp2 = getmHp2();
    
    double mHref_2 = mH1sq;
        
        
        
    double  R11 = 0.0;
    double  R12 = 0.0;
    double  R13 = 0.0;
    double  R21 = 0.0;
    double  R22 = 0.0;
    double  R23 = 0.0;
    double R31 = cosalpha2*sinalpha2*cosalpha3 + sinalpha1*sinalpha3;
    double R32 = sinalpha1*sinalpha2*cosalpha3 - cosalpha1*sinalpha3;
    double R33 = cosalpha2*cosalpha3;
    
    
    if(getSMHiggs()){
        R11 = cosalpha1*cosalpha2;
        R12 = sinalpha1*cosalpha2;
        R13 = -sinalpha2;
        R21 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
        R22 = sinalpha1*sinalpha2*sinalpha3 + cosalpha1*cosalpha3;
        R23 = cosalpha2*sinalpha3;
    }
    else{
        
        R21 = cosalpha1*cosalpha2; //R21 is R11 in this case
        R22 = sinalpha1*cosalpha2; //R22 is R12 in this case
        R23 = -sinalpha2;          //R23 is R13 in this case
        R11 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3; //R11 is R21 in this case
        R12 = sinalpha1*sinalpha2*sinalpha3 + cosalpha1*cosalpha3; //R12 is R22 in this case
        R13 = cosalpha2*sinalpha3; //R13 is R23 in this case
    }
    
        double  R11_2 = R11*R11;
        double  R21_2 = R21*R21;
        double  R31_2 = R31*R31;
     
    double MZ=getMz();
    double MZ2 = MZ*MZ;
    
    double MW=Mw_tree();
    double MW2 = MW*MW;
    
    double s_W2 = 1.0-c02();

    gslpp::complex B0_MZ2_0_MZ2_mH1_2;
    gslpp::complex B0_MZ2_0_MZ2_mH2_2;
    gslpp::complex B0_MZ2_0_MZ2_mH3_2;
    gslpp::complex B0_MZ2_0_MZ2_mHref_2;


    gslpp::complex B0_MZ2_0_MW2_mH1_2;
    gslpp::complex B0_MZ2_0_MW2_mH2_2;
    gslpp::complex B0_MZ2_0_MW2_mH3_2;
    gslpp::complex B0_MZ2_0_MW2_mHref_2;


    B0_MZ2_0_MZ2_mH1_2 = myGTHDMcache->B0_MZ2_0_MZ2_mHh2(MZ2,mH1_2);
    B0_MZ2_0_MZ2_mH2_2 = myGTHDMcache->B0_MZ2_0_MZ2_mHl2(MZ2,mH2_2);
    B0_MZ2_0_MZ2_mH3_2 = myGTHDMcache->B0_MZ2_0_MZ2_mHl2(MZ2,mH3_2);
    B0_MZ2_0_MZ2_mHref_2 = myGTHDMcache->B0_MZ2_0_MZ2_mHl2(MZ2,mHref_2);

    
    B0_MZ2_0_MW2_mH1_2 = myGTHDMcache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mH1_2);
    B0_MZ2_0_MW2_mH2_2 = myGTHDMcache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mH2_2);
    B0_MZ2_0_MW2_mH3_2 = myGTHDMcache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mH3_2);
    B0_MZ2_0_MW2_mHref_2 = myGTHDMcache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHref_2);

      
    return 1. / 16. / M_PI / MW2 / s_W2*((R12 + I*R13).abs2()*F(mHp2, mH1_2)+(R22 + I*R23).abs2()*F(mHp2, mH2_2)+
             (R32 + I*R33).abs2()*F(mHp2, mH3_2)- R11_2*F(mH2_2, mH3_2) - R21_2*F(mH1_2, mH3_2) - R31_2*F(mH1_2, mH2_2)
            + R11_2*(F(MW2, mH1_2) - F(MZ2, mH1_2) - 4.0*MW2*B0_MZ2_0_MW2_mH1_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mH1_2.real())
            + R21_2*(F(MW2, mH2_2) - F(MZ2, mH2_2) - 4.0*MW2*B0_MZ2_0_MW2_mH2_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mH2_2.real())
            + R31_2*(F(MW2, mH3_2) - F(MZ2, mH3_2) - 4.0*MW2*B0_MZ2_0_MW2_mH3_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mH3_2.real())
            -(F(MW2, mHref_2) - F(MZ2, mHref_2) - 4.0*MW2*B0_MZ2_0_MW2_mHref_2.real() + 4.0*MZ2*B0_MZ2_0_MZ2_mHref_2.real()));

}



double GeneralTHDM::GTHDMDeltaU() const
{
    
    gslpp::complex I = gslpp::complex::i();
    double mH1_2 = mH1sq;
    double mH2_2 = mH2sq;
    double mH3_2 = mH3sq;
    double mHp2 = getmHp2();
    
    double mHref_2 = mH1sq;
        
        
    
          
    double  R11 = 0.0;
    double  R12 = 0.0;
    double  R13 = 0.0;
    double  R21 = 0.0;
    double  R22 = 0.0;
    double  R23 = 0.0;
    double  R31 = cosalpha2*sinalpha2*cosalpha3 + sinalpha1*sinalpha3;
    double  R32 = sinalpha1*sinalpha2*cosalpha3 - cosalpha1*sinalpha3;
    double  R33 = cosalpha2*cosalpha3;

    
    
    if(getSMHiggs()){
        R11 = cosalpha1*cosalpha2;
        R12 = sinalpha1*cosalpha2;
        R13 = -sinalpha2;
        R21 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3;
        R22 = sinalpha1*sinalpha2*sinalpha3 + cosalpha1*cosalpha3;
        R23 = cosalpha2*sinalpha3;
    }
    else{
        
        R21 = cosalpha1*cosalpha2; //R21 is R11 in this case
        R22 = sinalpha1*cosalpha2; //R22 is R12 in this case
        R23 = -sinalpha2;          //R23 is R13 in this case
        R11 = cosalpha1*sinalpha2*sinalpha3 - sinalpha1*cosalpha3; //R11 is R21 in this case
        R12 = sinalpha1*sinalpha2*sinalpha3 + cosalpha1*cosalpha3; //R12 is R22 in this case
        R13 = cosalpha2*sinalpha3; //R13 is R23 in this case
    }
    
    
    double  R11_2 = R11*R11;
    double  R21_2 = R21*R21;
    double  R31_2 = R31*R31;

    double MZ=getMz();
    double MZ2 = MZ*MZ;
    
    double MW=Mw_tree();
    double MW2 = MW*MW;


    gslpp::complex B00prime_MZ2_MW2_mH1_2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mH2_2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mH3_2_mHp2;
    gslpp::complex B00prime_MZ2_MW2_mHref_2_mHp2;


    
    gslpp::complex B00prime_MZ2_MW2_mHp2_mHp2;
    
    
    gslpp::complex B00prime_MZ2_MW2_MW2_mH1_2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mH2_2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mH3_2;
    gslpp::complex B00prime_MZ2_MW2_MW2_mHref_2;

    
    gslpp::complex B0prime_MZ2_MW2_MW2_mH1_2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mH2_2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mH3_2;
    gslpp::complex B0prime_MZ2_MW2_MW2_mHref_2;



    B00prime_MZ2_MW2_mH1_2_mHp2 = myGTHDMcache->B00_MZ2_MW2_mHl2_mHp2(MZ2,MW2,mH1_2,mHp2) - myGTHDMcache->B00_MZ2_0_mHl2_mHp2(MZ2,mH1_2,mHp2);
    B00prime_MZ2_MW2_mH2_2_mHp2 = myGTHDMcache->B00_MZ2_MW2_mHh2_mHp2(MZ2,MW2,mH2_2,mHp2) - myGTHDMcache->B00_MZ2_0_mHh2_mHp2(MZ2,mH2_2,mHp2);
    B00prime_MZ2_MW2_mH3_2_mHp2 = myGTHDMcache->B00_MZ2_MW2_mHh2_mHp2(MZ2,MW2,mH3_2,mHp2) - myGTHDMcache->B00_MZ2_0_mHh2_mHp2(MZ2,mH3_2,mHp2);
    B00prime_MZ2_MW2_mHref_2_mHp2 = myGTHDMcache->B00_MZ2_MW2_mHl2_mHp2(MZ2,MW2,mHref_2,mHp2) - myGTHDMcache->B00_MZ2_0_mHl2_mHp2(MZ2,mHref_2,mHp2);

    B00prime_MZ2_MW2_mHp2_mHp2 = myGTHDMcache->B00_MZ2_MW2_mHp2_mHp2(MZ2,MW2,mHp2) - myGTHDMcache->B00_MZ2_0_mHp2_mHp2(MZ2,mHp2);

    B00prime_MZ2_MW2_MW2_mH1_2 = myGTHDMcache->B00_MZ2_MW2_MW2_mHl2(MZ2,MW2,mH1_2) - myGTHDMcache->B00_MZ2_0_MW2_mHl2(MZ2,MW2,mH1_2);
    B00prime_MZ2_MW2_MW2_mH2_2 = myGTHDMcache->B00_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH2_2) - myGTHDMcache->B00_MZ2_0_MW2_mHh2(MZ2,MW2,mH2_2);
    B00prime_MZ2_MW2_MW2_mH3_2 = myGTHDMcache->B00_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH3_2) - myGTHDMcache->B00_MZ2_0_MW2_mHh2(MZ2,MW2,mH3_2);
    B00prime_MZ2_MW2_MW2_mHref_2 = myGTHDMcache->B00_MZ2_MW2_MW2_mHl2(MZ2,MW2,mHref_2) - myGTHDMcache->B00_MZ2_0_MW2_mHl2(MZ2,MW2,mHref_2);

    B0prime_MZ2_MW2_MW2_mH1_2 = myGTHDMcache->B0_MZ2_MW2_MW2_mHl2(MZ2,MW2,mH1_2) - myGTHDMcache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mH1_2);
    B0prime_MZ2_MW2_MW2_mH2_2 = myGTHDMcache->B0_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH2_2) - myGTHDMcache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mH2_2);
    B0prime_MZ2_MW2_MW2_mH3_2 = myGTHDMcache->B0_MZ2_MW2_MW2_mHh2(MZ2,MW2,mH3_2) - myGTHDMcache->B0_MZ2_0_MW2_mHh2(MZ2,MW2,mH3_2);
    B0prime_MZ2_MW2_MW2_mHref_2 = myGTHDMcache->B0_MZ2_MW2_MW2_mHl2(MZ2,MW2,mHref_2) - myGTHDMcache->B0_MZ2_0_MW2_mHl2(MZ2,MW2,mHref_2);

    
    return - GTHDMDeltaS() + 1. / M_PI / MZ2 * (- R11_2*MW2*B0prime_MZ2_MW2_MW2_mH1_2.real() - R21_2*MW2*B0prime_MZ2_MW2_MW2_mH2_2.real() - R31_2*MW2*B0prime_MZ2_MW2_MW2_mH3_2.real() 
    + MW2*B0prime_MZ2_MW2_MW2_mHref_2.real() - B00prime_MZ2_MW2_MW2_mHref_2.real()
    + R11_2* B00prime_MZ2_MW2_MW2_mH1_2.real() +R21_2* B00prime_MZ2_MW2_MW2_mH2_2.real() +R31_2* B00prime_MZ2_MW2_MW2_mH3_2.real()
    + (R12+I*R13).abs2()*B00prime_MZ2_MW2_mH1_2_mHp2.real() + (R22+I*R23).abs2()*B00prime_MZ2_MW2_mH2_2_mHp2.real()
    + (R32+I*R33).abs2()*B00prime_MZ2_MW2_mH3_2_mHp2.real() - 2. * B00prime_MZ2_MW2_mHp2_mHp2.real());

    
    
}





//Here we had before the at tree level but it would be better to include the full expression (including the oblique corrections).
double GeneralTHDM::Mw() const{


    //std::cout<<"\033[1;33m GTHDMDeltaS() =  \033[0m "<<  GTHDMDeltaS() <<std::endl;
    //std::cout<<"\033[1;32m GTHDMDeltaT() =  \033[0m "<<  GTHDMDeltaT() <<std::endl;
    //std::cout<<"\033[1;31m GTHDMDeltaU() =  \033[0m "<<  GTHDMDeltaU() <<std::endl;
    
    //std::cout<<"\033[1;31m trueSM.Mw() =  \033[0m "<<  trueSM.Mw() <<std::endl;
    
    double myMw = trueSM.Mw();

    double alpha = trueSM.alphaMz();
    double c2 = trueSM.cW2();
    double s2 = trueSM.sW2();

    myMw *= 1.0 - alpha / 4.0 / (c2 - s2)
            *(GTHDMDeltaS() - 2.0 * c2 * GTHDMDeltaT() - (c2 - s2) * GTHDMDeltaU() / 2.0 / s2)
            - s2 / 2.0 / (c2 - s2) * DeltaGF();
    
    
    //std::cout<<"\033[1;31m myMw =  \033[0m "<<  myMw <<std::endl;
     
    return (myMw);
//    return (trueSM.Mw());
//    return (Mw_tree());
}
