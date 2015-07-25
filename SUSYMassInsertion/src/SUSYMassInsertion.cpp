/* 
 * Copyright (C) 2012 HEPfit Collaboration
 * All rights reserved.
 *
 * For the licensing terms see doc/COPYING.
 */

#include "SUSYMassInsertion.h"

const std::string SUSYMassInsertion::SusyMIvars[NSusyMIvars] = {"Msq", "m3", "MuM",
      
        "rDULL11", "rDULL12", "rDULL13", "rDULL21", "rDULL22", "rDULL23", 
        "rDULL31", "rDULL32", "rDULL33", "iDULL11", "iDULL12", "iDULL13", 
        "iDULL21", "iDULL22", "iDULL23", "iDULL31", "iDULL32", "iDULL33",
    
        "rDURL11", "rDURL12", "rDURL13", "rDURL21", "rDURL22", "rDURL23", 
        "rDURL31", "rDURL32", "rDURL33", "iDURL11", "iDURL12", "iDURL13", 
        "iDURL21", "iDURL22", "iDURL23", "iDURL31", "iDURL32", "iDURL33",
    
        "rDULR11", "rDULR12", "rDULR13", "rDULR21", "rDULR22", "rDULR23", 
        "rDULR31", "rDULR32", "rDULR33", "iDULR11", "iDULR12", "iDULR13", 
        "iDULR21", "iDULR22", "iDULR23", "iDULR31", "iDULR32", "iDULR33",
    
        "rDURR11", "rDURR12", "rDURR13", "rDURR21", "rDURR22", "rDURR23", 
        "rDURR31", "rDURR32", "rDURR33", "iDURR11", "iDURR12", "iDURR13", 
        "iDURR21", "iDURR22", "iDURR23", "iDURR31", "iDURR32", "iDURR33",
    
        "rDDLL11", "rDDLL12", "rDDLL13", "rDDLL21", "rDDLL22", "rDDLL23", 
        "rDDLL31", "rDDLL32", "rDDLL33", "iDDLL11", "iDDLL12", "iDDLL13", 
        "iDDLL21", "iDDLL22", "iDDLL23", "iDDLL31", "iDDLL32", "iDDLL33",
    
        "rDDRL11", "rDDRL12", "rDDRL13", "rDDRL21", "rDDRL22", "rDDRL23", 
        "rDDRL31", "rDDRL32", "rDDRL33", "iDDRL11", "iDDRL12", "iDDRL13", 
        "iDDRL21", "iDDRL22", "iDDRL23", "iDDRL31", "iDDRL32", "iDDRL33",
    
        "rDDLR11", "rDDLR12", "rDDLR13", "rDDLR21", "rDDLR22", "rDDLR23", 
        "rDDLR31", "rDDLR32", "rDDLR33", "iDDLR11", "iDDLR12", "iDDLR13", 
        "iDDLR21", "iDDLR22", "iDDLR23", "iDDLR31", "iDDLR32", "iDDLR33",
    
        "rDDRR11", "rDDRR12", "rDDRR13", "rDDRR21", "rDDRR22", "rDDRR23", 
        "rDDRR31", "rDDRR32", "rDDRR33", "iDDRR11", "iDDRR12", "iDDRR13", 
        "iDDRR21", "iDDRR22", "iDDRR23", "iDDRR31", "iDDRR32", "iDDRR33"
    
};

SUSYMassInsertion::SUSYMassInsertion()
:  StandardModel(),
        Du_LL(3,0.), Du_LR(3,0.), Du_RL(3,0.), Du_RR(3,0.),
        Dd_LL(3,0.), Dd_LR(3,0.), Dd_RL(3,0.), Dd_RR(3,0.){
    
        rDULL11=0.; rDULL12=0.; rDULL13=0.; 
        rDULL21=0.; rDULL22=0.; rDULL23=0.; 
        rDULL31=0.; rDULL32=0.; rDULL33=0.;
        iDULL11=0.; iDULL12=0.; iDULL13=0.; 
        iDULL21=0.; iDULL22=0.; iDULL23=0.; 
        iDULL31=0.; iDULL32=0.; iDULL33=0.;
    
        rDURL11=0.; rDURL12=0.; rDURL13=0.; 
        rDURL21=0.; rDURL22=0.; rDURL23=0.; 
        rDURL31=0.; rDURL32=0.; rDURL33=0.;
        iDURL11=0.; iDURL12=0.; iDURL13=0.; 
        iDURL21=0.; iDURL22=0.; iDURL23=0.; 
        iDURL31=0.; iDURL32=0.; iDURL33=0.;
    
        rDULR11=0.; rDULR12=0.; rDULR13=0.; 
        rDULR21=0.; rDULR22=0.; rDULR23=0.; 
        rDULR31=0.; rDULR32=0.; rDULR33=0.;
        iDULR11=0.; iDULR12=0.; iDULR13=0.; 
        iDULR21=0.; iDULR22=0.; iDULR23=0.; 
        iDULR31=0.; iDULR32=0.; iDULR33=0.;
    
        rDURR11=0.; rDURR12=0.; rDURR13=0.; 
        rDURR21=0.; rDURR22=0.; rDURR23=0.; 
        rDURR31=0.; rDURR32=0.; rDURR33=0.;
        iDURR11=0.; iDURR12=0.; iDURR13=0.; 
        iDURR21=0.; iDURR22=0.; iDURR23=0.; 
        iDURR31=0.; iDURR32=0.; iDURR33=0.;
    
        rDDLL11=0.; rDDLL12=0.; rDDLL13=0.; 
        rDDLL21=0.; rDDLL22=0.; rDDLL23=0.; 
        rDDLL31=0.; rDDLL32=0.; rDDLL33=0.;
        iDDLL11=0.; iDDLL12=0.; iDDLL13=0.; 
        iDDLL21=0.; iDDLL22=0.; iDDLL23=0.; 
        iDDLL31=0.; iDDLL32=0.; iDDLL33=0.;
    
        rDDRL11=0.; rDDRL12=0.; rDDRL13=0.; 
        rDDRL21=0.; rDDRL22=0.; rDDRL23=0.; 
        rDDRL31=0.; rDDRL32=0.; rDDRL33=0.;
        iDDRL11=0.; iDDRL12=0.; iDDRL13=0.; 
        iDDRL21=0.; iDDRL22=0.; iDDRL23=0.; 
        iDDRL31=0.; iDDRL32=0.; iDDRL33=0.;
    
        rDDLR11=0.; rDDLR12=0.; rDDLR13=0.; 
        rDDLR21=0.; rDDLR22=0.; rDDLR23=0.; 
        rDDLR31=0.; rDDLR32=0.; rDDLR33=0.;
        iDDLR11=0.; iDDLR12=0.; iDDLR13=0.; 
        iDDLR21=0.; iDDLR22=0.; iDDLR23=0.; 
        iDDLR31=0.; iDDLR32=0.; iDDLR33=0.;
    
        rDDRR11=0.; rDDRR12=0.; rDDRR13=0.; 
        rDDRR21=0.; rDDRR22=0.; rDDRR23=0.; 
        rDDRR31=0.; rDDRR32=0.; rDDRR33=0.;
        iDDRR11=0.; iDDRR12=0.; iDDRR13=0.; 
        iDDRR21=0.; iDDRR22=0.; iDDRR23=0.; 
        iDDRR31=0.; iDDRR32=0.; iDDRR33=0.;
    
}

SUSYMassInsertion::~SUSYMassInsertion(){    
}

///////////////////////////////////////////////////////////////////////////
// Initialization

bool SUSYMassInsertion::InitializeModel()
{
    mySUSYMIA = new SUSYMassInsertionMatching(*this);
    setModelInitialized(true);
    return(true);
}

///////////////////////////////////////////////////////////////////////////
// Parameters 

bool SUSYMassInsertion::Init(const std::map<std::string, double>& DPars) {
    
    if ((DPars.find("Msq") == DPars.end()) || (DPars.find("m3") == DPars.end()) ||
        (DPars.find("MuM") == DPars.end())) {
        std::cout << "missing mandatory Susy MIA parameter " << std::endl;
        std::cout << "the gluino mass, the mean value of the squark masses "
                     "and the matching scale must be initialized " << std::endl;
        return false;
    }

    return(StandardModel::Init(DPars));
}


bool SUSYMassInsertion::PreUpdate()
{    
    if(!StandardModel::PreUpdate()) return (false);
    return (true);
}

bool SUSYMassInsertion::Update(const std::map<std::string, double>& DPars)
{    
    if(!PreUpdate()) return (false);
    
    UpdateError = false;
    
    for (std::map<std::string, double>::const_iterator it = DPars.begin(); it != DPars.end(); it++)
        setParameter(it->first, it->second);
    
    if (UpdateError) return (false);
    
    if(!PostUpdate()) return (false);
    
    return (true);
}

bool SUSYMassInsertion::PostUpdate()
{
    if (!StandardModel::PostUpdate()) return (false);
    mySUSYMIA->StandardModelMatching::updateSMParameters();
    return (true);
}

void SUSYMassInsertion::setParameter(const std::string name, const double& value) {
    if (name.compare("Msq") == 0)
        Msq = value;
    else if (name.compare("m3") == 0)
        m3 = value;
    else if (name.compare("MuM") == 0)
        MuM = value;
    
    //FILLS THE DELTAS, UP SECTOR
    else if (name.compare("rDULL11") == 0){
        rDULL11 = value;
        Du_LL.assignre(0,0,rDULL11);}
    else if (name.compare("rDULL12") == 0){
        rDULL12 = value;
        Du_LL.assignre(0,1,rDULL12);}
    else if (name.compare("rDULL13") == 0){
        rDULL13 = value;
        Du_LL.assignre(0,2,rDULL13);}
    else if (name.compare("rDULL21") == 0){
        rDULL21 = value;
        Du_LL.assignre(1,0,rDULL21);}
    else if (name.compare("rDULL22") == 0){
        rDULL22 = value;
        Du_LL.assignre(1,1,rDULL22);}
    else if (name.compare("rDULL23") == 0){
        rDULL23 = value;
        Du_LL.assignre(1,2,rDULL23);}
    else if (name.compare("rDULL31") == 0){
        rDULL31 = value;
        Du_LL.assignre(2,0,rDULL31);}
    else if (name.compare("rDULL32") == 0){
        rDULL32 = value;
        Du_LL.assignre(2,1,rDULL32);}
    else if (name.compare("rDULL33") == 0){
        rDULL33 = value;
        Du_LL.assignre(2,2,rDULL33);}
    
    else if (name.compare("iDULL11") == 0){
        iDULL11 = value;
        Du_LL.assignim(0,0,iDULL11);}
    else if (name.compare("iDULL12") == 0){
        iDULL12 = value;
        Du_LL.assignim(0,1,iDULL12);}
    else if (name.compare("iDULL13") == 0){
        iDULL13 = value;
        Du_LL.assignim(0,2,iDULL13);}
    else if (name.compare("iDULL21") == 0){
        iDULL21 = value;
        Du_LL.assignim(1,0,iDULL21);}
    else if (name.compare("iDULL22") == 0){
        iDULL22 = value;
        Du_LL.assignim(1,1,iDULL22);}
    else if (name.compare("iDULL23") == 0){
        iDULL23 = value;
        Du_LL.assignim(1,2,iDULL23);}
    else if (name.compare("iDULL31") == 0){
        iDULL31 = value;
        Du_LL.assignim(2,0,iDULL31);}
    else if (name.compare("iDULL32") == 0){
        iDULL32 = value;
        Du_LL.assignim(2,1,iDULL32);}
    else if (name.compare("iDULL33") == 0){
        iDULL33 = value;
        Du_LL.assignim(2,2,iDULL33);}
    
   
    else if (name.compare("rDURL11") == 0){
        rDURL11 = value;
        Du_RL.assignre(0,0,rDURL11);}
    else if (name.compare("rDURL12") == 0){
        rDURL12 = value;
        Du_RL.assignre(0,1,rDURL12);}
    else if (name.compare("rDURL13") == 0){
        rDURL13 = value;
        Du_RL.assignre(0,2,rDURL13);}
    else if (name.compare("rDURL21") == 0){
        rDURL21 = value;
        Du_RL.assignre(1,0,rDURL21);}
    else if (name.compare("rDURL22") == 0){
        rDURL22 = value;
        Du_RL.assignre(1,1,rDURL22);}
    else if (name.compare("rDURL23") == 0){
        rDURL23 = value;
        Du_RL.assignre(1,2,rDURL23);}
    else if (name.compare("rDURL31") == 0){
        rDURL31 = value;
        Du_RL.assignre(2,0,rDURL31);}
    else if (name.compare("rDURL32") == 0){
        rDURL32 = value;
        Du_RL.assignre(2,1,rDURL32);}
    else if (name.compare("rDURL33") == 0){
        rDURL33 = value;
        Du_RL.assignre(2,2,rDURL33);}
    
    else if (name.compare("iDURL11") == 0){
        iDURL11 = value;
        Du_RL.assignim(0,0,iDURL11);}
    else if (name.compare("iDURL12") == 0){
        iDURL12 = value;
        Du_RL.assignim(0,1,iDURL12);}
    else if (name.compare("iDURL13") == 0){
        iDURL13 = value;
        Du_RL.assignim(0,2,iDURL13);}
    else if (name.compare("iDURL21") == 0){
        iDURL21 = value;
        Du_RL.assignim(1,0,iDURL21);}
    else if (name.compare("iDURL22") == 0){
        iDURL22 = value;
        Du_RL.assignim(1,1,iDURL22);}
    else if (name.compare("iDURL23") == 0){
        iDURL23 = value;
        Du_RL.assignim(1,2,iDURL23);}
    else if (name.compare("iDURL31") == 0){
        iDURL31 = value;
        Du_RL.assignim(2,0,iDURL31);}
    else if (name.compare("iDURL32") == 0){
        iDURL32 = value;
        Du_RL.assignim(2,1,iDURL32);}
    else if (name.compare("iDURL33") == 0){
        iDURL33 = value;
        Du_RL.assignim(2,2,iDURL33);}
    
    
    else if (name.compare("rDULR11") == 0){
        rDULR11 = value;
        Du_LR.assignre(0,0,rDULR11);}
    else if (name.compare("rDULR12") == 0){
        rDULR12 = value;
        Du_LR.assignre(0,1,rDULR12);}
    else if (name.compare("rDULR13") == 0){
        rDULR13 = value;
        Du_LR.assignre(0,2,rDULR13);}
    else if (name.compare("rDULR21") == 0){
        rDULR21 = value;
        Du_LR.assignre(1,0,rDULR21);}
    else if (name.compare("rDULR22") == 0){
        rDULR22 = value;
        Du_LR.assignre(1,1,rDULR22);}
    else if (name.compare("rDULR23") == 0){
        rDULR23 = value;
        Du_LR.assignre(1,2,rDULR23);}
    else if (name.compare("rDULR31") == 0){
        rDULR31 = value;
        Du_LR.assignre(2,0,rDULR31);}
    else if (name.compare("rDULR32") == 0){
        rDULR32 = value;
        Du_LR.assignre(2,0,rDULR31);}
    else if (name.compare("rDULR33") == 0){
        rDULR33 = value;
        Du_LR.assignre(2,2,rDULR33);}
    
    else if (name.compare("iDULR11") == 0){
        iDULR11 = value;
        Du_LR.assignim(0,0,iDULR11);}
    else if (name.compare("iDULR12") == 0){
        iDULR12 = value;
        Du_LR.assignim(0,1,iDULR12);}
    else if (name.compare("iDULR13") == 0){
        iDULR13 = value;
        Du_LR.assignim(0,2,iDULR13);}
    else if (name.compare("iDULR21") == 0){
        iDULR21 = value;
        Du_LR.assignim(1,0,iDULR21);}
    else if (name.compare("iDULR22") == 0){
        iDULR22 = value;
        Du_LR.assignim(1,1,iDULR22);}
    else if (name.compare("iDULR23") == 0){
        iDULR23 = value;
        Du_LR.assignim(1,2,iDULR23);}
    else if (name.compare("iDULR31") == 0){
        iDULR31 = value;
        Du_LR.assignim(2,0,iDULR31);}
    else if (name.compare("iDULR32") == 0){
        iDULR32 = value;
        Du_LR.assignim(2,1,iDULR32);}
    else if (name.compare("iDULR33") == 0){
        iDULR33 = value;
        Du_LR.assignim(2,2,iDULR33);}
    
    
    else if (name.compare("rDURR11") == 0){
        rDURR11 = value;
        Du_RR.assignre(0,0,rDURR11);}
    else if (name.compare("rDURR12") == 0){
        rDURR12 = value;
        Du_RR.assignre(0,1,rDURR12);}
    else if (name.compare("rDURR13") == 0){
        rDURR13 = value;
        Du_RR.assignre(0,2,rDURR13);}
    else if (name.compare("rDURR21") == 0){
        rDURR21 = value;
        Du_RR.assignre(1,0,rDURR21);}
    else if (name.compare("rDURR22") == 0){
        rDURR22 = value;
        Du_RR.assignre(1,1,rDURR22);}
    else if (name.compare("rDURR23") == 0){
        rDURR23 = value;
        Du_RR.assignre(1,2,rDURR23);}
    else if (name.compare("rDURR31") == 0){
        rDURR31 = value;
        Du_RR.assignre(2,0,rDURR31);}
    else if (name.compare("rDURR32") == 0){
        rDURR32 = value;
        Du_RR.assignre(2,1,rDURR32);}
    else if (name.compare("rDURR33") == 0){
        rDURR33 = value;
        Du_RR.assignre(2,2,rDURR33);}
    
    else if (name.compare("iDURR11") == 0){
        iDURR11 = value;
        Du_RR.assignim(0,0,iDURR11);}
    else if (name.compare("iDURR12") == 0){
        iDURR12 = value;
        Du_RR.assignim(0,1,iDURR12);}
    else if (name.compare("iDURR13") == 0){
        iDURR13 = value;
        Du_RR.assignim(0,2,iDURR13);}
    else if (name.compare("iDURR21") == 0){
        iDURR21 = value;
        Du_RR.assignim(1,0,iDURR21);}
    else if (name.compare("iDURR22") == 0){
        iDURR22 = value;
        Du_RR.assignim(1,1,iDURR22);}
    else if (name.compare("iDURR23") == 0){
        iDURR23 = value;
        Du_RR.assignim(1,2,iDURR23);}
    else if (name.compare("iDURR31") == 0){
        iDURR31 = value;
        Du_RR.assignim(2,0,iDURR31);}
    else if (name.compare("iDURR32") == 0){
        iDURR32 = value;
        Du_RR.assignim(2,1,iDURR32);}
    else if (name.compare("iDURR33") == 0){
        iDURR33 = value;
        Du_RR.assignim(2,2,iDURR33);}
    
    // FILLS THE DELTAS, DOWN SECTOR
    else if (name.compare("rDDLL11") == 0){
        rDDLL11 = value;
        Dd_LL.assignre(0,0,rDDLL11);}
    else if (name.compare("rDDLL12") == 0){
        rDDLL12 = value;
        Dd_LL.assignre(0,1,rDDLL12);}
    else if (name.compare("rDDLL13") == 0){
        rDDLL13 = value;
        Dd_LL.assignre(0,2,rDDLL13);}
    else if (name.compare("rDDLL21") == 0){
        rDDLL21 = value;
        Dd_LL.assignre(1,0,rDDLL21);}
    else if (name.compare("rDDLL22") == 0){
        rDDLL22 = value;
        Dd_LL.assignre(1,1,rDDLL22);}
    else if (name.compare("rDDLL23") == 0){
        rDDLL23 = value;
        Dd_LL.assignre(1,2,rDDLL23);}
    else if (name.compare("rDDLL31") == 0){
        rDDLL31 = value;
        Dd_LL.assignre(2,0,rDDLL31);}
    else if (name.compare("rDDLL32") == 0){
        rDDLL32 = value;
        Dd_LL.assignre(2,1,rDDLL32);}
    else if (name.compare("rDDLL33") == 0){
        rDDLL33 = value;
        Dd_LL.assignre(2,2,rDDLL33);}
    
    else if (name.compare("iDDLL11") == 0){
        iDDLL11 = value;
        Dd_LL.assignim(0,0,iDDLL11);}
    else if (name.compare("iDDLL12") == 0){
        iDDLL12 = value;
        Dd_LL.assignim(0,1,iDDLL12);}
    else if (name.compare("iDDLL13") == 0){
        iDDLL13 = value;
        Dd_LL.assignim(0,2,iDDLL13);}
    else if (name.compare("iDDLL21") == 0){
        iDDLL21 = value;
        Dd_LL.assignim(1,0,iDDLL21);}
    else if (name.compare("iDDLL22") == 0){
        iDDLL22 = value;
        Dd_LL.assignim(1,1,iDDLL22);}
    else if (name.compare("iDDLL23") == 0){
        iDDLL23 = value;
        Dd_LL.assignim(1,2,iDDLL23);}
    else if (name.compare("iDDLL31") == 0){
        iDDLL31 = value;
        Dd_LL.assignim(2,0,iDDLL31);}
    else if (name.compare("iDDLL32") == 0){
        iDDLL32 = value;
        Dd_LL.assignim(2,1,iDDLL32);}
    else if (name.compare("iDDLL33") == 0){
        iDDLL33 = value;
        Dd_LL.assignim(2,2,iDDLL33);}
    
   
    else if (name.compare("rDDRL11") == 0){
        rDDRL11 = value;
        Dd_RL.assignre(0,0,rDDRL11);}
    else if (name.compare("rDDRL12") == 0){
        rDDRL12 = value;
        Dd_RL.assignre(0,1,rDDRL12);}
    else if (name.compare("rDDRL13") == 0){
        rDDRL13 = value;
        Dd_RL.assignre(0,2,rDDRL13);}
    else if (name.compare("rDDRL21") == 0){
        rDDRL21 = value;
        Dd_RL.assignre(1,0,rDDRL21);}
    else if (name.compare("rDDRL22") == 0){
        rDDRL22 = value;
        Dd_RL.assignre(1,1,rDDRL22);}
    else if (name.compare("rDDRL23") == 0){
        rDDRL23 = value;
        Dd_RL.assignre(1,2,rDDRL23);}
    else if (name.compare("rDDRL31") == 0){
        rDDRL31 = value;
        Dd_RL.assignre(2,0,rDDRL31);}
    else if (name.compare("rDDRL32") == 0){
        rDDRL32 = value;
        Dd_RL.assignre(2,1,rDDRL32);}
    else if (name.compare("rDDRL33") == 0){
        rDDRL33 = value;
        Dd_RL.assignre(2,2,rDDRL33);}
    
    else if (name.compare("iDDRL11") == 0){
        iDDRL11 = value;
        Dd_RL.assignim(0,0,iDDRL11);}
    else if (name.compare("iDDRL12") == 0){
        iDDRL12 = value;
        Dd_RL.assignim(0,1,iDDRL12);}
    else if (name.compare("iDDRL13") == 0){
        iDDRL13 = value;
        Dd_RL.assignim(0,2,iDDRL13);}
    else if (name.compare("iDDRL21") == 0){
        iDDRL21 = value;
        Dd_RL.assignim(1,0,iDDRL21);}
    else if (name.compare("iDDRL22") == 0){
        iDDRL22 = value;
        Dd_RL.assignim(1,1,iDDRL22);}
    else if (name.compare("iDDRL23") == 0){
        iDDRL23 = value;
        Dd_RL.assignim(1,2,iDDRL23);}
    else if (name.compare("iDDRL31") == 0){
        iDDRL31 = value;
        Dd_RL.assignim(2,0,iDDRL31);}
    else if (name.compare("iDDRL32") == 0){
        iDDRL32 = value;
        Dd_RL.assignim(2,1,iDDRL32);}
    else if (name.compare("iDDRL33") == 0){
        iDDRL33 = value;
        Dd_RL.assignim(2,2,iDDRL33);}
    
    
    else if (name.compare("rDDLR11") == 0){
        rDDLR11 = value;
        Dd_LR.assignre(0,0,rDDLR11);}
    else if (name.compare("rDDLR12") == 0){
        rDDLR12 = value;
        Dd_LR.assignre(0,1,rDDLR12);}
    else if (name.compare("rDDLR13") == 0){
        rDDLR13 = value;
        Dd_LR.assignre(0,2,rDDLR13);}
    else if (name.compare("rDDLR21") == 0){
        rDDLR21 = value;
        Dd_LR.assignre(1,0,rDDLR21);}
    else if (name.compare("rDDLR22") == 0){
        rDDLR22 = value;
        Dd_LR.assignre(1,1,rDDLR22);}
    else if (name.compare("rDDLR23") == 0){
        rDDLR23 = value;
        Dd_LR.assignre(1,2,rDDLR23);}
    else if (name.compare("rDDLR31") == 0){
        rDDLR31 = value;
        Dd_LR.assignre(2,0,rDDLR31);}
    else if (name.compare("rDDLR32") == 0){
        rDDLR32 = value;
        Dd_LR.assignre(2,0,rDDLR31);}
    else if (name.compare("rDDLR33") == 0){
        rDDLR33 = value;
        Dd_LR.assignre(2,2,rDDLR33);}
    
    else if (name.compare("iDDLR11") == 0){
        iDDLR11 = value;
        Dd_LR.assignim(0,0,iDDLR11);}
    else if (name.compare("iDDLR12") == 0){
        iDDLR12 = value;
        Dd_LR.assignim(0,1,iDDLR12);}
    else if (name.compare("iDDLR13") == 0){
        iDDLR13 = value;
        Dd_LR.assignim(0,2,iDDLR13);}
    else if (name.compare("iDDLR21") == 0){
        iDDLR21 = value;
        Dd_LR.assignim(1,0,iDDLR21);}
    else if (name.compare("iDDLR22") == 0){
        iDDLR22 = value;
        Dd_LR.assignim(1,1,iDDLR22);}
    else if (name.compare("iDDLR23") == 0){
        iDDLR23 = value;
        Dd_LR.assignim(1,2,iDDLR23);}
    else if (name.compare("iDDLR31") == 0){
        iDDLR31 = value;
        Dd_LR.assignim(2,0,iDDLR31);}
    else if (name.compare("iDDLR32") == 0){
        iDDLR32 = value;
        Dd_LR.assignim(2,0,iDDLR31);}
    else if (name.compare("iDDLR33") == 0){
        iDDLR33 = value;
        Dd_LR.assignim(2,2,iDDLR33);}
    
    
    else if (name.compare("rDDRR11") == 0){
        rDDRR11 = value;
        Dd_RR.assignre(0,0,rDDRR11);}
    else if (name.compare("rDDRR12") == 0){
        rDDRR12 = value;
        Dd_RR.assignre(0,1,rDDRR12);}
    else if (name.compare("rDDRR13") == 0){
        rDDRR13 = value;
        Dd_RR.assignre(0,2,rDDRR13);}
    else if (name.compare("rDDRR21") == 0){
        rDDRR21 = value;
        Dd_RR.assignre(1,0,rDDRR21);}
    else if (name.compare("rDDRR22") == 0){
        rDDRR22 = value;
        Dd_RR.assignre(1,1,rDDRR22);}
    else if (name.compare("rDDRR23") == 0){
        rDDRR23 = value;
        Dd_RR.assignre(1,2,rDDRR23);}
    else if (name.compare("rDDRR31") == 0){
        rDDRR31 = value;
        Dd_RR.assignre(2,0,rDDRR31);}
    else if (name.compare("rDDRR32") == 0){
        rDDRR32 = value;
        Dd_RR.assignre(2,1,rDDRR32);}
    else if (name.compare("rDDRR33") == 0){
        rDDRR33 = value;
        Dd_RR.assignre(2,2,rDDRR33);}
    
    else if (name.compare("iDDRR11") == 0){
        iDDRR11 = value;
        Dd_RR.assignim(0,0,iDDRR11);}
    else if (name.compare("iDDRR12") == 0){
        iDDRR12 = value;
        Dd_RR.assignim(0,1,iDDRR12);}
    else if (name.compare("iDDRR13") == 0){
        iDDRR13 = value;
        Dd_RR.assignim(0,2,iDDRR13);}
    else if (name.compare("iDDRR21") == 0){
        iDDRR21 = value;
        Dd_RR.assignim(1,0,iDDRR21);}
    else if (name.compare("iDDRR22") == 0){
        iDDRR22 = value;
        Dd_RR.assignim(1,1,iDDRR22);}
    else if (name.compare("iDDRR23") == 0){
        iDDRR23 = value;
        Dd_RR.assignim(1,2,iDDRR23);}
    else if (name.compare("iDDRR31") == 0){
        iDDRR31 = value;
        Dd_RR.assignim(2,0,iDDRR31);}
    else if (name.compare("iDDRR32") == 0){
        iDDRR32 = value;
        Dd_RR.assignim(2,1,iDDRR32);}
    else if (name.compare("iDDRR33") == 0){
        iDDRR33 = value;
        Dd_RR.assignim(2,2,iDDRR33);}
    
    else 
        StandardModel::setParameter(name, value);
}

bool SUSYMassInsertion::CheckParameters(const std::map<std::string, double>& DPars)
{
    for(int i=0; i<NSusyMIvars; i++)
        if(DPars.find(SusyMIvars[i])==DPars.end()) {
            std::cout << "missing mandatory SUSY MIA parameter " << SusyMIvars[i] << std::endl;
            return false;
        }
    return(StandardModel::CheckParameters(DPars));
}
