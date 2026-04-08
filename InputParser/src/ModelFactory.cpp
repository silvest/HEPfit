/* 
 * Copyright (C) 2014 HEPfit Collaboration
 *
 *
 * For the licensing terms see doc/COPYING.
 */

#include "ModelFactory.h"
#include "NPSTU.h"
#include "NPSTUZbbbarLR.h"
#include "NPEpsilons.h"
#include "NPEpsilons_pureNP.h"
#include "NPHiggs.h"
#include "NPZbbbar.h"
#include "NPZbbbarLinearized.h"
#include "NPSMEFTd6.h"
#include "NPSMEFTd6General.h"
#include "NPSMEFTd6U2qU1le.h"
#include "NPSMEFTd6U2.h"
#include "NPSMEFTd6U3.h"
#include "NPSMEFTd6MFV.h"

#include "NPd6SILH.h"

#include "NPSMEFT6dtopquark.h"
#include "SigmaBR.h"
#include "HiggsKigen.h"
#include "HiggsChiral.h"
#include "GeneralSUSY.h"
#include "SUSYMassInsertion.h"
#include "SUSY.h"
#include "THDM.h"
#include "FlavourWilsonCoefficient.h"
#include "FlavourWilsonCoefficient_DF2.h"
#include "RealWeakEFTLFV.h"
#include "RealWeakEFTCC.h"
#include "LoopMediators.h"
#include "NPDF2.h"
#include "CMFV.h"
#include "GeorgiMachacek.h"
#include "LeftRightSymmetricModel.h"
/* BEGIN: REMOVE FROM THE PACKAGE */
#include "pMSSM.h"
#include "MFV.h"
#include "NPSTUVWXY.h"
#include "GeneralTHDM.h"
#include "GeneralTHDMZ2.h"
#include "THDMW.h"
/* END: REMOVE FROM THE PACKAGE */
        
ModelFactory::ModelFactory()
{
    modelFactory["StandardModel"] = []() -> StandardModel* { return new StandardModel(); };
    modelFactory["NPSTU"] = []() -> StandardModel* { return new NPSTU(); };
    modelFactory["NPSTUZbbbarLR"] = []() -> StandardModel* { return new NPSTUZbbbarLR(); };
    modelFactory["NPEpsilons"] = []() -> StandardModel* { return new NPEpsilons(); };
    modelFactory["NPEpsilons_pureNP"] = []() -> StandardModel* { return new NPEpsilons_pureNP(); };
    modelFactory["NPHiggs"] = []() -> StandardModel* { return new NPHiggs(); };
    modelFactory["NPZbbbar"] = []() -> StandardModel* { return new NPZbbbar(false); };
    modelFactory["NPZbbbarLR"] = []() -> StandardModel* { return new NPZbbbar(true); };
    modelFactory["NPZbbbarLinearized"] = []() -> StandardModel* { return new NPZbbbarLinearized(false); };
    modelFactory["NPZbbbarLinearizedLR"] = []() -> StandardModel* { return new NPZbbbarLinearized(true); };
    modelFactory["NPSMEFTd6"] = []() -> StandardModel* { return new NPSMEFTd6(false, false); };
    modelFactory["NPSMEFTd6General"] = []() -> StandardModel* { return new NPSMEFTd6General(); };
    modelFactory["NPSMEFTd6U2qU1le"] = []() -> StandardModel* { return new NPSMEFTd6U2qU1le(); };
    modelFactory["NPSMEFTd6U2"] = []() -> StandardModel* { return new NPSMEFTd6U2(); };
    modelFactory["NPSMEFTd6U3"] = []() -> StandardModel* { return new NPSMEFTd6U3(); };
    modelFactory["NPSMEFTd6MFV"] = []() -> StandardModel* { return new NPSMEFTd6MFV(); };
    
    modelFactory["NPd6SILH"] = []() -> StandardModel* { return new NPd6SILH(); };
    
    modelFactory["NPSMEFT6dtopquark"] = []() -> StandardModel* { return new NPSMEFT6dtopquark(); };

    modelFactory["SigmaBR"] = []() -> StandardModel* { return new SigmaBR(); };
    modelFactory["NPSMEFTd6_LFU_QFU"] = []() -> StandardModel* { return new NPSMEFTd6(true, true); };
    modelFactory["HiggsKigen"] = []() -> StandardModel* { return new HiggsKigen(); };
    modelFactory["HiggsChiral"] = []() -> StandardModel* { return new HiggsChiral(); };
    modelFactory["NPDF2"] = []() -> StandardModel* { return new NPDF2(); };
    modelFactory["CMFV"] = []() -> StandardModel* { return new CMFV(); };
    modelFactory["FlavourWilsonCoefficient"] = []() -> StandardModel* { return new FlavourWilsonCoefficient(); };
    modelFactory["FlavourWilsonCoefficient_DF2"] = []() -> StandardModel* { return new FlavourWilsonCoefficient_DF2(); };
    modelFactory["RealWeakEFTLFV"] = []() -> StandardModel* { return new RealWeakEFTLFV(); };
    modelFactory["RealWeakEFTCC"] = []() -> StandardModel* { return new RealWeakEFTCC(0); };
    modelFactory["RealWeakEFTCCPM"] = []() -> StandardModel* { return new RealWeakEFTCC(1); };
    modelFactory["LoopMediators"] = []() -> StandardModel* { return new LoopMediators(); };
    modelFactory["SUSYMassInsertion"] = []() -> StandardModel* { return new SUSYMassInsertion(); };
    modelFactory["THDM"] = []() -> StandardModel* { return new THDM(); };
    modelFactory["GeneralSUSY"] = []() -> StandardModel* { return new GeneralSUSY(); };
    modelFactory["GeorgiMachacek"] = []() -> StandardModel* { return new GeorgiMachacek(); };
    modelFactory["LeftRightSymmetricModel"] = []() -> StandardModel* { return new LeftRightSymmetricModel(); };
/* BEGIN: REMOVE FROM THE PACKAGE */
    modelFactory["MFV"] = []() -> StandardModel* { return new MFV(); };
    modelFactory["pMSSM"] = []() -> StandardModel* { return new pMSSM(); };
    modelFactory["NPSTUVWXY"] = []() -> StandardModel* { return new NPSTUVWXY(); };
    modelFactory["GeneralTHDM"] = []() -> StandardModel* { return new GeneralTHDM(); };
    modelFactory["GeneralTHDMZ2"] = []() -> StandardModel* { return new GeneralTHDMZ2(); };
    modelFactory["THDMW"] = []() -> StandardModel* { return new THDMW(); };
/* END: REMOVE FROM THE PACKAGE */
}

void ModelFactory::addModelToFactory(const std::string name, std::function<StandardModel*() > funct)
{
    modelFactory[name] = funct;
}

StandardModel* ModelFactory::CreateModel(const std::string& name)
{
    if (modelFactory.find(name) == modelFactory.end())
        throw std::runtime_error("ERROR: Wrong model " + name + " passed to ModelFactory.\n");
    return (modelFactory[name]());
}
