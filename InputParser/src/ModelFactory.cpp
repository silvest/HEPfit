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
#include "THDMW.h"
/* END: REMOVE FROM THE PACKAGE */
#include <boost/bind/bind.hpp>
using namespace boost::placeholders;
        
ModelFactory::ModelFactory()
{
    modelFactory["StandardModel"] = boost::factory<StandardModel*>();
    modelFactory["NPSTU"] = boost::factory<NPSTU*>();
    modelFactory["NPSTUZbbbarLR"] = boost::factory<NPSTUZbbbarLR*>();
    modelFactory["NPEpsilons"] = boost::factory<NPEpsilons*>();
    modelFactory["NPEpsilons_pureNP"] = boost::factory<NPEpsilons_pureNP*>();
    modelFactory["NPHiggs"] = boost::factory<NPHiggs*>();
    modelFactory["NPZbbbar"] = bind(boost::factory<NPZbbbar*>(), false);
    modelFactory["NPZbbbarLR"] = bind(boost::factory<NPZbbbar*>(), true);
    modelFactory["NPZbbbarLinearized"] = bind(boost::factory<NPZbbbarLinearized*>(), false);
    modelFactory["NPZbbbarLinearizedLR"] = bind(boost::factory<NPZbbbarLinearized*>(), true);
    modelFactory["NPSMEFTd6"] = bind(boost::factory<NPSMEFTd6*>(), false, false);
    
    modelFactory["NPSMEFT6dtopquark"] = boost::factory<NPSMEFT6dtopquark*>();

    //modelFactory["NPSMEFT6dtopquark"] = bind(boost::factory<NPSMEFT6dtopquark*>(), false, false);
    modelFactory["SigmaBR"] = boost::factory<SigmaBR*>();
    modelFactory["NPSMEFTd6_LFU_QFU"] = bind(boost::factory<NPSMEFTd6*>(), true, true);
    modelFactory["HiggsKigen"] = boost::factory<HiggsKigen*>();
    modelFactory["HiggsChiral"] = boost::factory<HiggsChiral*>();
    modelFactory["NPDF2"] = boost::factory<NPDF2*>();
    modelFactory["CMFV"] = boost::factory<CMFV*>();
    modelFactory["FlavourWilsonCoefficient"] = boost::factory<FlavourWilsonCoefficient*>();
    modelFactory["FlavourWilsonCoefficient_DF2"] = boost::factory<FlavourWilsonCoefficient_DF2*>();
    modelFactory["RealWeakEFTLFV"] = boost::factory<RealWeakEFTLFV*>();
    modelFactory["RealWeakEFTCC"] = bind(boost::factory<RealWeakEFTCC*>(), 0);
    modelFactory["RealWeakEFTCCPM"] = bind(boost::factory<RealWeakEFTCC*>(), 1);
    modelFactory["LoopMediators"] = boost::factory<LoopMediators*>();
    modelFactory["SUSYMassInsertion"] = boost::factory<SUSYMassInsertion*>();
    modelFactory["THDM"] = boost::factory<THDM*>();
    modelFactory["GeneralSUSY"] = boost::factory<GeneralSUSY*>();
    modelFactory["GeorgiMachacek"] = boost::factory<GeorgiMachacek*>();
    modelFactory["LeftRightSymmetricModel"] = boost::factory<LeftRightSymmetricModel*>();
/* BEGIN: REMOVE FROM THE PACKAGE */
    modelFactory["MFV"] = boost::factory<MFV*>();
    modelFactory["pMSSM"] = boost::factory<pMSSM*>();
    modelFactory["NPSTUVWXY"] = boost::factory<NPSTUVWXY*>();
    modelFactory["GeneralTHDM"] = boost::factory<GeneralTHDM*>();
    modelFactory["THDMW"] = boost::factory<THDMW*>();
/* END: REMOVE FROM THE PACKAGE */
}

void ModelFactory::addModelToFactory(const std::string name, boost::function<StandardModel*() > funct)
{
    modelFactory[name] = funct;
}

StandardModel* ModelFactory::CreateModel(const std::string& name)
{
    if (modelFactory.find(name) == modelFactory.end())
        throw std::runtime_error("ERROR: Wrong model " + name + " passed to ModelFactory.\n");
    return (modelFactory[name]());
}
